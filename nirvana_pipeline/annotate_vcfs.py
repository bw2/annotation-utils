import hailtop.fs as hfs
import logging
import os
import re

from step_pipeline import pipeline, Backend, Localize, Delocalize, files_exist

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO)
logger = logging.getLogger(__name__)

DOCKER_IMAGE = "annotation/nirvana:3.14"

#REFERENCE_FASTA_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
#REFERENCE_FASTA_FAI_PATH = "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai"

NIRVANA_REF_DATA_BUCKET = "gs://ttn-neb-analysis"

def main():

    bp = pipeline("run Nirvana", backend=Backend.HAIL_BATCH_SERVICE, config_file_path="~/.step_pipeline")
    parser = bp.get_config_arg_parser()
    #parser.add_argument("--reference-fasta", default=REFERENCE_FASTA_PATH)
    #parser.add_argument("--reference-fasta-fai", default=REFERENCE_FASTA_FAI_PATH)
    parser.add_argument("--output-dir", help="Default is the same directory as the input vcf")
    parser.add_argument("-n", type=int, help="Only process the first n inputs. Useful for testing.")
    parser.add_argument("vcf_list", action="append",
                        help="Either a list of gs:// paths, or a text file containing a list of gs:// paths")
    args = bp.parse_known_args()



    if len(args.vcf_list) == 1 and os.path.isfile(args.vcf_list[0]):
        with open(args.vcf_list[0]) as f:
            args.vcf_list = f.read().strip().split("\n")

    if args.n:
        args.vcf_list = args.vcf_list[:args.n]

    for i, vcf_path in enumerate(args.vcf_list):
        if not vcf_path.startswith("gs://"):
            parser.error(f"VCF path #{i+1} doesn't start with gs://: '{vcf_path}'")
        if not hfs.exists(vcf_path):
            parser.error(f"VCF file doesn't exist: {vcf_path}")

    for vcf_path in args.vcf_list:
        print("Annotating ", vcf_path)

        filename_prefix = re.sub(".vcf(.gz)?$", "", os.path.basename(vcf_path))
        output_dir = args.output_dir or os.path.dirname(vcf_path)

        num_cpu = 2
        s1 = bp.new_step(
            f"run Nirvana: {filename_prefix}",
            arg_suffix=f"dv",
            image=DOCKER_IMAGE,
            step_number=1,
            cpu=num_cpu,
            storage="100Gi",
            #localize_by=Localize.HAIL_BATCH_CLOUDFUSE,
            localize_by=Localize.COPY,
            delocalize_by=Delocalize.COPY,
            output_dir=output_dir)

        # set job inputs & outputs
        nirvana_bucket = s1.input(NIRVANA_REF_DATA_BUCKET, localize_by=Localize.HAIL_BATCH_CLOUDFUSE)
        local_input_vcf = s1.input(vcf_path)

        s1.command("set -ex")
        s1.command(f"""dotnet /opt/nirvana/Nirvana.dll \
            -c {nirvana_bucket}/Nirvana/Data/Cache/GRCh38/Both \
            -r {nirvana_bucket}/Nirvana/Data/References/Homo_sapiens.GRCh38.Nirvana.dat \
            --sd {nirvana_bucket}/Nirvana/Data/SupplementaryAnnotation/GRCh38 \
            -i {local_input_vcf} \
            -o {filename_prefix}
        """)

        s1.output(f"{filename_prefix}.json.gz")
        s1.output(f"{filename_prefix}.json.gz.jsi")
        s1.command(f"date")

    bp.run()


if __name__ == "__main__":
    main()
