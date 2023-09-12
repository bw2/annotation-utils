import gzip
import pandas as pd
import requests

from annotations.cache_utils import cache_data_table

URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.ensembl_genomic.gtf.gz"


@cache_data_table
def get_MANE_ensembl_transcript_table(feature_type=None):
    """Download the MANE ensembl gtf and return it as a pandas DataFrame

    Args:
        feature_type (str): if not None, only keep features of this type. Allowed values are:
            'gene', 'transcript', 'CDS', 'UTR', 'exon', 'start_codon', 'stop_codon'

    Return:
        pandas.DataFrame: MANE ensembl gtf as a pandas DataFrame
    """

    r = requests.get(URL)
    if not r.ok:
        raise Exception(f"Failed to download {URL}: {r}")

    table_rows = []
    for line in gzip.decompress(r.content).decode("UTF-8").strip().split("\n"):
        fields = line.split("\t")

        if len(fields) < 2:
            print(f"WARNING: unable to parse line: {line}")
            continue

        if feature_type is not None and fields[2] != feature_type:
            continue

        record = {
            "chrom": fields[0],
            "source": fields[1],
            "feature": fields[2],
            "start": int(fields[3]),
            "end": int(fields[4]),
            "strand": fields[6],
        }

        info = fields[8]
        record.update({k: v.strip(';" ') for k, v in [x.split(" ") for x in info.split("; ")]})

        record["is_MANE_select"] = record.get("tag") == "MANE_Select"
        table_rows.append(record)

    return pd.DataFrame(table_rows)
