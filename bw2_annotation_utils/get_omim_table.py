import base64
from cache_utils import cache_data_table
import datetime
import json
import pandas as pd
import requests

INPUT_TABLE_HEADER = [
    'mim_number',                   # 0
    'phenotype_mim_number',         # 1
    'phenotypic_series_number',     # 2
    'phenotype_inheritance',        # 3

    'locus_size',                   # 4
    'locus',                        # 5
    'locus_hg19',                   # 6

    'cyto',                         # 7
    'gene_symbols',                 # 8
    'gene_id',                      # 9
    'gene_description',             # 10
    'phenotype_description',        # 11
    'date_created',                 # 12
    'date_updated',                 # 13
    'mouse_gene_id',                # 14

    'oe_lof_upper',                 # 15
    'pLI',                          # 16
    'mis_z',                        # 17

    'text',                         # 15
    'comments',                     # 16

    'xstart',
    'xend',
    'xstart_hg19',
    'xend_hg19',
    'phenotype_map_method',
    #'liftover_to_hg19_failed',
]

OUTPUT_COLUMNS = [
    'chrom',
    'start',
    'end',
    'mim_number',
    'phenotype_mim_number',
    'phenotypic_series_number',
    'phenotype_inheritance',
    'gene_symbols',
    'gene_id',
    'gene_description',
    'phenotype_description',
    'date_created',
    'date_updated',
    'mouse_gene_id',
    'oe_lof_upper',
    'pLI',
    'mis_z',
    'text',
    'comments',
]
MAX_GENE_SIZE = 5*10**6  # 5 Mbases (dystrophin is 2.3Mb)

"""
Example row:
$1                      chrom : 1
$2                      start : 7784284
$3                        end : 7845180
$4                 mim_number : 603427
$5       phenotype_mim_number : 616882
$6   phenotypic_series_number :
$7      phenotype_inheritance : Autosomal dominant
$8               gene_symbols : FASPS3, PER3
$9                    gene_id : ENSG00000049246
$10          gene_description : Period circadian regulator 3
$11     phenotype_description : ?Advanced sleep phase syndrome, familial, 3
$12              date_created :
$13              date_updated :
$14             mouse_gene_id : Per3 (MGI:1277134)
$15              oe_lof_upper : 8.4100e-01
$16                       pLI : 1.3661e-17
$17                     mis_z : -1.0855e-01
$18                      text :
$19                  comments : mutation identified in 1 FASPS3 family
"""


@cache_data_table
def get_omim_table():
    """Retrieves the latest OMIM table"""

    r = requests.get("https://broadinstitute.github.io/omim-search-p/d")
    if not r.ok:
        raise Exception(f"Failed to download latest OMIM json from omim-search-p: {r}")

    json_string = r.json()
    json_obj = json.loads(base64.b64decode(json_string[0]))
    omim_df = pd.DataFrame(columns=INPUT_TABLE_HEADER, data=json_obj["data"])

    #print(omim_df.columns)
    #print(omim_df["locus"])

    omim_df[["chrom", "interval"]] = omim_df["locus"].str.split(":", expand=True)
    omim_df[["start", "end"]] = omim_df["interval"].str.split("-", expand=True).astype("int32")

    omim_df = omim_df[omim_df["locus_size"] < MAX_GENE_SIZE]
    omim_df = omim_df[omim_df["start"] > 1]

    omim_df["pLI"] = omim_df["pLI"].replace("", float('nan')).replace("NA", float('nan')).astype(float)
    omim_df["mis_z"] = omim_df["mis_z"].replace("", float('nan')).replace("NA", float('nan')).astype(float)
    omim_df["oe_lof_upper"] = omim_df["oe_lof_upper"].replace("", float('nan')).replace("NA", float('nan')).astype(float)

    """
    Autosomal recessive                                          3333
    Autosomal dominant                                           2480
    X-linked recessive                                            217
    Autosomal dominant, Autosomal recessive                       209
    X-linked                                                       75
    X-linked dominant                                              69
    Somatic mutation, Autosomal dominant                           60
    Multifactorial                                                 31
    Multifactorial, Autosomal dominant, Autosomal recessive        12
    Digenic recessive, Autosomal recessive                         11
    Somatic mutation                                               10
    Autosomal dominant, Multifactorial                              9
    Digenic recessive                                               8
    Isolated cases, Autosomal dominant                              7
    Digenic dominant, Autosomal recessive                           6
    X-linked dominant, X-linked recessive                           5
    Isolated cases                                                  5
    Digenic dominant                                                5
    ?Autosomal dominant                                             4
    Digenic dominant, Autosomal dominant, Autosomal recessive       4
    Somatic mutation, Autosomal dominant, Autosomal recessive       4
    Digenic dominant, Autosomal dominant                            3
    X-linked, Isolated cases, Multifactorial                        3
    Y-linked                                                        3
    Multifactorial, Autosomal recessive                             2
    Somatic mosaicism, Autosomal recessive                          2
    Pseudoautosomal dominant                                        2
    Pseudoautosomal recessive                                       2
    Isolated cases, Somatic mutation                                1
    Mitochondrial                                                   1
    ?Autosomal dominant, Autosomal recessive                        1
    X-linked, Somatic mosaicism                                     1
    Somatic mosaicism, Autosomal dominant                           1
    X-linked dominant, Somatic mosaicism                            1
    """    

    # rename inheritance
    omim_df["phenotype_inheritance"] = omim_df["phenotype_inheritance"].map({
        "Autosomal recessive": "AR",
        "Autosomal dominant": "AD",
        "X-linked recessive": "XR",
        "X-linked dominant": "XD",
        "X-linked": "XR",
        "Autosomal dominant, Autosomal recessive": "AD/AR",
        "X-linked dominant, X-linked recessive": "XD/XR",
        "Somatic mutation, Autosomal dominant": "Somatic/AD",
        "Somatic mutation, Autosomal recessive": "Somatic/AR",
        "Somatic mutation, Autosomal dominant, Autosomal recessive": "Somatic/AD/AR",
        "Somatic mosaicism, Autosomal dominant": "Somatic/AD",
        "Somatic mosaicism, Autosomal recessive": "Somatic/AR",
        "Somatic mutation": "Somatic",
        "Mitochondrial": "MITO",
    })

    return omim_df[OUTPUT_COLUMNS]


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_omim_table()
    print(df)
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d")
    output_path = f"omim_{timestamp}.tsv"
    df.to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(df)} records to {output_path}")


    print(df["phenotype_inheritance"].value_counts())
"""
5/19/2025 OMIM Table:

Columns:

0: chrom
1: start
2: end
3: mim_number
4: phenotype_mim_number
5: phenotypic_series_number
6: phenotype_inheritance
7: gene_symbols
8: gene_id
9: gene_description
10: phenotype_description
11: date_created
12: date_updated
13: mouse_gene_id
14: oe_lof_upper
15: pLI
16: mis_z
17: text
18: comments


(19525, 19)

Example:

chrom                                                   1
start                                             1013497
end                                               1014540
mim_number                                         147571
phenotype_mim_number                               616126
phenotypic_series_number
phenotype_inheritance                 Autosomal recessive
gene_symbols                  G1P2,  IFI15,  IMD38, ISG15
gene_id                                   ENSG00000187608
gene_description            ISG15 ubiquitin-like modifier
phenotype_description                 Immunodeficiency 38
date_created
date_updated
mouse_gene_id                         Isg15 (MGI:1855694)
oe_lof_upper                                        1.691
pLI                                               0.40527
mis_z                                           -0.044129
text
comments
"""
