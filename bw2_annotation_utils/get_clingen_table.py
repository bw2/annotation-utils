import pandas as pd
import requests
from io import StringIO

from bw2_annotation_utils.cache_utils import cache_data_table


def _get_clingen_table(url):
    """Download one of the ClinGen .csv tables and return it as a pandas DataFrame

    Args:
        url (str): for example "https://search.clinicalgenome.org/kb/gene-validity/download"
    Return:
        pandas DataFrame
    """
    r = requests.get(url)
    if not r.ok:
        raise Exception(f"Failed to download {url}: {r}")

    table_contents = r.content.decode('UTF-8')
    lines = table_contents.split("\n")
    header_line = lines[4]
    #print(header_line)
    table_contents = "\n".join([header_line] + lines[6:])
    return pd.read_csv(StringIO(table_contents))


@cache_data_table
def get_clingen_gene_disease_validity_table():
    "Download ClinGen gene-disease validity table and return it as a pandas DataFrame"
    return _get_clingen_table("https://search.clinicalgenome.org/kb/gene-validity/download")

@cache_data_table
def get_clingen_dosage_sensitivity_table():
    return _get_clingen_table("https://search.clinicalgenome.org/kb/gene-dosage/download")


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_clingen_gene_disease_validity_table()
    print("Gene-Disease Validity Table columns:")
    print(df)

    df = get_clingen_dosage_sensitivity_table()
    print("Dosage Sensitivity Table columns:")
    print(df)


"""
5/19/2025 ClinGen Gene-Disease Validity Table:

Columns:
0: GENE SYMBOL
1: GENE ID (HGNC)
2: DISEASE LABEL
3: DISEASE ID (MONDO)
4: MOI
5: SOP
6: CLASSIFICATION
7: ONLINE REPORT
8: CLASSIFICATION DATE
9: GCEP

Example:

GENE SYMBOL                                                        AARS1
GENE ID (HGNC)                                                   HGNC:20
DISEASE LABEL                 Charcot-Marie-Tooth disease axonal type 2N
DISEASE ID (MONDO)                                         MONDO:0013212
MOI                                                                   AD
SOP                                                                SOP10
CLASSIFICATION                                                Definitive
ONLINE REPORT          https://search.clinicalgenome.org/kb/gene-vali...
CLASSIFICATION DATE                             2024-03-14T16:00:00.000Z
GCEP                   Charcot-Marie-Tooth Disease Gene Curation Expe...


Unique GeneIds: 2,564


CLASSIFICATION
Definitive                       1932
Limited                           461
Moderate                          367
Disputed                          177
Strong                             70
No Known Disease Relationship      48
Refuted                            40


MOI
AR    1519
AD    1260
XL     182
MT      54
SD      51
UD      29
"""
