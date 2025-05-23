import pandas as pd
import requests
from io import StringIO

from bw2_annotation_utils.cache_utils import cache_data_table

URL = "https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_ensembl_id&col=gd_mgd_id&col=gd_pubmed_ids&col=gd_locus_type&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"

@cache_data_table
def get_hgnc_table():
    """Download the HGNC table from https://www.genenames.org/download/custom/ and return it as a pandas DataFrame"""
    r = requests.get(URL)
    if not r.ok:
        raise Exception(f"Failed to download {URL}: {r}")

    table_contents = StringIO(r.content.decode('UTF-8'))
    return pd.read_table(table_contents)


if __name__ == "__main__":
    pd.set_option('display.max_columns', 500)

    df = get_hgnc_table()
    print(df)


"""
5/19/2025HGNC Table:

0: HGNC ID
1: Approved symbol
2: Approved name
3: Status
4: Previous symbols
5: Alias symbols
6: Chromosome
7: Accession numbers
8: RefSeq IDs
9: Ensembl gene ID
10: Mouse genome database ID
11: Pubmed IDs
12: Locus type

(44109, 13)

Example:

HGNC ID                                        HGNC:5
Approved symbol                                  A1BG
Approved name                  alpha-1-B glycoprotein
Status                                       Approved
Previous symbols                                  NaN
Alias symbols                                     NaN
Chromosome                                   19q13.43
Accession numbers                                 NaN
RefSeq IDs                                  NM_130786
Ensembl gene ID                       ENSG00000121410
Mouse genome database ID                  MGI:2152878
Pubmed IDs                                    2591067
Locus type                  gene with protein product

"""