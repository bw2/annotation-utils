import pandas as pd


MANE_SUMMARY_TABLE_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.summary.txt.gz"


def get_MANE_ensembl_transcript_table(mane_summary_table_url=MANE_SUMMARY_TABLE_URL):
    """Download the MANE summary table and return it as a pandas DataFrame. Columns and example values are:

        $1   #NCBI_GeneID : GeneID:15
        $2   Ensembl_Gene : ENSG00000129673.10
        $3        HGNC_ID : HGNC:19
        $4         symbol : AANAT
        $5           name : aralkylamine N-acetyltransferase
        $6     RefSeq_nuc : NM_001088.3
        $7    RefSeq_prot : NP_001079.1
        $8    Ensembl_nuc : ENST00000392492.8
        $9   Ensembl_prot : ENSP00000376282.2
        $10   MANE_status : MANE Select
        $11    GRCh38_chr : NC_000017.11
        $12     chr_start : 76467603
        $13       chr_end : 76470117
        $14    chr_strand : +
    """
    return pd.read_table(mane_summary_table_url)
