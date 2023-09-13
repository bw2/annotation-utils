import collections
import pymysql


# NOTE:   SQL db schema @ http://useast.ensembl.org/info/docs/api/core/core_schema.html
# To list the latest database name, run: mysql -h useastdb.ensembl.org -u anonymous -e "show databases;" | grep homo_sapiens_core
# To list available tables, run: mysql -D homo_sapiens_core_110_38 -u anonymous -h useastdb.ensembl.org -e "show tables"
# or use the mysql shell and SHOW DATABASES; SHOW TABLES; DESCRIBE <table>



CURRENT_ENSEMBL_DATABASE = "homo_sapiens_core_110_38"
ENSEMBL_HOST = "useastdb.ensembl.org"

"""
homo_sapiens_cdna_102_38
homo_sapiens_cdna_103_38
homo_sapiens_core_102_38
homo_sapiens_core_103_38
homo_sapiens_funcgen_102_38
homo_sapiens_funcgen_103_38
homo_sapiens_otherfeatures_102_38
homo_sapiens_otherfeatures_103_38
homo_sapiens_rnaseq_102_38
homo_sapiens_rnaseq_103_38
homo_sapiens_variation_102_38
homo_sapiens_variation_103_38
"""


def get_gene_name_to_gene_id(
        database=CURRENT_ENSEMBL_DATABASE,
        only_protein_coding=False,
        only_canonical_transcripts=False):
    """Retrieves a dictionary containing gene_name => ENSG gene_id.

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): If True, only return protein-coding genes and protein-coding transcripts
        only_canonical_transcripts (bool): If True, only return canonical transcripts

    Return:
        dict: mapping gene name to ENSG id string
    """

    gene_id_to_transcript_id = collections.defaultdict(list)
    with pymysql.connect(host=ENSEMBL_HOST, user="anonymous", database=database) as conn:
        with conn.cursor() as cursor:
            if only_canonical_transcripts:
                join_clause = "canonical_transcript_id = transcript_id"
            else:
                join_clause = "transcript.gene_id = gene.gene_id"

            columns = [
                # Gene fields
                "gene.stable_id",
                "gene.biotype",
                "gene.created_date",
                "gene.modified_date",
                # Transcript fields
                "transcript.stable_id",
                "transcript.biotype",
                "transcript.created_date",
                "transcript.modified_date",
            ]

            columns_str = ", ".join(columns)
            query_string = f"SELECT {columns_str} FROM gene LEFT JOIN transcript ON {join_clause}"
            if only_protein_coding:
                query_string += " WHERE gene.biotype = 'protein_coding' AND transcript.biotype = 'protein_coding'"

            cursor.execute(query_string)

            for row in cursor:
                gene_and_transcript_info = dict(zip(columns, row))
                gene_id = gene_and_transcript_info['gene.stable_id']
                gene_id_to_transcript_id[gene_id].append(gene_and_transcript_info)

    return gene_id_to_transcript_id


def get_gene_id_to_transcript_metadata(
        database=CURRENT_ENSEMBL_DATABASE,
        only_protein_coding=False,
        only_canonical_transcripts=False):
    """Retrieves a dictionary containing gene_id => a list of dictionaries each of which
    contains information about one transcript that belongs to that gene.

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): If True, only return protein-coding genes and protein-coding transcripts
        only_canonical_transcripts (bool): If True, only return canonical transcripts

    Return:
        dict: mapping ENSG id string to a list of dictionaries where each dictionary contains metadata fields
    """

    gene_id_to_transcript_id = collections.defaultdict(list)
    with pymysql.connect(host=ENSEMBL_HOST, user="anonymous", database=database) as conn:
        with conn.cursor() as cursor:
            if only_canonical_transcripts:
                join_clause = "canonical_transcript_id = transcript_id"
            else:
                join_clause = "transcript.gene_id = gene.gene_id"

            columns = [
                # Gene fields
                "gene.stable_id",
                "gene.biotype",
                "gene.created_date",
                "gene.modified_date",
                # Transcript fields
                "transcript.stable_id",
                "transcript.biotype",
                "transcript.created_date",
                "transcript.modified_date",
            ]

            columns_str = ", ".join(columns)
            query_string = f"SELECT {columns_str} FROM gene LEFT JOIN transcript ON {join_clause}"
            if only_protein_coding:
                query_string += " WHERE gene.biotype = 'protein_coding' AND transcript.biotype = 'protein_coding'"

            cursor.execute(query_string)

            for row in cursor:
                gene_and_transcript_info = dict(zip(columns, row))
                gene_id = gene_and_transcript_info['gene.stable_id']
                gene_id_to_transcript_id[gene_id].append(gene_and_transcript_info)

    return gene_id_to_transcript_id


def get_gene_id_to_transcript_ids(
        database=CURRENT_ENSEMBL_DATABASE,
        only_protein_coding=False,
        only_canonical_transcripts=False):
    """Returns a dictionary mapping each Ensembl gene_id => a list of transcript ids for that gene.

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): If True, only return protein coding genes
        only_canonical_transcripts (bool): If True, only return canonical transcripts
    Return:
        dict: mapping ENSG id string to a list of ENST id strings
    """

    gene_id_to_transcript_metadata_list = get_gene_id_to_transcript_metadata(
        database=database,
        only_canonical_transcripts=only_canonical_transcripts,
        only_protein_coding=only_protein_coding)

    return {
        gene_id: [
            transcript_metadata["transcript.stable_id"] for transcript_metadata in transcript_metadata_list
        ]
        for gene_id, transcript_metadata_list in gene_id_to_transcript_metadata_list.items()
    }


def get_gene_id_to_canonical_transcript_id(database=CURRENT_ENSEMBL_DATABASE, only_protein_coding=False):
    """Returns a dictionary mapping each Ensembl gene_id => canonical transcript id

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): Only include protein-coding genes

    Return:
        dict: mapping ENSG id string to the canonical ENST id string
    """

    gene_id_to_transcript_metadata_list = get_gene_id_to_transcript_metadata(
        database=database,
        only_canonical_transcripts=True,
        only_protein_coding=only_protein_coding)

    gene_id_to_canonical_transcript_id = {}
    for gene_id, transcript_metadata_list in gene_id_to_transcript_metadata_list.items():
        if len(transcript_metadata_list) > 1:
            raise Exception(f"{gene_id} has more than 1 canonical transcript")
        if len(transcript_metadata_list) == 0:
            raise Exception(f"{gene_id} has 0 canonical transcripts")

        gene_id_to_canonical_transcript_id[gene_id] = transcript_metadata_list[0]["transcript.stable_id"]

    return gene_id_to_canonical_transcript_id


def get_gene_created_modified_dates(
        database=CURRENT_ENSEMBL_DATABASE,
        only_protein_coding=False,
        only_canonical_transcripts=False):
    """Returns a dictionary mapping each Ensembl gene_id => a 2-tuple containing the created date and the modified date
    for that gene.

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): If True, only return protein coding genes
        only_canonical_transcripts (bool): If True, only return canonical transcripts
    Return:
        dict: mapping ENSG id string to a 2-tuple containing the created date and the modified date for that gene.
    """

    gene_id_to_transcript_metadata_list = get_gene_id_to_transcript_metadata(
        database=database,
        only_canonical_transcripts=only_canonical_transcripts,
        only_protein_coding=only_protein_coding)

    return {
        gene_id: (transcript_metadata_list[0]["gene.created_date"], transcript_metadata_list[0]["gene.modified_date"])
        for gene_id, transcript_metadata_list in gene_id_to_transcript_metadata_list.items()
    }


def get_transcript_created_modified_dates(
        database=CURRENT_ENSEMBL_DATABASE,
        only_protein_coding=False,
        only_canonical_transcripts=False):
    """Returns a dictionary mapping each Ensembl gene_id => a list of 3-tuples each containing the
    transcript id, created date, and modified date of a transcript for that gene.

    Args:
        database (str): The Ensembl database name (eg. "homo_sapiens_core_107_38")
        only_protein_coding (bool): If True, only return protein coding genes
        only_canonical_transcripts (bool): If True, only return canonical transcripts
    Return:
        dict: a dictionary mapping each Ensembl gene_id => a list of 3-tuples each containing the transcript id,
        created date, and modified date of a transcript for that gene.
    """

    gene_id_to_transcript_metadata_list = get_gene_id_to_transcript_metadata(
        database=database,
        only_canonical_transcripts=only_canonical_transcripts,
        only_protein_coding=only_protein_coding)

    return {
        gene_id: [
            (
                transcript_metadata['transcript.stable_id'],
                transcript_metadata['transcript.created_date'],
                transcript_metadata['transcript.modified_date']
            )
            for transcript_metadata in transcript_metadata_list
        ]
        for gene_id, transcript_metadata_list in gene_id_to_transcript_metadata_list.items()
    }

def get_ensembl_ENST_to_RefSeq_ids(database=CURRENT_ENSEMBL_DATABASE):

    db = pymysql.connect(host=ENSEMBL_HOST, user="anonymous", database=database)
    cursor = db.cursor()
    cursor.execute(" ".join([
        "SELECT",
            "transcript.stable_id, xref.display_label",
        "FROM",
            "transcript, object_xref, xref,external_db",
        "WHERE",
            "transcript.transcript_id = object_xref.ensembl_id",
            "AND object_xref.ensembl_object_type = 'Transcript'",
            "AND object_xref.xref_id = xref.xref_id",
            "AND xref.external_db_id = external_db.external_db_id",
            "AND external_db.db_name = 'RefSeq_mRNA'",
    ]))

    ensembl_ENST_to_RefSeq_ids = collections.defaultdict(list)
    for ensembl_ESNT_id_without_version, refseq_NM_id in cursor.fetchall():
        ensembl_ENST_to_RefSeq_ids[ensembl_ESNT_id_without_version].append(refseq_NM_id)

    cursor.close()
    db.close()

    return ensembl_ENST_to_RefSeq_ids

