from datetime import datetime
import os
import pandas as pd
from bw2_annotation_utils.get_panel_app_table import get_panel_app_table
from bw2_annotation_utils.get_omim_table import get_omim_table
from bw2_annotation_utils.get_clingen_table import get_clingen_gene_disease_validity_table
from bw2_annotation_utils.get_hgnc_table import get_hgnc_table
from bw2_annotation_utils.get_panel_app_table import get_panel_app_table
from bw2_annotation_utils.get_ensembl_db_info import get_transcript_id_to_gene_id

#df = get_clingen_gene_disease_validity_table()

separtor = "; "

def normalize_nulls(x):
    return str(x) if not pd.isna(x) else ""

df_hgnc = get_hgnc_table()
HGNC_to_ENSG_map = dict(zip(df_hgnc["HGNC ID"], df_hgnc["Ensembl gene ID"]))
ENSG_to_gene_name_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["Approved symbol"]))
ENSG_to_gene_name_aliases_map = dict(zip(df_hgnc["Ensembl gene ID"], df_hgnc["Alias symbols"]))

"""
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
"""

df_omim = get_omim_table()
print(f"Got {len(df_omim):,d} rows from OMIM, containing {len(df_omim['gene_id'].unique()):,d} unique genes")
df_omim = df_omim[[
    "gene_id",   # ENSG id
    "mim_number", 
    "phenotype_mim_number", 
    "phenotypic_series_number", 
    "phenotype_inheritance", 
    "phenotype_description", 
    #"gene_symbols", 
    #"gene_description", 
    #"mouse_gene_id", 
    "oe_lof_upper", 
    "pLI", 
    "mis_z",
]]

df_omim = df_omim.rename(columns={
    "gene_id": "OMIM_gene_id",
    "mim_number": "OMIM_mim_number",
    "phenotype_mim_number": "OMIM_phenotype_mim_number",
    "phenotypic_series_number": "OMIM_phenotypic_series_number",
    "phenotype_inheritance": "OMIM_inheritance",
    "phenotype_description": "OMIM_phenotype_description",
})


before = len(df_omim)
df_omim = df_omim[(df_omim["OMIM_gene_id"] != '') & ((df_omim["OMIM_phenotype_mim_number"] != '') | (df_omim["OMIM_phenotype_description"] != ''))]
print("\t", f"Kept {len(df_omim):,d} out of {before:,d} ({(len(df_omim) / before):.1%}) rows which had both a gene and a phenotype")


# group by gene_id and combine the other fields using ; as a separator
df_omim = df_omim.groupby("OMIM_gene_id").agg({
    "OMIM_mim_number": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "OMIM_phenotype_mim_number": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "OMIM_phenotypic_series_number": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "OMIM_inheritance": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "OMIM_phenotype_description": lambda x: separtor.join(normalize_nulls(v) for v in x),
}).reset_index()


df_omim.set_index("OMIM_gene_id", inplace=True)


"""
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
"""
df_clingen = get_clingen_gene_disease_validity_table()
print(f"Got {len(df_clingen):,d} rows from ClinGen, containing {len(df_clingen['GENE ID (HGNC)'].unique()):,d} unique genes")

df_clingen = df_clingen.rename(columns={
    "DISEASE LABEL": "CLINGEN_disease_label",
    "DISEASE ID (MONDO)": "CLINGEN_disease_mondo_id",
    "MOI": "CLINGEN_inheritance",
    "CLASSIFICATION": "CLINGEN_classification",
})

df_clingen["CLINGEN_gene_id"] = df_clingen["GENE ID (HGNC)"].map(HGNC_to_ENSG_map)
hgnc_ids_with_missing_esng = df_clingen[df_clingen['CLINGEN_gene_id'].isna()]['GENE ID (HGNC)'].unique()
assert len(hgnc_ids_with_missing_esng) == 0, f"Could not convert the following HGNC ids to ENSG: {', '.join(hgnc_ids_with_missing_esng)}"

before = len(df_clingen)
df_clingen = df_clingen[df_clingen["CLINGEN_classification"].isin({
    "Definitive", "Limited", "Moderate", "Strong"
})]
print("\t", f"Kept {len(df_clingen):,d} out of {before:,d} ({(len(df_clingen) / before):.1%}) rows which had a Definitive, Limited, Moderate, or Strong classification")


df_clingen = df_clingen[[
    "CLINGEN_gene_id",
    "CLINGEN_disease_label",
    "CLINGEN_disease_mondo_id",
    "CLINGEN_inheritance",
    "CLINGEN_classification",
    #"GENE_SYMBOL",
    #"SOP",
    #"ONLINE_REPORT",
    #"CLASSIFICATION_DATE",
    #"GCEP",
]]

# group by CLINGEN_gene_id and combine the other fields using ; as a separator
df_clingen = df_clingen.groupby("CLINGEN_gene_id").agg({
    "CLINGEN_disease_label": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "CLINGEN_disease_mondo_id": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "CLINGEN_inheritance": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "CLINGEN_classification": lambda x: separtor.join(normalize_nulls(v) for v in x),
}).reset_index()

df_clingen.set_index("CLINGEN_gene_id", inplace=True)


"""
source                                               PanelApp UK
hgnc                                                  HGNC:24641
gene_name                    chromosome 16 open reading frame 62
biotype                                           protein_coding
gene_id                                          ENSG00000103544
confidence                                                     2
penetrance                                                  None
mode_of_pathogenicity                                       None
mode_of_inheritance      BIALLELIC, autosomal or pseudoautosomal
publications                                            31712251
evidence                         Expert Review Amber, Literature
phenotypes                    3C/Ritscher-Schinzel-like syndrome
panel_name                             Chondrodysplasia punctata
"""

df_panel_app = get_panel_app_table()
PANEL_APP_UK_LABEL = "PanelApp UK"
PANEL_APP_AU_LABEL = "PanelApp Australia"
assert set(df_panel_app["source"]) == {
    PANEL_APP_UK_LABEL, PANEL_APP_AU_LABEL,
}
print(f"Got {len(df_panel_app):,d} rows from PanelApp, containing {len(df_panel_app['gene_id'].unique()):,d} unique genes")

before = len(df_panel_app)
df_panel_app["gene_id"] = df_panel_app["gene_id"].fillna(df_panel_app["hgnc"].map(HGNC_to_ENSG_map))
df_panel_app = df_panel_app[df_panel_app["gene_id"].notna() & (df_panel_app["gene_id"] != "")]
print("\t", f"Kept {len(df_panel_app):,d} out of {before:,d} ({(len(df_panel_app) / before):.1%}) rows which had a gene id")

before = len(df_panel_app)
df_panel_app = df_panel_app[~df_panel_app["evidence"].apply(lambda x: not isinstance(x, str) or "Expert Review Red" in x)]
print("\t", f"Kept {len(df_panel_app):,d} out of {before:,d} ({(len(df_panel_app) / before):.1%}) rows which had evidence other than 'Expert Review Red'")

df_panel_app = df_panel_app[[
    "gene_id",
    "source",
    #"hgnc",
    #"gene_name",
    #"biotype",
    "confidence",
    "penetrance",
    "mode_of_pathogenicity",
    "mode_of_inheritance",
    #"publications",
    "evidence",
    "phenotypes",
    "panel_name",
]]

df_panel_app = df_panel_app.groupby(["gene_id", "source"]).agg({
    "confidence": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "penetrance": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "mode_of_pathogenicity": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "mode_of_inheritance": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "evidence": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "phenotypes": lambda x: separtor.join(normalize_nulls(v) for v in x),
    "panel_name": lambda x: separtor.join(normalize_nulls(v) for v in x),
}).reset_index()

# drop column "source"
df_panel_app_uk = df_panel_app[df_panel_app["source"] == PANEL_APP_UK_LABEL].drop("source", axis=1)
df_panel_app_au = df_panel_app[df_panel_app["source"] == PANEL_APP_AU_LABEL].drop("source", axis=1)
for panel_app_label, df_pannel_app in [
    ("UK", df_panel_app_uk),
    ("AU", df_panel_app_au),
]:
    df_pannel_app.rename(columns={
        "gene_id": f"PANEL_APP_{panel_app_label}_gene_id",
        "confidence": f"PANEL_APP_{panel_app_label}_confidence",
        "penetrance": f"PANEL_APP_{panel_app_label}_penetrance",
        "mode_of_pathogenicity": f"PANEL_APP_{panel_app_label}_mode_of_pathogenicity",
        "mode_of_inheritance": f"PANEL_APP_{panel_app_label}_inheritance",
        "evidence": f"PANEL_APP_{panel_app_label}_evidence",
        "phenotypes": f"PANEL_APP_{panel_app_label}_phenotypes",
        "panel_name": f"PANEL_APP_{panel_app_label}_panel_name",
    }, inplace=True)

    df_pannel_app.set_index(f"PANEL_APP_{panel_app_label}_gene_id", inplace=True)

    print("\t", f"PanelApp {panel_app_label} contains {len(df_pannel_app):,d} gene ids")
    # group by gene_id and combine the other fields using ; as a separator

# do an outer join of the 2 tables
df_panel_app = pd.merge(df_panel_app_uk, df_panel_app_au, how="outer", left_index=True, right_index=True)
print("\t", f"Merged PanelApp table contains {len(df_panel_app):,d} gene ids")


transcript_id_to_gene_id = get_transcript_id_to_gene_id()

fridman_path = "bw2_annotation_utils/data/AR_genes_from_Fridman_2025.tsv"
if os.path.exists(fridman_path):
    df_fridman = pd.read_table(fridman_path)
    df_fridman["FRIDMAN_gene_id"] = df_fridman["Transcripts"].apply(
        lambda transcript_list: ", ".join(transcript_id_to_gene_id[t] for t in transcript_list.split(",") if t in transcript_id_to_gene_id)
    )
    assert sum(df_fridman["FRIDMAN_gene_id"].str.contains(",")) == 0, "Some rows had multiple gene ids: " + str(df_fridman[df_fridman["FRIDMAN_gene_id"].str.contains(",")])

    df_fridman = df_fridman[df_fridman["FRIDMAN_gene_id"].notna() & (df_fridman["FRIDMAN_gene_id"] != "")]
    df_fridman = df_fridman[[
        "FRIDMAN_gene_id",
        #"Transcripts",
        "OMIM phenotype ID",
        "Disorder group",
        "Inheritance mode (AR/AR-AD)"
    ]]

    df_fridman.rename(columns={
        "OMIM phenotype ID": "FRIDMAN_omim_phenotype_id",
        "Disorder group": "FRIDMAN_phenotype_category",
        "Inheritance mode (AR/AR-AD)": "FRIDMAN_inheritance",
    }, inplace=True)


    df_fridman = df_fridman.groupby("FRIDMAN_gene_id").agg({
        "FRIDMAN_omim_phenotype_id": lambda x: separtor.join(normalize_nulls(v) for v in x),
        "FRIDMAN_phenotype_category": lambda x: separtor.join(normalize_nulls(v) for v in x),
        "FRIDMAN_inheritance": lambda x: separtor.join(normalize_nulls(v) for v in x),
    }).reset_index()

    df_fridman.set_index("FRIDMAN_gene_id", inplace=True)


# merge df_omim, df_clingen, df_panel_app, df_fridman
print(f"Merging OMIM ({len(df_omim):,d} rows), ClinGen ({len(df_clingen):,d} rows), PanelApp ({len(df_panel_app):,d} rows), and Fridman ({len(df_fridman):,d} rows)")
df_combined = pd.merge(df_omim, df_clingen, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging OMIM and ClinGen"
df_combined = pd.merge(df_combined, df_panel_app, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with PanelApp"
df_combined = pd.merge(df_combined, df_fridman, how="outer", left_index=True, right_index=True)
assert df_combined.index.is_unique, "The merged dataframe has duplicate gene ids after merging with Fridman"
df_combined.reset_index(inplace=True)
df_combined.rename(columns={
    "index": "gene_id",
}, inplace=True)

df_combined["gene_name"] = df_combined["gene_id"].map(ENSG_to_gene_name_map).str.upper()
df_combined["gene_aliases"] = df_combined["gene_id"].map(ENSG_to_gene_name_aliases_map).str.upper()

timestamp = datetime.now().strftime("%Y_%m_%d")
output_path = f"combined_mendelian_gene_disease_table.{len(df_combined)}_genes.{timestamp}.tsv"
df_combined.to_csv(output_path, sep="\t", index=False)
print(f"Wrote {len(df_combined):,d} genes to {output_path}")