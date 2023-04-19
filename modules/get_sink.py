import os
import logging

import requests
from bs4 import BeautifulSoup

import csv
import pandas as pd


ORGANISM = "tez"
DATA_DIR = "../data/retropath/"
KEGG_FILEPATH = os.path.join(DATA_DIR, f"kegg_{ORGANISM}_rules.csv")


# Read MetaNetX reac_prop.tsv file
metanetx_reac_prop = pd.read_table(
    "../data/metanetx/reac_prop.tsv",
    comment='#', # Skip comment rows
    header=None,
    names=[
        "ID",
        "mnx_equation",
        "reference",
        "classifs",
        "is_balanced",
        "is_transport"
    ]
)

# Read MetaNetX chem_prop.tsv file
metanetx_chem_prop = pd.read_table(
    "../data/metanetx/chem_prop.tsv",
    comment='#', # Skip comment rows
    header=None,
    names=[
        "Name", # Format for sink.csv
        "compound",
        "reference",
        "formula",
        "charge",
        "mass",
        "InChI",
        "InChIKey",
        "SMILES"
    ]
)

# Load rules extracted by mapping KEGG's ECs to RetroRules DB (get_rules.py)
rules_df = pd.read_csv(
    KEGG_FILEPATH,
    sep=","
)

# Get all MetaNetX reaction IDs
rules_df[["MNXR", "MNXM"]]  = rules_df["Rule ID"].str.split("_", expand=True)

# Drop diameter and rule usage duplicates
unique_rules_df = rules_df.drop_duplicates(subset=["MNXR", "Rule usage"])

# Get all reaction information by MNXR
reactions_df = pd.merge(
    left=metanetx_reac_prop,
    right=unique_rules_df,
    left_on="ID",
    right_on="MNXR",
    how="outer",
    indicator=True
)

# Remove "left_only" reactions (not present in the organism rules)
reactions_df = reactions_df[reactions_df["_merge"] != "left_only"]

# Get all reactions that could not be mapped in METANETX
not_in_metanetx = reactions_df[reactions_df["_merge"] == "right_only"]
print(
    "[WARNING] #reactions not present in METANETX: ",
    len(not_in_metanetx), "/", len(reactions_df)
)

# Finally, get reactions in the organism that could be mapped to METANETX
reactions_df = reactions_df[reactions_df["_merge"] == "both"]

# NOTE: METANETX reactions are in forward format. Therefore, we should get the
# substrates for the sink.
compounds_forward = reactions_df["mnx_equation"]\
    .str.split("=").str[0]\
    .str.extractall(r"(MNXM\d+)")\
    .reset_index(drop=True)\
    .rename(columns={0: "Name"})

# Add products of reversible reactions
compounds_reverse = reactions_df[
        reactions_df["Rule usage"] == "both"
    ]["mnx_equation"]\
    .str.split("=").str[-1]\
    .str.extractall(r"(MNXM\d+)")\
    .reset_index(drop=True)\
    .rename(columns={0: "Name"})

# Merge substrates and products (reversible reactions)
compounds_df = pd.concat(
    [compounds_forward, compounds_reverse],
    axis=0,
    ignore_index=True
)

# Drop potential duplicates (compounds appearing in more than one reation)
compounds_df = compounds_df.drop_duplicates()

# Get InChIs for the extracted compound IDs
compounds_df = pd.merge(
    left=compounds_df,
    right=metanetx_chem_prop,
    on="Name",
    how="left"
)
compounds_df = compounds_df[["Name", "InChI"]]

# WARNING: there are compounds without InChI
no_inchi = compounds_df[compounds_df["InChI"].isnull()]
print(
    "[WARNING] #compounds without InChI: ",
    len(no_inchi), "/", len(compounds_df)
)

# Fill missing InChIs and match sink.csv format
compounds_df = compounds_df.fillna("None")

# Write to file
results_filepath = os.path.join(DATA_DIR, f"kegg_{ORGANISM}_sink.csv")
compounds_df.to_csv(
    results_filepath,
    header=True,
    index=False,
    sep=",",
    mode="w",
    quotechar='"',
    quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
)
