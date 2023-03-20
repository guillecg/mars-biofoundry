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

# Get all reaction information by MNXR
reactions_df = metanetx_reac_prop[
    metanetx_reac_prop["ID"].isin(rules_df["MNXR"].unique())
]

# Extract compound IDs from mnx_equation
compounds_df = reactions_df["mnx_equation"]\
    .str.extractall(r"(MNXM\d+)")\
    .reset_index(drop=True)\
    .rename(columns={0: "Name"})

# Drop duplicates
compounds_df = compounds_df.drop_duplicates()

# Get InChIs for the extracted compound IDs
compounds_df = pd.merge(
    left=compounds_df,
    right=metanetx_chem_prop,
    on="Name",
    how="left"
)
compounds_df = compounds_df[["Name", "InChI"]]

compounds_df[compounds_df["InChI"].isnull()]

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
