import os

import csv
import pandas as pd


ORGANISM = "tez"
DATA_DIR = "../data/"
RETRORULES_DIR = "../../../papers/RetroRules/retrorules_rr01_rp2/"
KEGG_FILEPATH = os.path.join(DATA_DIR, f"kegg_{ORGANISM}.csv")


kegg_df = pd.read_csv(
    KEGG_FILEPATH,
    sep=","
)

retrorules_df = pd.read_csv(
    os.path.join(RETRORULES_DIR, "retrorules_rr01_rp2_flat_all.csv"),
    sep=","
)

# Explode multiple EC numbers per row
kegg_df["EC_NUMBER"] = kegg_df["EC_NUMBER"].str.split(";")
kegg_df = kegg_df.explode("EC_NUMBER").reset_index(drop=True)

retrorules_df["EC number"] = retrorules_df["EC number"].str.split(";")
retrorules_df = retrorules_df.explode("EC number").reset_index(drop=True)

merged_df = pd.merge(
    left=kegg_df,
    right=retrorules_df,
    left_on="EC_NUMBER",
    right_on="EC number",
    how="inner"
)

merged_df.to_csv(
    KEGG_FILEPATH.replace(".csv", "_rules.csv"),
    header=True,
    index=False,
    sep=",",
    quotechar='"',
    quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
)
