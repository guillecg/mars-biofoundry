import os

import csv
import pandas as pd


DATA_DIR = "../data/retropath/"
RETRORULES_DIR = "../data/retrorules/retrorules_rr01_rp2/"
EC_NUM_FILEPATH = os.path.join(DATA_DIR, "ec_numbers_all.csv")


ec_num_df = pd.read_csv(
    EC_NUM_FILEPATH,
    sep=","
)

retrorules_df = pd.read_csv(
    os.path.join(RETRORULES_DIR, "retrorules_rr01_rp2_flat_forward.csv"),
    sep=","
)

# Drop potential duplicates (more than one species with the same EC)
ec_numbers = ec_num_df["ec_numbers"].unique()

retrorules_df["EC number"] = retrorules_df["EC number"].str.split(";")
retrorules_df = retrorules_df.explode("EC number").reset_index(drop=True)

merged_df = retrorules_df[
    retrorules_df["EC number"].isin(ec_numbers)
]

merged_df.to_csv(
    EC_NUM_FILEPATH.replace(".csv", "_rules.csv"),
    header=True,
    index=False,
    sep=",",
    quotechar='"',
    quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
)
