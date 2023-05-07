import os

import json

import csv
import pandas as pd


# Read model
model_path = "../data/modelseedpy/aci_formatted.json"


def get_ec_from_model(model_path: str) -> pd.DataFrame:

    with open(model_path, mode="r") as fh:
        model_dict = json.loads(fh.read())

    # Get reaction IDs
    reaction_ids = [item["id"].split("_")[0] for item in model_dict["reactions"]]

    # Get EC numbers from reactions in ModelSEEDDatabase
    modelseed_df = pd.read_table("../data/modelseed/reactions.tsv")

    ec_numbers = modelseed_df[modelseed_df["id"].isin(reaction_ids)]["ec_numbers"]

    # Explode series since there may be multiple EC numbers per reaction
    ec_numbers = ec_numbers.str.split("|").explode()

    # Convert to frame
    ec_numbers = ec_numbers.to_frame()

    # Add model ID
    ec_numbers["ID"] = model_dict["id"]

    # TODO: log statistics

    return ec_numbers


MODELS_DIR = "../data/modelseedpy/"


metadata_df = pd.read_csv("../data/genomes/genomes-metadata.csv")

# Drop rercords with missing data
metadata_df = metadata_df.dropna(axis=0)

results_df = pd.DataFrame()

for organism in metadata_df["Code"]:

    model_path = os.path.join(
        MODELS_DIR, f"{organism}.json"
    )

    results_df = pd.concat(
        [results_df, get_ec_from_model(model_path)],
        axis=0,
        ignore_index=True
    )

results_df.to_csv(
    "../data/retropath/ec_numbers_all.csv",
    header=True,
    index=False,
    sep=",",
    mode="w",
    quotechar='"',
    quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
)
