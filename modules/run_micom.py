from typing import List

import os
import json
import string
from functools import reduce

import pandas as pd

from micom import Community
from micom.workflows import build


def rename_compartments(model_text: str) -> str:
    return model_text\
        .replace("_c0", "_c")\
        .replace("_e0", "_e")\
        .replace('"c0"', '"c"')\
        .replace('"e0"', '"e"')\
        .replace("\n", "")\
        .replace(
            '"c":""', '"c":"cytosol"'
        )\
        .replace(
            '"e":""', '"e":"extracellular"'
        )


def rename_metabolites(
    model_text: str,
    modelseed_cpd_path: str = "../data/modelseed/compounds.tsv"
) -> str:

    modelseed_cpd = pd.read_table(modelseed_cpd_path)

    # Avoid errors in COBRApy due to punctuation and whitespaces in names
    modelseed_cpd["abbreviation"] = modelseed_cpd["abbreviation"]\
        .str.replace(
            pat=r"[{}|\s]".format(string.punctuation),
            repl="-",
            regex=True
        )

    # Create new column to avoid abbreviation duplicates
    modelseed_cpd["abbreviation_id"] = \
        modelseed_cpd["id"] + "=" + modelseed_cpd["abbreviation"]

    return reduce(
        lambda a, kv: \
            a.replace(*kv),
            modelseed_cpd[["id", "abbreviation_id"]].to_numpy(dtype=str),
            model_text
    )


def format_model(model_path: str) -> str:

    model_path_new = model_path\
        .replace(".json", "_formatted.json")

    with open(model_path, "r") as fh:
        model_text = fh.read()
        model_text = rename_compartments(model_text)
        model_text = rename_metabolites(model_text)

    with open(model_path_new, "w") as fh:
        json.dump(
            obj=json.loads(model_text),
            fp=fh,
            indent=4,
            sort_keys=False
        )

    return model_path_new


MODELS_DIR = "../data/modelseedpy/"
MODEL_DIR = "../data/micom/rio_tinto/amils_2023/"

# WARNING: Not working for more than one thread, caution is advised!
N_THREADS = 1


metadata_df = pd.read_csv("../data/genomes/genomes-metadata.csv")

# Drop rercords with missing data
metadata_df = metadata_df.dropna(axis=0)

taxonomy = []

for _, row in metadata_df.iterrows():
    species, organism = row[["Species", "Code"]]

    # NOTE: ModelSEEDpy creates compartments with numbers (e.g. "c0" instead
    # of "c") which cannot be recognized by COBRApy (used by MICOM).
    model_path = os.path.join(
        MODELS_DIR, f"{organism}.json"
    )
    model_path_new = format_model(model_path)

    taxonomy += [
        {
            "id": organism,
            "genus": species.split(" ")[0],
            "species": species,
            "reactions": None,
            "metabolites": None,
            "sample_id": "amils_2023",
            "file": model_path_new
        }
    ]

taxonomy = pd.DataFrame.from_records(taxonomy)

# ---------------------------------------------------------------------------- #
# Communities

# com = Community(taxonomy)
# sol = com.cooperative_tradeoff()

# ---------------------------------------------------------------------------- #
# Manifest

manifest = build(
    taxonomy=taxonomy,
    out_folder=MODEL_DIR,
    model_db=None,
    cutoff=1e-2,
    threads=N_THREADS
)

manifest = pd.read_csv(
    os.path.join(
        MODEL_DIR,
        "manifest.csv"
    )
)
