import os

import json

import pandas as pd

from micom import Community


def fix_compartments(model_path: str) -> str:

    model_path_new = model_path\
        .replace(".json", "_corrected.json")

    with open(model_path, "r") as fh:
        model_text = fh.read()
        model_text = model_text\
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

    with open(model_path_new, "w") as fh:
        json.dump(
            obj=json.loads(model_text),
            fp=fh,
            indent=4,
            sort_keys=False
        )

    return model_path_new


DATA_DIR = "../data/modelseedpy/"

# WARNING: Not working for more than one thread, caution is advised!
# N_THREADS = 1


SPECIES_DICT = {
    "tez": "Tessaracoccus sp. T2.5-30"
}

taxonomy = []

for organism, species in SPECIES_DICT.items():

    # NOTE: ModelSEEDpy creates compartments with numbers (e.g. "c0" instead
    # of "c") which cannot be recognized by COBRApy (used by MICOM).
    model_path = os.path.join(
        DATA_DIR, f"{organism}.json"
    )
    model_path_new = fix_compartments(model_path)

    taxonomy += [
        {
            "id": organism,
            "genus": species.split(" ")[0],
            "species": species,
            "reactions": None,
            "metabolites": None,
            "file": model_path_new
        }
    ]

taxonomy = pd.DataFrame.from_records(taxonomy)

com = Community(taxonomy)
