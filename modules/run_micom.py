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

    modelseed_cpd = pd.read_table(
        modelseed_cpd_path,
        dtype=object # Avoid warnings
    )

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


def get_counts(
    model: dict,
    element: str
):
    if type not in ("metabolites", "reactions"):
        raise NotImplementedError

    return len(model[element])


def format_model(model_path: str) -> str:

    model_path_new = model_path\
        .replace(".json", "_formatted.json")

    # Load model as text and rename compartments and metabolites
    with open(model_path, "r") as fh:
        model_text = fh.read()
        model_text = rename_compartments(model_text)
        model_text = rename_metabolites(model_text)

    # Dump formatted model
    with open(model_path_new, "w") as fh:
        json.dump(
            obj=json.loads(model_text),
            fp=fh,
            indent=4,
            sort_keys=False
        )

    # Load model as dictionary and count reactions and metabolites
    # NOTE: do this after formatting to check for any errors in the process
    with open(model_path_new, "r") as fh:
        # See https://stackoverflow.com/questions/64268575/how-can-i-import-a-json-as-a-dict#comment113647180_64268609
        model_dict, _ = json.load(fh)

        n_reactions = get_counts(model=model_dict, element="reactions")
        n_metabolites = get_counts(model=model_dict, element="metabolites")

    return model_path_new, n_reactions, n_metabolites


MODELS_DIR = "../data/modelseedpy/"
OUT_DIR = "../data/micom/rio_tinto/amils_2023/"

# WARNING: Not working for more than one thread, caution is advised!
N_THREADS = 1


metadata_df = pd.read_csv("../data/genomes/genomes-metadata.csv")

# Drop rercords with missing data
metadata_df = metadata_df.dropna(axis=0)

taxonomy = []

for _, row in metadata_df.iterrows():
    species, organism = row[["Species", "Code"]].values

    # NOTE: ModelSEEDpy creates compartments with numbers (e.g. "c0" instead
    # of "c") which cannot be recognized by COBRApy (used by MICOM).
    model_path = os.path.join(
        MODELS_DIR, f"{organism}.json"
    )
    model_path_new, n_reactions, n_metabolites = format_model(model_path)

    taxonomy += [
        {
            "id": organism,
            "genus": species.split(" ")[0],
            "species": species,
            "reactions": n_reactions,
            "metabolites": n_metabolites,
            "file": model_path_new
        }
    ]

taxonomy = pd.DataFrame.from_records(taxonomy)

# ---------------------------------------------------------------------------- #
# Abundances

# Read 454 and Illumina abundances from supplementary materials
PAPER_DIR = "../data/amils2023/"

abundances_illumina_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
        "emi16291-sup-0005-datasets4.xlsx"
    ),
    sheet_name="Filtered OTUs",
    skiprows=11
)
abundances_roche_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
        "emi16291-sup-0006-datasets5.xlsx"
    ),
    sheet_name="DW Filtered OTUs",
    skiprows=11
)

# Filter out species also present in the drilling water (possible contamination)
# abundances_illumina_df = abundances_illumina_df[
#     (abundances_illumina_df["DW_RG"] == 0) & \
#     (abundances_illumina_df["IC"] == 0)
# ]
# abundances_roche_df = abundances_roche_df[
#     abundances_roche_df["DWδ"] == 0
# ]

# Group samples by genus to avoid repeats
abundances_illumina_df = abundances_illumina_df\
    .groupby("Genus", as_index=False)\
    .sum()
abundances_roche_df = abundances_roche_df\
    .groupby("Genus", as_index=False)\
    .sum()

# Get samples in long format for Illumina reads
abundances_illumina_df = pd.wide_to_long(
    df=abundances_illumina_df,
    stubnames="BH10",
    i="Genus",
    j="Sample",
    sep="-",
    suffix="\\d+"
)
abundances_illumina_df = abundances_illumina_df\
    .reset_index()\
    [["Genus", "Sample", "BH10"]]\
    .rename(columns={
        "BH10": "abundance",
        "Genus": "genus",
        "Sample": "sample_id"
    })

abundances_illumina_df["sample_id"] = abundances_illumina_df["sample_id"]\
    .apply(lambda row: f"BH10-{str(row)}-Illumina")\
    .astype(str)

# Get samples in long format for Roche reads
abundances_roche_df = pd.wide_to_long(
    df=abundances_roche_df,
    stubnames="BH10",
    i="Genus",
    j="Sample",
    sep="-",
    suffix="\\d+"
)
abundances_roche_df = abundances_roche_df\
    .reset_index()\
    [["Genus", "Sample", "BH10"]]\
    .rename(columns={
        "BH10": "abundance",
        "Genus": "genus",
        "Sample": "sample_id"
    })

abundances_roche_df["sample_id"] = abundances_roche_df["sample_id"]\
    .apply(lambda row: f"BH10-{str(row)}-Roche")

# Concatenate
abundances_df = pd.concat(
    [abundances_illumina_df, abundances_roche_df],
    axis=0,
    ignore_index=True
)

# TODO: divide counts for Tessaracoccus and Rhizobium!!!
taxonomy_abundances = pd.merge(
    left=taxonomy,
    right=abundances_df,
    on="genus",
    how="left"
)

# Combine ID column with sample ID to avoid errors in MICOM due to ID duplicates
taxonomy_abundances["id"] = \
    taxonomy_abundances["id"] + "_" + taxonomy_abundances["sample_id"]

# ---------------------------------------------------------------------------- #
# Communities

com = Community(taxonomy_abundances)
sol = com.cooperative_tradeoff()

# ---------------------------------------------------------------------------- #
# Manifest

manifest = build(
    taxonomy=taxonomy,
    out_folder=OUT_DIR,
    model_db=None,
    cutoff=1e-2,
    threads=N_THREADS
)

manifest = pd.read_csv(
    os.path.join(
        OUT_DIR,
        "manifest.csv"
    )
)
