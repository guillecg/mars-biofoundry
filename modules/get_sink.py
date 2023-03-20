import os
import logging

import requests
from bs4 import BeautifulSoup

import csv
import pandas as pd


ORGANISM = "tez"
DATA_DIR = "../data/retropath/"
KEGG_FILEPATH = os.path.join(DATA_DIR, f"kegg_{ORGANISM}_rules.csv")
FILE = os.path.basename(__file__).replace(".py", "")


# Configure logging
logging.basicConfig(
    filename=f"{FILE}_{ORGANISM}.log",
    filemode="w",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO
)
logger = logging.getLogger(FILE)
logger.info(f"Starting with organism {ORGANISM}")


results_columns = [
    "Name",
    "InChI"
]

# Write empty CSV to later append to
results_filepath = os.path.join(DATA_DIR, f"kegg_{ORGANISM}_sink.csv")
pd.DataFrame(columns=results_columns).to_csv(
    results_filepath,
    header=True,
    index=False,
    sep=",",
    mode="w",
    quotechar='"',
    quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
)

# Load rules extracted by mapping KEGG's ECs to RetroRules DB (get_rules.py)
rules_df = pd.read_csv(
    KEGG_FILEPATH,
    sep=","
)

# Get all MetaNetX reaction IDs
rules_df[["MNXR", "MNXM"]]  = rules_df["Rule ID"].str.split("_", expand=True)

reaction_id_list = rules_df["MNXR"].unique()

logger.info(f"Number of unique reaction IDs: {len(reaction_id_list)}")

# NOTE: compound IDs need to be extracted first from each reaction entry. Then,
# information fo each compound can be retrieved.

# Get all compounds for the given reaction IDs
for reaction_id in reaction_id_list:

    logger.info(f"Processing reaction {reaction_id}")

    html = requests.get(f"https://www.metanetx.org/equa_info/{reaction_id}")

    if "Not found (or deleted):" in html.text:
        logger.warning(f"No information found for {reaction_id}")
        continue

    soup = BeautifulSoup(html.text, features="lxml")

    # Get equation text by finding the equation td and extracting the next td
    equation_text = soup\
        .find("td", string="equation")\
        .find_next_sibling("td")\
        .text

    # Extract compound names by prefix MNXM
    compound_list = [
        elem.split("@")[0] # Remove compartment
        for elem in equation_text.split(" ")
        if elem.startswith("MNXM")
    ]

    logger.info(f"Compounds found for {reaction_id}: {compound_list}")

    # Create empty dataframe to append
    compounds_df = pd.DataFrame(columns=results_columns)

    # Get information for each compound
    for compound_id in compound_list:
        html = requests.get(f"https://www.metanetx.org/chem_info/{compound_id}")

        if "Not found (or deleted):" in html.text:
            logger.warning(f"No information found for {compound_id}")
            continue

        soup = BeautifulSoup(html.text, features="lxml")

        # Get InChI text by finding the td and extracting the next td
        inchi = soup\
            .find("td", string="InChI")\
            .find_next_sibling("td")\
            .text

        if not inchi.startswith("InChI"):
            logger.warning(
                f"InChI for compound {compound_id} is malformed: {inchi}"
            )
            continue

        # Concatenate results to final dataframe
        row = pd.Series({
            "Name": compound_id,
            "InChI": inchi
        })
        compounds_df = pd.concat(
            [compounds_df, row.to_frame().T],
            axis=0,
            ignore_index=True
        )

    logger.info(f"Writing reaction {reaction_id} to file")

    compounds_df.to_csv(
        results_filepath,
        header=False,
        index=False,
        sep=",",
        mode="a",
        quotechar='"',
        quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
    )


# Finally, drop duplicates
compounds_df = pd.read_csv(
    results_filepath,
    sep=","
)

original_len = len(compounds_df)
compounds_df = compounds_df.drop_duplicates()
final_len = len(compounds_df)

logger.info(f"Dropping duplicates: from {original_len} to {final_len} records")

logger.info("Overwriting existing sink file...")
compounds_df.to_csv(
    results_filepath,
    header=False,
    index=False,
    sep=",",
    mode="w",
    quotechar='"',
    quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
)
logger.info("Finished!")
