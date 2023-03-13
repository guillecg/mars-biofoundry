import os
import logging

import pandas as pd

from bioservices.kegg import KEGG


ORGANISM = "tez"
FILE = os.path.basename(__file__).replace(".py", "")


# Configure logging
logging.basicConfig(
    filename=f"{FILE}-{ORGANISM}.log",
    filemode="w",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO
)
logger = logging.getLogger(FILE)
logger.info(f"Starting with organism {ORGANISM}")


s = KEGG()
s.organism = ORGANISM

# KEGG retrieves all enzyme IDs although it filters pathway IDs by organism
logger.info(f"Number of enzymes: {len(s.enzymeIds)}")
logger.info(f"Number of pathways: {len(s.pathwayIds)}")


results_columns = [
    "ORGANISM",
    "PATHWAY_NAME",
    "KO_PATHWAY",
    "MODULE",
    "MODULE_NAME",
    "REACTION",
    "REACTION_SCHEMA"
]

# Write empty CSV to later append to
results_filepath = f"{FILE}-{ORGANISM}.csv"
pd.DataFrame(columns=results_columns).to_csv(
    results_filepath,
    header=True,
    index=False,
    sep=";",
    mode="w"
)

# To avoid errors (e.g. missing modules in pathways) follow this process:
# Pathway > KO pathway > modules
for pathway in s.pathwayIds:

    logger.info(f"Starting with pathway {pathway}")

    # Retrieve info for the given pathway
    pathway_res = s.get(pathway)
    pathway_res = s.parse(pathway_res)

    # Retrieve info for the given KO pathway
    ko_res = s.get(pathway_res["KO_PATHWAY"])
    ko_res = s.parse(ko_res)

    pathway_df = pd.DataFrame(columns=results_columns)

    # Retrieve info for all the modules in the KO pathway
    for module in ko_res["MODULE"].keys():

        logger.info(f"Processing module {module}")

        module_res = s.get(module)
        module_res = s.parse(module_res)

        # TODO: try extracting EC from module_res["ORTHOLOGY"]
        if "REACTION" not in module_res.keys():
            logger.warning(f"No reactions found for {module}")
            continue

        reaction_keys = module_res["REACTION"].keys()

        for key in reaction_keys:
            # Concatenate results to final dataframe
            row_df = pd.DataFrame({
                "ORGANISM": [pathway_res["ORGANISM"]],
                "PATHWAY_NAME": [pathway_res["NAME"]],
                "KO_PATHWAY": [pathway_res["KO_PATHWAY"]],
                "MODULE": [module],
                "MODULE_NAME": [module_res["NAME"]],
                "REACTION": [key],
                "REACTION_SCHEMA": [module_res["REACTION"][key]],
            })
            pathway_df = pd.concat(
                [pathway_df, row_df],
                axis=0,
                ignore_index=True
            )

    logger.info(f"Writing pathway {pathway} to file")

    pathway_df.to_csv(
        results_filepath,
        header=False,
        index=False,
        sep=";",
        mode="a"
    )
