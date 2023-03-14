import os
import re
import logging

import csv
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

# KEGG retrieves all enzyme and modules IDs although it filters pathway IDs by
# organism
logger.info(f"Number of enzymes: {len(s.enzymeIds)}")
logger.info(f"Number of pathways: {len(s.pathwayIds)}")


results_columns = [
    "GENE_ID",
    "GENE_NAME",
    "KO_NUMBER",
    "EC_NUMBER",
    "PATHWAY_ID",
    "PATHWAY_NAME",
    "KO_PATHWAY",
    "ORGANISM"
]

# Write empty CSV to later append to
results_filepath = f"{FILE}_{ORGANISM}.csv"
pd.DataFrame(columns=results_columns).to_csv(
    results_filepath,
    header=True,
    index=False,
    sep=",",
    mode="w"
)

# NOTE: In order to retrieve information from RetroRules, we need either the
# reaction numbers or the ECs. The bioservices package retrieves, among other
# information, modules for the given pathway in case they are available.
# Modules are important because they contain reaction numbers. However, some 
# pathways may lack any modules and, thus, the information regarding the 
# reactions. Genes (EC numbers) information is not present either.

for pathway in s.pathwayIds:

    logger.info(f"Starting with pathway {pathway}")

    # Retrieve info for the given pathway
    pathway_res = s.get(pathway)
    pathway_res = s.parse(pathway_res)

    pathway_df = pd.DataFrame(columns=results_columns)

    # There may be pathways without any genes (generic/top-level pathways)
    if "GENE" not in pathway_res.keys():
        logger.warning(f"No genes found for {pathway}")
        continue

    # Retrieve info for all the modules in the KO pathway
    for gene_id, gene_desc in pathway_res["GENE"].items():

        logger.info(f"Processing gene {gene_id}")

        # TODO: get information for each gene
        # gene_res = s.get(f"{ORGANISM}:{gene_id}")
        # gene_res = s.parse(gene_res)

        gene_desc = gene_desc.split(" [")
        gene_desc = [x.replace("]", "") for x in gene_desc]

        # Get KO and EC numbers by regex
        gene_ko = list(
            filter(lambda x: re.match(pattern="KO\:", string=x), gene_desc)
        )
        gene_ec = list(
            filter(lambda x: re.match(pattern="EC\:", string=x), gene_desc)
        )

        # Remove prefix
        gene_ko = [x.replace("KO:", "") for x in gene_ko]
        gene_ec = [x.replace("EC:", "") for x in gene_ec]

        # Join multiple numbers
        gene_ko = ";".join(gene_ko).replace(" ", ";")
        gene_ec = ";".join(gene_ec).replace(" ", ";")

        # Concatenate results to final dataframe
        row = pd.Series({
            "GENE_ID": gene_id,
            "GENE_NAME": gene_desc[0],
            "KO_NUMBER": gene_ko,
            "EC_NUMBER": gene_ec,
            "PATHWAY_ID": pathway,
            "PATHWAY_NAME": pathway_res["NAME"][0],
            "KO_PATHWAY": pathway_res["KO_PATHWAY"],
            "ORGANISM": pathway_res["ORGANISM"]
        })
        pathway_df = pd.concat(
            [pathway_df, row.to_frame().T],
            axis=0,
            ignore_index=True
        )

    logger.info(f"Writing pathway {pathway} to file")

    pathway_df.to_csv(
        results_filepath,
        header=False,
        index=False,
        sep=",",
        mode="a",
        quotechar='"',
        quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
    )
