import os
import logging

from io import StringIO

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


# TODO: check https://www.genome.jp/dbget-bin/get_linkdb?-t+genes+gn:T04747

# Get all genes for the given organism
genes = s.list(ORGANISM)

gene_cols = [
    "GENE_ID",
    "GENE_TYPE",
    "COORDINATES",
    "DEFINITION"
]
genes_df = pd.read_csv(
    StringIO(genes),
    header=None,
    sep="\t",
    names=gene_cols
)

# TODO: not working - find ECs
genes_df[genes_df["DESCRIPTION"].str.contains("\[EC\:")]
