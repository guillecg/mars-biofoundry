import logging

import os

from modelseedpy import MSBuilder, MSGenome

import cobra
from cobra.io import load_json_model


FILE = os.path.basename(__file__).replace(".py", "")

# Configure logging
logging.basicConfig(
    filename=f"{FILE}.log",
    filemode="w",
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO
)
LOGGER = logging.getLogger(FILE)


def build_model(genome_path: str) -> cobra.Model:
    LOGGER.info(f"Starting with genome")

    # Load annotated genome
    genome = MSGenome.from_fasta(
        genome_path,
        split=" "
    )

    LOGGER.info(f"Number of features: {len(genome.features)}")

    # Build metabolic model from annotated genome
    model = MSBuilder.build_metabolic_model(
        model_id=os.path.basename(genome_path),
        genome=genome,
        gapfill_media=None,
        allow_all_non_grp_reactions=True,
        annotate_with_rast=True
    )

    return model


def validate_model(model_path: str) -> None:
    try:
        load_json_model(model_path)
        LOGGER.info(f"Correctly loaded {model_path}")
    except FileNotFoundError:
        LOGGER.exception(f"Cannot load {model_path}")
    else:
        LOGGER.exception(f"Unhandled exception while loading {model_path}")
