import logging

import os
import time

import pandas as pd

from modelseedpy import MSBuilder, MSGenome

import cobra
from cobra.io import load_json_model, save_json_model

from biofoundry.base import BaseModelBuilder, BaseModelValidator


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


class ModelValidator(BaseModelValidator):
    """
    Auxiliary class for validating generated GEMs.

    Parameters
    ----------
    None

    Examples
    --------
    None

    """

    @staticmethod
    def validate_loading(model_path: str) -> None:
        """
        Validates whether the model can be loaded.

        Parameters
        ----------
        model_path : str
            The path to the model's JSON file.

        Returns
        -------
        None

        Examples
        --------
        None

        """
        try:
            load_json_model(model_path)
            LOGGER.info(f"Correctly loaded {model_path}")
        except FileNotFoundError:
            LOGGER.exception(f"Cannot load {model_path}")
        else:
            LOGGER.exception(f"Unhandled exception while loading {model_path}")

    def validate(self, model_path: str) -> None:
        """
        Validates the model according to a set of predefined checks.

        Parameters
        ----------
        model_path : str
            The path to the model's JSON file.

        Returns
        -------
        None

        Examples
        --------
        None

        """
        self.validate_loading(model_path)


class ModelBuilder(BaseModelBuilder):
    """
    Auxiliary class for validating generated GEMs.

    Parameters
    ----------
    config : dict
        The configuration dictionary.
    model_validator : ModelValidator
        The instance of ModelValidator for validating the models.

    Examples
    --------
    None

    """

    def __init__(self, config: dict, model_validator: ModelValidator) -> None:
        super().__init__()

        self.config = config
        self.model_validator = model_validator

    @staticmethod
    def build_model(genome_path: str) -> cobra.Model:
        """
        Builds a model from its annotated genome.

        Parameters
        ----------
        genome_path : str
            Path pointing to the annotated genome.

        Returns
        -------
        model : cobra.Model
            The generated model using ModelSEEDpy

        Examples
        --------
        None

        """

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

    def build(self, metadata_df: pd.DataFrame) -> None:
        """
        Builds all models from their respectives annotated genomes in 
        metadata_df.

        Parameters
        ----------
        metadata_df : pandas.DataFrame
            Metadata dataframe with all paths to the annotated genomes.

        Returns
        -------
        None

        Examples
        --------
        None

        """

        for _, row in metadata_df.iterrows():
            organism = row["Code"]

            # Build model
            genome_path = os.path.join(
                self.config["paths"]["genomes"],
                organism,
                row["Protein annotation file"]
            )
            model = self.build_model(genome_path)

            # Save model
            model_path = os.path.join(
                self.config["paths"]["models"],
                f"{organism}.json"
            )
            save_json_model(
                model=model,
                filename=model_path
            )

            # Test saved model
            self.model_validator.validate(model_path)

            # Wait between calls
            time.sleep(0.5)
