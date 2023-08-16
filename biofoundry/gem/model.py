import logging

import os
import time

import json
import string
from functools import reduce

import pandas as pd

from modelseedpy import MSBuilder, MSGenome

import cobra
from cobra.io import load_json_model, save_json_model

from biofoundry.base import BaseModelBuilder, BaseModelValidator


# Configure logging
logging.basicConfig(
    filename=os.path.basename(__file__).replace(".py", ".log"),
    filemode="w",
    format="%(asctime)s - %(filename)s:%(lineno)s - %(funcName)s - " + \
        "%(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO
)
LOGGER = logging.getLogger(__name__)


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
        Validate whether the model can be loaded.

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
        Validate the model according to a set of predefined checks.

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

        LOGGER.info(f"Validating model {model_path}")
        self.validate_loading(model_path)


class ModelBuilder(BaseModelBuilder):
    """
    Auxiliary class for building COBRA GEMs from annotated genomes.

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
        Build a model from its annotated genome.

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

        # Load annotated genome
        genome = MSGenome.from_fasta(
            genome_path,
            split=" "
        )

        LOGGER.info(f"Genome created from {genome_path}")
        LOGGER.info(f"Number of features: {len(genome.features)}")

        # Build metabolic model from annotated genome
        model = MSBuilder.build_metabolic_model(
            model_id=os.path.basename(genome_path),
            genome=genome,
            gapfill_media=None,
            template=None,
            index="0",
            allow_all_non_grp_reactions=True,
            annotate_with_rast=True,
            gapfill_model=True
        )

        LOGGER.info("Model successfully built")

        return model

    @staticmethod
    def rename_compartments(model_text: str) -> str:
        """
        Rename compartments to fit COBRA conventions.

        Parameters
        ----------
        model_text : str
            The full model as string.

        Returns
        -------
        _ : str
            The model's text with compartments renamed.

        Notes
        -----
        ModelSEEDpy creates compartments with numbers (e.g. "c0" instead of 
        "c") which cannot be recognized by COBRApy (which is used by MICOM).

        Examples
        --------
        None

        """

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

    def rename_metabolites(self, model_text: str) -> str:
        """
        Rename metabolites to add ModelSEED compound ID and abbreviation for
        better tracking of compounds.

        Parameters
        ----------
        model_text : str
            The full model as string.

        Returns
        -------
        _ : str
            The model's text with metabolites renamed.

        Examples
        --------
        None

        """

        modelseed_cpd = pd.read_table(
            os.path.join(
                self.config["paths"]["modelseed"],
                "compounds.tsv"
            ),
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

    def format_model(self, model_path: str) -> str:
        """
        Format model to fit COBRA conventions and better track compounds.

        Parameters
        ----------
        model_path : str
            The path to the model's JSON file.

        Returns
        -------
        model_path_new : str
            The formatted model's path.

        Examples
        --------
        None

        """

        LOGGER.info(f"Formatting model {model_path}")

        # Create new path for the formatted model
        model_path_new = model_path\
            .replace(".json", "_formatted.json")

        # Load model as text and rename compartments and metabolites
        with open(model_path, "r") as fh:
            model_text = fh.read()
            LOGGER.debug(f"Model string length: {len(model_text)}")

            LOGGER.info("Renaming compartments...")
            model_text = self.rename_compartments(model_text)
            LOGGER.debug(f"Compartments renamed: {len(model_text)}")

            LOGGER.info("Renaming metabolites...")
            model_text = self.rename_metabolites(model_text)
            LOGGER.debug(f"Metabolites renamed: {len(model_text)}")

        # Dump formatted model
        with open(model_path_new, "w") as fh:
            json.dump(
                obj=json.loads(model_text),
                fp=fh,
                indent=4,
                sort_keys=False
            )

        LOGGER.info(f"Saved formatted model to {model_path_new}")

        return model_path_new

    def build(self, metadata_df: pd.DataFrame) -> None:
        """
        Build all models from their respectives annotated genomes in 
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

        LOGGER.info("Building models from:")
        LOGGER.debug("\t" + metadata_df.to_string().replace("\n", "\n\t"))

        for _, row in metadata_df.iterrows():
            organism = row["Code"]

            LOGGER.info(f"Starting with organism {organism}")

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

            model_path_formatted = self.format_model(model_path)

            # Test saved model
            self.model_validator.validate(model_path_formatted)

            # Wait between calls
            time.sleep(0.5)
