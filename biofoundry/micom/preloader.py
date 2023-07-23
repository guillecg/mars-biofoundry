import os
import json

import pandas as pd

from biofoundry.base import BaseMICOMPreloader


class MICOMPreloader(BaseMICOMPreloader):
    """
    Auxiliary class for generating MICOM inputs.

    Parameters
    ----------
    config : dict
        The configuration dictionary.

    Examples
    --------
    None

    """

    def __init__(self, config: dict) -> None:
        super().__init__()

        self.config = config

    @staticmethod
    def get_counts(model: dict, element: str) -> int:
        """
        Count metabolites and/or compartments in COBRA models.

        Parameters
        ----------
        model : dict
            The selected model as dictionary.
        element : str
            The selected element ("metabolites" or "reactions").

        Returns
        -------
        _ : int
            The number of elements found in the model.

        Examples
        --------
        None

        """

        if type not in ("metabolites", "reactions"):
            raise NotImplementedError

        return len(model[element])

    def get_taxonomy(self, metadata_df: pd.DataFrame) -> pd.DataFrame:
        """
        Get taxonomies from model paths.

        Parameters
        ----------
        metadata_df : pandas.DataFrame
            Metadata dataframe with all paths to the annotated genomes.

        Returns
        -------
        _ : pandas.DataFrame
            The taxonomy dataframe for MICOM.

        Examples
        --------
        None

        """

        for _, row in metadata_df.iterrows():
            species, organism = row[["Species", "Code"]].values

            model_path = os.path.join(
                self.config["paths"]["models"],
                f"{organism}_formatted.json"
            )

            # Load model as dictionary and count reactions and metabolites
            with open(model_path, "r") as fh:
                model_dict = json.load(fh)

                n_reactions = self.get_counts(
                    model=model_dict,
                    element="reactions"
                )
                n_metabolites = self.get_counts(
                    model=model_dict,
                    element="metabolites"
                )

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

        return pd.DataFrame.from_records(taxonomy)
