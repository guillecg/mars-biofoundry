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
                f"{organism}.json"
            )

            # Load model as dictionary and count reactions and metabolites
            with open(model_path, "r") as fh:
                # See https://stackoverflow.com/questions/64268575/how-can-i-import-a-json-as-a-dict#comment113647180_64268609
                model_dict, _ = json.load(fh)

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
