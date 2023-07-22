import os

import json

import pandas as pd


class RetroPathPreloader(object):

    @staticmethod
    def get_ec_from_model(
        model_path: str,
        modelseed_dir: str
    ) -> pd.DataFrame:
        """
        Get the EC numbers from the GEM model in JSON format.

        Parameters
        ----------
        model_path : str
            The path to the model's JSON file.
        modelseed_dir : str
            The path to the directory containing the ModelSEED database.

        Returns
        -------
        ec_numbers : DataFrame
            Dataframe containing each all EC numbers found for the model.

        Examples
        --------
        None

        """

        with open(model_path, mode="r") as fh:
            model_dict = json.loads(fh.read())

        # Get reaction IDs
        reaction_ids = [
            item["id"].split("_")[0]
            for item in model_dict["reactions"]
        ]

        # Get EC numbers from reactions in ModelSEEDDatabase
        modelseed_df = pd.read_table(
            os.path.join(modelseed_dir, "reactions.tsv")
        )

        ec_numbers = modelseed_df[
            modelseed_df["id"].isin(reaction_ids)
        ]["ec_numbers"]

        # Explode series since there may be multiple EC numbers per reaction
        ec_numbers = ec_numbers.str.split("|").explode()

        # Convert to frame
        ec_numbers = ec_numbers.to_frame()

        # Add model ID
        ec_numbers["ID"] = model_dict["id"]

        # TODO: log statistics

        return ec_numbers

    def get_rules(self) -> None:
        raise NotImplementedError

    def get_sink(self) -> None:
        raise NotImplementedError

    def get_sources(self) -> None:
        raise NotImplementedError


class RetroPathLauncher(object):

    def __init__(self) -> None:
        raise NotImplementedError

    def launch(self) -> None:
        raise NotImplementedError
