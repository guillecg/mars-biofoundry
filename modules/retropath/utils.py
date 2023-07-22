import os

import json

import pandas as pd


class RetroPathPreloader(object):

    def __init__(self, config: dict) -> None:
        self.config = config

    def get_ec_from_model(self, model_path: str) -> pd.DataFrame:
        """
        Get the EC numbers from the GEM model in JSON format.

        Parameters
        ----------
        model_path : str
            The path to the model's JSON file.

        Returns
        -------
        ec_numbers : DataFrame
            Dataframe containing each all EC numbers found for the model.

        Examples
        --------
        None

        """

        # Load model
        with open(model_path, mode="r") as fh:
            model_dict = json.loads(fh.read())

        # Get reaction IDs
        reaction_ids = [
            item["id"].split("_")[0]
            for item in model_dict["reactions"]
        ]

        # Get EC numbers from reactions in ModelSEEDDatabase
        modelseed_df = pd.read_table(
            os.path.join(
                self.config["paths"]["modelseed"],
                "reactions.tsv"
            )
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

    def get_rules(self) -> pd.DataFrame:
        """
        Get the rules by filtering RetroRules database by the EC numbers found
        in the community.

        Parameters
        ----------
        None

        Returns
        -------
        rules_df : DataFrame
            Dataframe containing all rules found in the community.

        Examples
        --------
        None

        """

        # Load EC numbers in the community
        ec_num_df = pd.read_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                self.config["retropath"]["files"]["ec_numbers"]
            ),
            sep=","
        )

        # Load RetroRules database
        retrorules_df = pd.read_csv(
            self.config["paths"]["retrorules"],
            sep=","
        )

        # Drop potential duplicates (more than one species with the same EC)
        ec_numbers = ec_num_df["ec_numbers"].unique()

        # Explode rules with multiple EC numbers
        retrorules_df["EC number"] = retrorules_df["EC number"].str.split(";")
        retrorules_df = retrorules_df\
            .explode("EC number")\
            .reset_index(drop=True)

        # Filter by found EC numbers in the community
        rules_df = retrorules_df[
            retrorules_df["EC number"].isin(ec_numbers)
        ]

        return rules_df

    def get_sink(self) -> None:
        raise NotImplementedError

    def get_sources(self) -> None:
        raise NotImplementedError


class RetroPathLauncher(object):

    def __init__(self) -> None:
        raise NotImplementedError

    def launch(self) -> None:
        raise NotImplementedError
