import os

import json

import csv
import pandas as pd


class RetroPathPreloader(object):

    def __init__(self, config: dict) -> None:
        self.config = config

    @staticmethod
    def get_ec_from_model(
        model_dict: dict,
        modelseed_reactions: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Get the EC numbers from a GEM model.

        Parameters
        ----------
        model_dict : dict
            The model as a dictionary.
        modelseed_reactions : DataFrame
            ModelSEED reactions mapping model IDs (RXNs) to EC numbers.

        Returns
        -------
        ec_numbers_df : DataFrame
            Dataframe containing each all EC numbers found for the model.

        Examples
        --------
        None

        """

        # Get reaction IDs
        reaction_ids = [
            item["id"].split("_")[0]
            for item in model_dict["reactions"]
        ]

        ec_numbers_df = modelseed_reactions[
            modelseed_reactions["id"].isin(reaction_ids)
        ]["ec_numbers"]

        # Explode series since there may be multiple EC numbers per reaction
        ec_numbers_df = ec_numbers_df.str.split("|").explode()

        # Convert to frame
        ec_numbers_df = ec_numbers_df.to_frame()

        # Add model ID
        ec_numbers_df["ID"] = model_dict["id"]

        # TODO: log statistics

        return ec_numbers_df

    def get_ec_numbers(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Get all the EC numbers for each model specified in the metadata.

        Parameters
        ----------
        metadata : str
            The dataframe containing the models' paths.

        Returns
        -------
        all_ec_numbers_df : DataFrame
            Dataframe containing each all EC numbers found for the model.

        Examples
        --------
        None

        """

        # Get EC numbers from reactions in ModelSEEDDatabase
        modelseed_reactions = pd.read_table(
            os.path.join(
                self.config["paths"]["modelseed"],
                "reactions.tsv"
            )
        )

        # Create empty dataframe
        all_ec_numbers_df = pd.DataFrame()

        for organism in metadata["Code"]:

            model_path = os.path.join(
                self.config["paths"]["models"],
                f"{organism}.json"
            )

            # Load model
            with open(model_path, mode="r") as fh:
                model_dict = json.loads(fh.read())

            # Append results to dataframe
            all_ec_numbers_df = pd.concat(
                [
                    all_ec_numbers_df,
                    self.get_ec_from_model(
                        model_dict=model_dict,
                        modelseed_reactions=modelseed_reactions
                    )
                ],
                axis=0,
                ignore_index=True
            )

        # Save to file
        all_ec_numbers_df.to_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                self.config["retropath"]["files"]["ec_numbers"]
            ),
            header=True,
            index=False,
            sep=",",
            mode="w",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
        )

        return all_ec_numbers_df

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

        # Save to file
        rules_df.to_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                self.config["retropath"]["files"]["rules"]
            ),
            header=True,
            index=False,
            sep=",",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
        )

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
