import os

import json

import csv
import pandas as pd

from rdkit import Chem


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

    def get_sink(self) -> pd.DataFrame:
        """
        Get the sink by extracting compounds from the rules when:
            - They appear in the reactants
            - They appear in the products if their rule is reversible

        Parameters
        ----------
        None

        Returns
        -------
        sink_df : DataFrame
            Dataframe containing as the compounds found for the given rules.

        Examples
        --------
        None

        """

        # Read MetaNetX reac_prop.tsv file
        metanetx_reac_prop = pd.read_table(
            os.path.join(
                self.config["paths"]["metanetx"],
                "reac_prop.tsv"
            ),
            comment="#", # Skip comment rows
            header=None,
            names=[
                "ID",
                "mnx_equation",
                "reference",
                "classifs",
                "is_balanced",
                "is_transport"
            ]
        )

        # Read MetaNetX chem_prop.tsv file
        metanetx_chem_prop = pd.read_table(
            os.path.join(
                self.config["paths"]["metanetx"],
                "chem_prop.tsv"
            ),
            comment="#", # Skip comment rows
            header=None,
            names=[
                "Name", # Format for sink.csv
                "compound",
                "reference",
                "formula",
                "charge",
                "mass",
                "InChI",
                "InChIKey",
                "SMILES"
            ]
        )

        # Load rules extracted by mapping ECs to RetroRules DB
        rules_df = pd.read_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                self.config["retropath"]["files"]["rules"]
            ),
            sep=","
        )

        # Get all MetaNetX reaction IDs
        rules_df[["MNXR", "MNXM"]]  = rules_df["Rule ID"]\
            .str.split("_", expand=True)

        # Drop diameter and rule usage duplicates
        unique_rules_df = rules_df\
            .drop_duplicates(subset=["MNXR", "Rule usage"])

        # Get all reaction information by MNXR
        reactions_df = pd.merge(
            left=metanetx_reac_prop,
            right=unique_rules_df,
            left_on="ID",
            right_on="MNXR",
            how="outer",
            indicator=True
        )

        # Remove "left_only" reactions (not present in the organism rules)
        reactions_df = reactions_df[reactions_df["_merge"] != "left_only"]

        # Get all reactions that could not be mapped in METANETX
        not_in_metanetx = reactions_df[reactions_df["_merge"] == "right_only"]
        print(
            "[WARNING] #reactions not present in METANETX: ",
            len(not_in_metanetx), "/", len(reactions_df)
        )

        # Finally, get reactions that could be mapped to METANETX
        reactions_df = reactions_df[reactions_df["_merge"] == "both"]

        # NOTE: METANETX reactions are in forward format. Therefore, we should 
        # get the substrates for the sink.
        compounds_forward = reactions_df["mnx_equation"]\
            .str.split("=").str[0]\
            .str.extractall(r"(MNXM\d+)")\
            .reset_index(drop=True)\
            .rename(columns={0: "Name"})

        # Add products of reversible reactions
        compounds_reverse = reactions_df[
                reactions_df["Rule usage"] == "both"
            ]["mnx_equation"]\
            .str.split("=").str[-1]\
            .str.extractall(r"(MNXM\d+)")\
            .reset_index(drop=True)\
            .rename(columns={0: "Name"})

        # Merge substrates and products (reversible reactions)
        sink_df = pd.concat(
            [compounds_forward, compounds_reverse],
            axis=0,
            ignore_index=True
        )

        # Drop potential duplicates (compounds in more than one reation)
        sink_df = sink_df.drop_duplicates()

        # Get InChIs for the extracted compound IDs
        sink_df = pd.merge(
            left=sink_df,
            right=metanetx_chem_prop,
            on="Name",
            how="left"
        )
        sink_df = sink_df[["Name", "InChI"]]

        # WARNING: there are compounds without InChI
        no_inchi = sink_df[sink_df["InChI"].isnull()]
        print(
            "[WARNING] #compounds without InChI: ",
            len(no_inchi), "/", len(sink_df)
        )

        # Fill missing InChIs and match sink.csv format
        sink_df = sink_df.fillna("None")

        # Save to file
        sink_df.to_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                self.config["retropath"]["files"]["sink"]
            ),
            header=True,
            index=False,
            sep=",",
            mode="w",
            quotechar='"',
            quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
        )

        return sink_df

    def get_sources(self) -> pd.DataFrame:
        """
        Get the sources by creating a folder and CSV file for each compound
        in the specified sources file.

        Parameters
        ----------
        None

        Returns
        -------
        sources_df : DataFrame
            Dataframe containing all the source compounds with their 
            corresponding InChI.

        Examples
        --------
        None

        """

        sources_df = pd.read_csv(
            os.path.join(
                self.config["paths"]["retropath"],
                "interesting_metabolites/",
                self.config["retropath"]["files"]["sources"]
            )
        )

        # Get InChIs
        sources_df["InChI"] = \
            sources_df["Smile"]\
            .dropna(how="all")\
            .apply(lambda row: Chem.MolToInchi(Chem.MolFromSmiles(row)))

        # TODO: log compounds with missing InChI

        # Drop metabolites withou InChI
        sources_df = sources_df.dropna(subset="InChI")

        # Rename to fit RetroPath2.0 format
        sources_df = sources_df\
            .rename(columns={"Molecule": "Name"})

        sources_df["Name"] = sources_df["Name"]\
            .str.lower()\
            .str.replace(" ", "_")

        # Create a source file for each compound
        for _, row in sources_df.iterrows():
            compound_name = row["Name"]

            row = row.to_frame().T[["Name", "InChI"]]

            row.to_csv(
                os.path.join(
                    self.config["paths"]["retropath"],
                    "interesting_metabolites/",
                    "sources/",
                    f"{compound_name}.csv"
                ),
                header=True,
                index=False
            )

        return sources_df
