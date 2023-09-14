import logging

from typing import Iterable

import os
import json

import csv
import pandas as pd

from rdkit import Chem

from biofoundry.base import BaseRetroPathPreloader


# Configure logging
logging.basicConfig(
    filename="retropath-" + os.path.basename(__file__).replace(".py", ".log"),
    filemode="w",
    format="%(asctime)s - %(filename)s:%(lineno)s - %(funcName)s - " + \
        "%(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.DEBUG
)
LOGGER = logging.getLogger("retropath-" + __name__)


class RetroPathPreloader(BaseRetroPathPreloader):
    """
    Auxiliary class for generating RetroPath2.0 inputs.

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
        modelseed_reactions : pandas.DataFrame
            ModelSEED reactions mapping model IDs (RXNs) to EC numbers.

        Returns
        -------
        ec_numbers_df : pandas.DataFrame
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

        LOGGER.info(f"Number of reaction IDs: {len(reaction_ids)}")

        ec_numbers_df = modelseed_reactions[
            modelseed_reactions["id"].isin(reaction_ids)
        ]["ec_numbers"]

        LOGGER.info(f"Number of EC numbers (raw): {len(ec_numbers_df)}")

        # Explode series since there may be multiple EC numbers per reaction
        ec_numbers_df = ec_numbers_df.str.split("|").explode()

        LOGGER.info(f"Number of EC numbers: {len(ec_numbers_df)}")

        # Convert to frame
        ec_numbers_df = ec_numbers_df.to_frame()

        # Add model ID
        ec_numbers_df["ID"] = model_dict["id"]

        return ec_numbers_df

    def get_ec_numbers(self, metadata: pd.DataFrame) -> pd.DataFrame:
        """
        Get all the EC numbers for each model specified in the metadata.

        Parameters
        ----------
        metadata : pandas.DataFrame
            The dataframe containing the models' paths.

        Returns
        -------
        all_ec_numbers_df : pandas.DataFrame
            Dataframe containing each all EC numbers found for the model.

        Examples
        --------
        None

        """

        # Get EC numbers from reactions in ModelSEED
        modelseed_path = os.path.join(
            self.config["paths"]["modelseed"],
            "reactions.tsv"
        )
        modelseed_reactions = pd.read_table(modelseed_path)

        LOGGER.debug(f"Loaded ModelSEED from {modelseed_path}")
        LOGGER.debug(
             f"Number of ModelSEED reactions: {len(modelseed_reactions)}"
        )

        # Create empty dataframe
        all_ec_numbers_df = pd.DataFrame()

        for organism in metadata["Code"]:

            LOGGER.info(f"Starting with organism {organism}")

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
        ec_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["ec_numbers"]
        )
        all_ec_numbers_df.to_csv(
            ec_path,
            header=True,
            index=False,
            sep=",",
            mode="w",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
        )

        LOGGER.info(f"Saved community EC numbers to {ec_path}")
        LOGGER.info(f"Total amount of EC numbers: {len(all_ec_numbers_df)}")

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
        rules_df : pandas.DataFrame
            Dataframe containing all rules found in the community.

        Examples
        --------
        None

        """

        # Load EC numbers in the community
        ec_numbers_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["ec_numbers"]
        )
        ec_num_df = pd.read_csv(ec_numbers_path, sep=",")

        LOGGER.debug(f"Loaded community EC numbers from {ec_numbers_path}")
        LOGGER.debug(f"Number of ECs: {len(ec_num_df)}")

        # Load RetroRules database
        retrorules_path = self.config["paths"]["retrorules"]
        retrorules_df = pd.read_csv(retrorules_path, sep=",")

        LOGGER.debug(f"Loaded RetroRules from {retrorules_path}")
        LOGGER.debug(f"Number of RetroRules: {len(retrorules_df)}")

        # Drop potential duplicates (more than one species with the same EC)
        ec_numbers = ec_num_df["ec_numbers"].unique()

        LOGGER.debug(f"Dropped duplicates in ECs: {len(ec_numbers)} unique")

        # Explode rules with multiple EC numbers
        retrorules_df["EC number"] = retrorules_df["EC number"].str.split(";")
        retrorules_df = retrorules_df\
            .explode("EC number")\
            .reset_index(drop=True)

        LOGGER.debug(f"Exploded ECs in RetroRules: {len(retrorules_df)} total")

        # Filter by found EC numbers in the community
        rules_df = retrorules_df[
            retrorules_df["EC number"].isin(ec_numbers)
        ]

        LOGGER.info(
            f"Retrieved rules by ECs present in the community: {len(rules_df)}"
        )

        # Save to file
        rules_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["rules"]
        )
        rules_df.to_csv(
            rules_path,
            header=True,
            index=False,
            sep=",",
            quotechar='"',
            quoting=csv.QUOTE_MINIMAL # Avoid errors with commas in names
        )

        LOGGER.info(f"Saved community rules to {rules_path}")

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
        sink_df : pandas.DataFrame
            Dataframe containing as the compounds found for the given rules.

        Examples
        --------
        None

        """

        # Read MetaNetX reac_prop.tsv file
        metanetx_reac_path = os.path.join(
            self.config["paths"]["metanetx"],
            "reac_prop.tsv"
        )
        metanetx_reac_prop = pd.read_table(
            metanetx_reac_path,
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

        LOGGER.debug(f"Loaded MetaNetX reactions from {metanetx_reac_path}")
        LOGGER.debug(f"Number of MetaNetX reactions: {len(metanetx_reac_prop)}")

        # Read MetaNetX chem_prop.tsv file
        metanetx_chem_path = os.path.join(
            self.config["paths"]["metanetx"],
            "chem_prop.tsv"
        )
        metanetx_chem_prop = pd.read_table(
            metanetx_chem_path,
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

        LOGGER.debug(f"Loaded MetaNetX compounds from {metanetx_chem_path}")
        LOGGER.debug(f"Number of MetaNetX compounds: {len(metanetx_chem_prop)}")

        # Load rules extracted by mapping ECs to RetroRules DB
        rules_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["rules"]
        )
        rules_df = pd.read_csv(rules_path, sep=",")

        LOGGER.debug(f"Loaded community rules from {rules_path}")
        LOGGER.debug(f"Number of community rules: {len(rules_df)}")

        # Get all MetaNetX reaction IDs
        rules_df[["MNXR", "MNXM"]]  = rules_df["Rule ID"]\
            .str.split("_", expand=True)

        # Drop diameter and rule usage duplicates
        unique_rules_df = rules_df\
            .drop_duplicates(subset=["MNXR", "Rule usage"])

        LOGGER.debug(
            "Dropped duplicated rules (diameter and usage columns): " + \
            str(len(unique_rules_df))
        )

        # Get all reaction information by MNXR
        reactions_df = pd.merge(
            left=metanetx_reac_prop,
            right=unique_rules_df,
            left_on="ID",
            right_on="MNXR",
            how="outer",
            indicator=True
        )

        LOGGER.debug(
            "Number of reactions only present in MetaNetX: " + \
            str(len(reactions_df[reactions_df["_merge"] == "right_only"]))
        )
        LOGGER.debug(
            "Number of reactions only present in the community: " + \
            str(len(reactions_df[reactions_df["_merge"] == "left_only"]))
        )
        LOGGER.debug(
            "Number of reactions shared: " + \
            str(len(reactions_df[reactions_df["_merge"] == "both"]))
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

        LOGGER.debug(
            "Number of compounds in forward rules: " + \
            str(len(compounds_forward))
        )

        # Add products of reversible reactions
        compounds_reverse = reactions_df[
                reactions_df["Rule usage"] == "both"
            ]["mnx_equation"]\
            .str.split("=").str[-1]\
            .str.extractall(r"(MNXM\d+)")\
            .reset_index(drop=True)\
            .rename(columns={0: "Name"})

        LOGGER.debug(
            "Number of compounds in reversible rules: " + \
            str(len(compounds_reverse))
        )

        # Merge substrates and products (reversible reactions)
        sink_df = pd.concat(
            [compounds_forward, compounds_reverse],
            axis=0,
            ignore_index=True
        )

        LOGGER.debug(f"Number of compounds in sink (raw): {len(sink_df)}")

        # Drop potential duplicates (compounds in more than one reation)
        sink_df = sink_df.drop_duplicates()

        LOGGER.debug(f"Dropped duplicates in sink (raw): {len(sink_df)}")

        # Get InChIs for the extracted compound IDs
        sink_df = pd.merge(
            left=sink_df,
            right=metanetx_chem_prop,
            on="Name",
            how="left"
        )
        sink_df = sink_df[["Name", "InChI"]]

        LOGGER.warning(
            "Number of compounds without InChI: " + \
            str(len(sink_df[sink_df["InChI"].isnull()])) + "/" + \
            str(len(sink_df))
        )

        # Fill missing InChIs and match sink.csv format
        sink_df = sink_df.fillna("None")

        # Save to file
        sink_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["sink"]
        )
        sink_df.to_csv(
            sink_path,
            header=True,
            index=False,
            sep=",",
            mode="w",
            quotechar='"',
            quoting=csv.QUOTE_ALL # Fit sink.csv format for RetroPath2.0
        )

        LOGGER.info(f"Saved community sink to {sink_path}")

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
        sources_df : pandas.DataFrame
            Dataframe containing all the source compounds with their 
            corresponding InChI.

        Examples
        --------
        None

        """

        sources_path = os.path.join(
            self.config["paths"]["retropath"],
            "interesting_metabolites/",
            self.config["retropath"]["files"]["sources"]
        )
        sources_df = pd.read_csv(sources_path)

        LOGGER.debug(f"Loaded sources from {sources_path}")

        # Get InChIs
        sources_df["InChI"] = \
            sources_df["Smile"]\
            .dropna(how="all")\
            .apply(
                lambda row: Chem.MolToInchi(
                    Chem.MolFromSmiles(row),
                    options="-SNon"
                )
            )

        LOGGER.warning(
            "Number of compounds without InChI: " + \
            str(len(sources_df[sources_df["InChI"].isnull()])) + "/" + \
            str(len(sources_df))
        )

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

            compound_path = os.path.join(
                self.config["paths"]["retropath"],
                "interesting_metabolites/",
                "sources/",
                f"{compound_name}.csv"
            )
            row.to_csv(
                compound_path,
                header=True,
                index=False
            )

            LOGGER.info(f"Saved compound to {compound_path}")

        return sources_df

    def get_orgs_with_source_in_sink(
        self,
        source_ids: Iterable[str]
    ) -> pd.DataFrame:
        """
        Get the organisms containing the source already in the sink.

        Parameters
        ----------
        source_ids : Iterable[str]
            The list of possible source's IDs (e.g. enantiomers have different
            IDs).

        Returns
        -------
        reactions_df : pandas.DataFrame
            Dataframe containing the organism code and the found ID.

        Examples
        --------
        None

        """

        source_str = r"|".join(source_ids)

        # Load EC numbers and rules in the community
        ec_numbers_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["ec_numbers"]
        )
        ec_num_df = pd.read_csv(ec_numbers_path, sep=",")

        LOGGER.debug(f"Loaded community EC numbers from {ec_numbers_path}")
        LOGGER.debug(f"Number of ECs: {len(ec_num_df)}")

        rules_path = os.path.join(
            self.config["paths"]["retropath"],
            self.config["retropath"]["files"]["rules"]
        )
        rules_df = pd.read_csv(rules_path, sep=",")

        LOGGER.debug(f"Loaded community rules from {rules_path}")
        LOGGER.debug(f"Number of community rules: {len(rules_df)}")

        # Drop potential duplicates (more than one species with the same EC)
        ec_numbers = ec_num_df["ec_numbers"].unique()

        LOGGER.debug(f"Dropped duplicates in ECs: {len(ec_numbers)} unique")

        # Avoid duplicates because of different rule usages and diameters
        rules_df = rules_df.drop(["Diameter", "Rule usage"], axis=1)

        LOGGER.debug(
            "Dropped duplicated rules (diameter and usage columns): " + \
            str(len(rules_df))
        )

        # Add organism to rules
        rules_df = pd.merge(
            left=rules_df,
            right=ec_num_df,
            left_on="EC number",
            right_on="ec_numbers",
            how="left"
        )
        rules_df = rules_df.rename(columns={"ID": "Organism"})

        # Read MetaNetX reac_prop.tsv file
        metanetx_reac_path = os.path.join(
            self.config["paths"]["metanetx"],
            "reac_prop.tsv"
        )
        metanetx_reac_prop = pd.read_table(
            metanetx_reac_path,
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

        LOGGER.debug(f"Loaded MetaNetX reactions from {metanetx_reac_path}")
        LOGGER.debug(f"Number of MetaNetX reactions: {len(metanetx_reac_prop)}")

        # Get all MetaNetX reaction IDs
        rules_df[["MNXR", "MNXM"]]  = rules_df["Rule ID"]\
            .str.split("_", expand=True)

        # Get all reaction information by MNXR
        reactions_df = pd.merge(
            left=metanetx_reac_prop,
            right=rules_df,
            left_on="ID",
            right_on="MNXR",
            how="outer",
            indicator=True
        )

        LOGGER.debug(
            "Number of reactions only present in MetaNetX: " + \
            str(len(reactions_df[reactions_df["_merge"] == "right_only"]))
        )
        LOGGER.debug(
            "Number of reactions only present in the community: " + \
            str(len(reactions_df[reactions_df["_merge"] == "left_only"]))
        )
        LOGGER.debug(
            "Number of reactions shared: " + \
            str(len(reactions_df[reactions_df["_merge"] == "both"]))
        )

        # Finally, get reactions that could be mapped to METANETX
        reactions_df = reactions_df[reactions_df["_merge"] == "both"]
        reactions_df["Source match"] = reactions_df["mnx_equation"]\
            .str.findall(source_str)\
            .apply(",".join)

        # Get only matches
        reactions_df = reactions_df[
            reactions_df["Source match"] != ""
        ]

        reactions_df = reactions_df[["Organism", "Source match"]]\
            .drop_duplicates()

        LOGGER.debug(f"Dropped duplicates in matches: {len(reactions_df)}")

        return reactions_df
