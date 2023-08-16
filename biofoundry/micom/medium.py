import logging

import os
import string

import pandas as pd

from micom.media import Community, minimal_medium

from biofoundry.base import BaseMICOMMediumManager


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


class MICOMMediumManager(BaseMICOMMediumManager):
    """
    Auxiliary class for creating the medium for MICOM.

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
    def get_min_medium(
        com: Community,
        growth_tolerance: float = 0.75
    ) -> pd.DataFrame:
        """
        Get minimal medium for the given community.

        Parameters
        ----------
        com : Community
            The selected MICOM community.
        growth_tolerance : float
            Fraction of both community and minimum growth as defined in
            MICOM's minimal_medium function.

        Returns
        -------
        min_medium : pandas.DataFrame
            The computed minimal medium.

        Examples
        --------
        None

        """

        LOGGER.info("Starting to solve the cooperative tradeoff...")

        sol = com.cooperative_tradeoff()

        # Extracellular medium has no growth rate
        rates = sol.members.growth_rate.drop("medium")

        LOGGER.info("Starting to get the minimal medium...")

        # Get minimal medium
        min_medium = minimal_medium(
            community=com,
            community_growth=growth_tolerance*sol.growth_rate,
            min_growth=growth_tolerance*rates,
            minimize_components=True
        )

        min_medium = min_medium\
            .sort_values()\
            .to_frame("flux")\
            .reset_index(names="reaction")

        LOGGER.debug(
            f"Number of fluxes in the minimal medium: {len(min_medium)}"
        )

        return min_medium


    def extrapolate_medium(
        self,
        medium_df: pd.DataFrame,
        min_medium: pd.DataFrame,
        shared_comp: str,
        shared_flux: str
    ) -> pd.DataFrame:
        """
        Extrapolate fluxes from raw medium concentrations using minimal medium
        computed fluxes.

        Parameters
        ----------
        medium_df : pandas.DataFrame
            The concentration of metabolites in the environment.
        min_medium : pandas.DataFrame
            The computed minimal medium.
        shared_comp : str
            The compound shared by both dataframes which will be used for the
            extrapolation.
        shared_flux : str
            The flux shared by both dataframes which will be used for the
            extrapolation.

        Returns
        -------
        medium_adjusted : pd.DataFrame
            The extrapolated fluxes.

        Examples
        --------
        None

        """

        LOGGER.info(f"Compound name of shared metabolite: {shared_comp}")
        LOGGER.info(f"Flux name of shared metabolite: {shared_flux}")

        # Get ratio for mapping concentrations to fluxes
        shared_conc = medium_df.loc[
            (medium_df["Species"] == shared_comp), "Concentration (ppm)"
        ].values[0]

        shared_conc_flux = min_medium.loc[
            (min_medium["reaction"] == shared_flux), "flux"
        ].values[0]

        LOGGER.debug(f"Concentration of shared metabolite: {shared_conc}")
        LOGGER.debug(f"Flux of shared metabolite: {shared_conc_flux}")
        LOGGER.info(
            "Ratio of flux to concentration: " + \
            str({shared_conc_flux / shared_conc})
        )

        medium_df["flux"] = \
            medium_df["Concentration (ppm)"] * shared_conc_flux / shared_conc

        # -------------------------------------------------------------------- #

        modelseed_path = os.path.join(
            self.config["paths"]["modelseed"],
            "compounds.tsv"
        )
        modelseed_cpd = pd.read_table(
            modelseed_path,
            dtype=object # Avoid warnings
        )

        LOGGER.debug(f"Loaded ModelSEED from {modelseed_path}")

        # Try to get the maximum of chemical species
        modelseed_map = modelseed_cpd[
            (modelseed_cpd["name"].isin(medium_df["Species"])) |
            (modelseed_cpd["formula"].isin(medium_df["Species"])) |
            (modelseed_cpd["abbreviation"].isin(medium_df["Species"]))
        ]

        # Get only primary databases sources
        modelseed_map = modelseed_map[
            modelseed_map["source"] == "Primary Database"
        ].copy()

        LOGGER.warning(
            f"Number of retrieved chemical species from ModelSEED: " + \
            f"{len(modelseed_map)}/{len(medium_df['Species'])}"
        )

        # Avoid errors in COBRApy due to punctuation and whitespaces in names
        modelseed_map["abbreviation"] = modelseed_map["abbreviation"]\
            .str.replace(
                pat=r"[{}|\s]".format(string.punctuation),
                repl="-",
                regex=True
            )

        # Create new column to avoid abbreviation duplicates
        modelseed_map["abbreviation_id"] = \
            "EX_" + modelseed_map["id"] + "=" + \
            modelseed_map["abbreviation"] + "_m"

        # Filter columns
        modelseed_map = modelseed_map[[
            "id",
            "name",
            "formula",
            "abbreviation",
            "abbreviation_id"
        ]]

        # Wide to long
        modelseed_map = pd.melt(
            modelseed_map,
            id_vars=["abbreviation_id"],
            value_vars=["name", "formula", "abbreviation"],
            value_name="Species"
        )
        modelseed_map = modelseed_map[["abbreviation_id", "Species"]]\
            .drop_duplicates()

        LOGGER.debug(f"Dropped duplicates in mapping: {len(modelseed_map)}")

        # Add ModelSEED mapping for matching fluxes in the models 
        medium_mapped = pd.merge(
            left=medium_df,
            right=modelseed_map,
            on="Species",
            how="inner"
        )

        # Fix format
        medium_mapped = medium_mapped.rename(columns={
            "Flux": "flux",
            "abbreviation_id": "reaction"
        })
        medium_mapped = medium_mapped[["reaction", "flux"]]

        # Add extrapolated fluxes to minimal media
        medium_adjusted = pd.concat(
            [medium_mapped, min_medium],
            axis=0,
            ignore_index=True
        )

        # Prioritize values for mapped fluxes instead of original ones
        medium_adjusted = medium_adjusted.drop_duplicates(keep="first")

        LOGGER.info(
            f"Number of fluxes in the adjusted medium: {len(medium_adjusted)}"
        )

        medium_adjusted = dict(zip(
            medium_adjusted["reaction"],
            medium_adjusted["flux"]
        ))

        return medium_adjusted
