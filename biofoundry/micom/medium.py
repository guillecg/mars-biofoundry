import os
import pandas as pd

from micom.media import Community, minimal_medium

from biofoundry.base import BaseMICOMMediumManager


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
        min_medium : pd.DataFrame
            The computed minimal medium.

        Examples
        --------
        None

        """

        sol = com.cooperative_tradeoff()

        # Extracellular medium has no growth rate
        rates = sol.members.growth_rate.drop("medium")

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

        return min_medium
