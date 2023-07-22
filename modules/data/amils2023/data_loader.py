import os

import pandas as pd

from tabula import read_pdf

from modules.base import BaseDataLoader


class Amils2023DataLoader(BaseDataLoader):

    def __init__(
        self,
        data_dir: str = "../data/papers/amils2023"
    ) -> None:
        super().__init__()

        self.data_dir = data_dir

    def get_elements(self) -> pd.DataFrame:
        """
        Get data from Dataset S2 - ICP-MS elemental analysis of core samples 
        (ppm).

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : DataFrame
            Dataframe in long format containing the data.

        Examples
        --------
        None

        """

        data_df = read_pdf(
            input_path=os.path.join(
                self.data_dir,
                "emi16291-sup-0003-datasets2.pdf"
            ),
            pages=1
        )

        # Extract the only dataframe since there's only one page
        data_df = data_df[0]

        # Convert to long for plotting
        data_df_long = pd.melt(
            data_df,
            id_vars=["Depth"],
            value_vars=[
                col for col in data_df.columns
                if col not in ["Depth"]
            ],
            var_name="Species",
            value_name="Concentration (ppm)"
        )

        # Convert numeric columns to float
        data_df_long[["Depth", "Concentration (ppm)"]] = \
            data_df_long[["Depth", "Concentration (ppm)"]]\
            .apply(lambda row: row.str.replace(",", "."))\
            .astype(float)

        return data_df_long

    def get_anions(self) -> pd.DataFrame:
        """
        Get data from Dataset S3 - Ionic chromatography of BH10 soluble organic 
        and inorganic anions (ppm).

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : DataFrame
            Dataframe in long format containing the data.

        Examples
        --------
        None

        """

        data_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0004-datasets3.xls"
            ),
            sheet_name="BH10_CI",
            skiprows=1
        )

        # Drop W90 sample to avoid duplicates
        data_df = data_df[data_df["SAMPLE"] != "BH10_W90"]\
            .reset_index(drop=True)

        # Create depth column
        data_df["Depth"] = data_df["SAMPLE"]\
            .str.split("_").str[-1]\
            .str.split("-").str[-1]

        data_df = data_df.drop("SAMPLE", axis=1)

        # Remove letters from numbers (e.g. depth W90)
        data_df["Depth"] = data_df["Depth"]\
            .str.strip(r"\W")

        # Convert to numeric
        data_df["Depth"] = data_df["Depth"]\
            .str.replace(",", ".")\
            .astype(float)

        data_df_long = pd.melt(
            data_df,
            id_vars=["Depth"],
            value_vars=[
                col for col in data_df.columns
                if col not in ["Depth"]
            ],
            var_name="Species",
            value_name="Concentration (ppm)"
        )

        # Capitalize to fit format
        data_df_long["Species"] = data_df_long["Species"]\
            .str.capitalize()

        # Replace capitalized pH
        data_df_long["Species"] = data_df_long["Species"]\
            .str.replace("Ph(?!\w+)", "pH", regex=True)

        return data_df_long

    def get_cations(self) -> pd.DataFrame:
        """
        Get data from Table S1 - Soluble cations (ppm).

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : DataFrame
            Dataframe in long format containing the data.

        Examples
        --------
        None

        """

        data_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0001-supinfo-tables1.ods"
            ),
            sheet_name="Sheet1",
            skiprows=1
        )

        data_df_long = pd.melt(
            data_df,
            id_vars=["Depth"],
            value_vars=[
                col for col in data_df.columns
                if col not in ["Depth"]
            ],
            var_name="Species",
            value_name="Concentration (ppm)"
        )

        return data_df_long

    def get_gases(self) -> pd.DataFrame:
        """
        Get data from Table S7 - Occluded gases and natural activities at 
        different depths (10/8/22).

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : DataFrame
            Dataframe in long format containing the data.

        Examples
        --------
        None

        """

        data_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0001-supinfo-tables7.ods"
            ),
            sheet_name="Sheet1",
            skiprows=1
        )

        # Drop last row containing the explanation
        data_df = data_df.iloc[:-1, :].copy()

        # Drop last three columns since they correspond to activate metabolism
        data_df = data_df.iloc[:, :-3].copy()

        # Rename depth column
        data_df = data_df.rename(columns={"Depth mbs": "Depth"})

        # Apply conversion of values
        h2_co2_symbol_map = {
            "+++": 2000,
            "++": 1000,
            "+": 250,
            "+/-": 40,
            "-": 0
        }
        ch4_symbol_map = {
            "+++": 35,
            "++": 20,
            "+": 10,
            "+/-": 5,
            "-": 0
        }

        data_df["H2"] = data_df["H2"].replace(h2_co2_symbol_map).astype(float)
        data_df["CO2"] = data_df["CO2"].replace(h2_co2_symbol_map).astype(float)
        data_df["CH4"] = data_df["CH4"].replace(ch4_symbol_map).astype(float)

        data_df_long = pd.melt(
            data_df,
            id_vars=["Depth"],
            value_vars=[
                col for col in data_df.columns
                if col not in ["Depth"]
            ],
            var_name="Species",
            value_name="Concentration (ppm)"
        )

        return data_df_long

    def get_data(self) -> pd.DataFrame:
        """
        Get data from Table S7 - Occluded gases and natural activities at 
        different depths (10/8/22).

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : DataFrame
            Dataframe in long format containing the data.

        Examples
        --------
        None

        """

        # Concatenate all data
        data_df_long = pd.concat(
            [
                self.get_elements(),
                self.get_anions(),
                self.get_cations(),
                self.get_gases()
            ],
            axis=0,
            ignore_index=True
        )

        # Round to fit Illumina and Roche datasets (Datasets S4 and S5)
        data_df_long["Depth"] = data_df_long["Depth"].astype(int)

        return data_df_long
