import os

import pandas as pd

from tabula import read_pdf

from biofoundry.base import BaseDataLoader


class Amils2023DataLoader(BaseDataLoader):

    def __init__(
        self,
        data_dir: str = "../data/papers/amils2023/"
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
        data_df_long : pandas.DataFrame
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
        data_df_long : pandas.DataFrame
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
        data_df_long : pandas.DataFrame
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
        data_df_long : pandas.DataFrame
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
        Get from all previous sources combined.

        Parameters
        ----------
        None

        Returns
        -------
        data_df_long : pandas.DataFrame
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

    def get_microbial_data(self) -> pd.DataFrame:
        """
        Get the microbial data from table S8 in Amils et al. 2023.

        Parameters
        ----------
        None

        Returns
        -------
        microbes_df : pandas.DataFrame
            Dataframe containing the microbial data from Amils et al. 2023.

        Examples
        --------
        >>> from biofoundry.data.amils2023 import Amils2023DataLoader
        >>> data_loader = Amils2023DataLoader()
        >>> microbes_df = data_loader.get_microbial_data()

        """

        microbes_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0001-supinfo-tables8-2.ods"
            ),
            sheet_name="Sheet1"
        )

        # Rename pathway column
        microbes_df = microbes_df.rename(columns={"Pathway/depth": "Pathway"})

        # Drop last row containing the explanation
        microbes_df = microbes_df.iloc[:-1, :].copy()

        # Drop rows containing the cycles
        microbes_df = microbes_df[
            ~microbes_df["Pathway"].str.endswith(" cycle")
        ]

        # Convert to numeric
        numeric_cols = [
            col for col in microbes_df.columns
            if col not in ["Pathway"]
        ]
        microbes_df[numeric_cols] = microbes_df[numeric_cols].apply(
            pd.to_numeric,
            errors="coerce"
        )

        return microbes_df

    def get_abundances(self) -> pd.DataFrame:
        """
        Get the abundance data from Amils et al. 2023.

        Parameters
        ----------
        None

        Returns
        -------
        abundances_df : pandas.DataFrame
            Dataframe containing the microbial data from Amils et al. 2023.

        Notes
        -----
        Abundances come from the two methods used in the paper: Illumin and
        Roche 454.

        Examples
        --------
        >>> from biofoundry.data.amils2023 import Amils2023DataLoader
        >>> data_loader = Amils2023DataLoader()
        >>> abundance_df = data_loader.get_abundance_data()

        """

        # Read 454 and Illumina abundances from supplementary materials
        abundances_illumina_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0005-datasets4.xlsx"
            ),
            sheet_name="Filtered OTUs",
            skiprows=11
        )
        abundances_roche_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0006-datasets5.xlsx"
            ),
            sheet_name="DW Filtered OTUs",
            skiprows=11
        )

        # Filter out species also present in the drilling water (possible contamination)
        # abundances_illumina_df = abundances_illumina_df[
        #     (abundances_illumina_df["DW_RG"] == 0) & \
        #     (abundances_illumina_df["IC"] == 0)
        # ]
        # abundances_roche_df = abundances_roche_df[
        #     abundances_roche_df["DWδ"] == 0
        # ]

        # Group samples by genus to avoid repeats
        abundances_illumina_df = abundances_illumina_df\
            .groupby("Genus", as_index=False)\
            .sum()
        abundances_roche_df = abundances_roche_df\
            .groupby("Genus", as_index=False)\
            .sum()

        # Get samples in long format for Illumina reads
        abundances_illumina_df = pd.wide_to_long(
            df=abundances_illumina_df,
            stubnames="BH10",
            i="Genus",
            j="Sample",
            sep="-",
            suffix="\\d+"
        )
        abundances_illumina_df = abundances_illumina_df\
            .reset_index()\
            [["Genus", "Sample", "BH10"]]\
            .rename(columns={
                "BH10": "abundance",
                "Genus": "genus",
                "Sample": "sample_id"
            })

        abundances_illumina_df["sample_id"] = abundances_illumina_df["sample_id"]\
            .apply(lambda row: f"BH10-{str(row)}-Illumina")\
            .astype(str)

        # Get samples in long format for Roche reads
        abundances_roche_df = pd.wide_to_long(
            df=abundances_roche_df,
            stubnames="BH10",
            i="Genus",
            j="Sample",
            sep="-",
            suffix="\\d+"
        )
        abundances_roche_df = abundances_roche_df\
            .reset_index()\
            [["Genus", "Sample", "BH10"]]\
            .rename(columns={
                "BH10": "abundance",
                "Genus": "genus",
                "Sample": "sample_id"
            })

        abundances_roche_df["sample_id"] = abundances_roche_df["sample_id"]\
            .apply(lambda row: f"BH10-{str(row)}-Roche")

        # Concatenate
        abundances_df = pd.concat(
            [abundances_illumina_df, abundances_roche_df],
            axis=0,
            ignore_index=True
        )

        return abundances_df
