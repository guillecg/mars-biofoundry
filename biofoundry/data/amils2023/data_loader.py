import logging

# Ignore openpyxl warnings: https://stackoverflow.com/a/75025242
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

import os

import pandas as pd

from tabula import read_pdf

from biofoundry.base import BaseDataLoader


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

        LOGGER.info("Loaded emi16291-sup-0003-datasets2.pdf")

        # Extract the only dataframe since there's only one page
        data_df = data_df[0]

        LOGGER.info(f"Number of entries (wide): {data_df.shape}")

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

        LOGGER.info(f"Number of entries (long): {data_df_long.shape}")

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

        LOGGER.info("Loaded emi16291-sup-0004-datasets3.xls")
        LOGGER.info(f"Number of entries (wide): {data_df.shape}")

        # Drop W90 sample to avoid duplicates
        data_df = data_df[data_df["SAMPLE"] != "BH10_W90"]\
            .reset_index(drop=True)
        LOGGER.debug(f"Dropped W90 sample: {data_df.shape}")

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

        LOGGER.info(f"Number of entries (long): {data_df_long.shape}")

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

        LOGGER.info("Loaded emi16291-sup-0001-supinfo-tables1.ods")
        LOGGER.info(f"Number of entries (wide): {data_df.shape}")

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

        LOGGER.info(f"Number of entries (long): {data_df_long.shape}")

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

        LOGGER.info("Loaded emi16291-sup-0001-supinfo-tables7.ods")
        LOGGER.info(f"Number of entries (wide): {data_df.shape}")

        # Drop last row containing the explanation
        data_df = data_df.iloc[:-1, :].copy()
        LOGGER.debug(f"Dropped last row: {data_df.shape}")

        # Drop last three columns since they correspond to activate metabolism
        data_df = data_df.iloc[:, :-3].copy()
        LOGGER.debug(f"Dropped last three columns: {data_df.shape}")

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

        LOGGER.info(f"Number of entries (long): {data_df_long.shape}")

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

        LOGGER.info(f"Number of entries (long): {data_df_long.shape}")

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

        LOGGER.info("Loaded emi16291-sup-0001-supinfo-tables8-2.ods")
        LOGGER.info(f"Number of entries (wide): {microbes_df.shape}")

        # Rename pathway column
        microbes_df = microbes_df.rename(columns={"Pathway/depth": "Pathway"})

        # Fix super- and subscripts
        microbes_df["Pathway"] = microbes_df["Pathway"]\
            .str.replace(
                "comp denitr",
                "Complete denitrification"
            )\
            .str.replace(
                "N2 fix",
                r"$N_{2} \; fix$",
                regex=True
            )\
            .str.replace(
                "SO32- → S2-",
                r"$SO_{3}^{2-} \\xrightarrow{} S^{2-}$",
                regex=True
            )\
            .str.replace(
                "S4O62- →  S2O32-",
                r"$S_{4}O_{6}^{2-} \\xrightarrow{} S_{2}O_{3}^{2-}$",
                regex=True
            )\
            .str.replace(
                "S2O32-→ S2-",
                r"$S_{2}O_{3}^{2-} \\xrightarrow{} S^{2-}$",
                regex=True
            )\
            .str.replace(
                "S2 → pS",
                r"$S^{2-} \\xrightarrow{} pS$",
                regex=True
            )\
            .str.replace(
                "S2O32- → S4O62-",
                r"$S_{2}O_{3}^{2-} \\xrightarrow{} S_{4}O_{6}^{2-}$",
                regex=True
            )\
            .str.replace(
                "S2O32- → SO42-",
                r"$S_{2}O_{3}^{2-} \\xrightarrow{} SO_{4}^{2-}$",
                regex=True
            )\
            .str.replace(
                "H2 ox",
                r"$H_{2} \; oxidation$",
                regex=True
            )\
            .str.replace(
                "Fe3\+ →  Fe2\+",
                r"$Fe^{3+} \\xrightarrow{} Fe^{2+}$",
                regex=True
            )\
            .str.replace(
                "TCA CO2",
                r"$TCA \; CO_{2}$",
                regex=True
            )\
            .str.replace(
                "fm CO2",
                r"$fm \; CO_{2}$",
                regex=True
            )\
            .str.replace(
                "fm H2",
                r"$fm \; H_{2}$",
                regex=True
            )\
            .str.replace(
                "CO2fix",
                r"$CO_{2} \; fix$",
                regex=True
            )\
            .str.replace(
                "CO → CO2",
                r"$CO \\xrightarrow{} CO_{2}$",
                regex=True
            )\
            .str.replace(
                "nº compl cyc",
                "Nº complete cycles"
            )

        LOGGER.debug("Formatted super- and subscripts:")
        LOGGER.debug("\t" + microbes_df.to_string().replace("\n", "\n\t"))

        # Drop last row containing the explanation
        microbes_df = microbes_df.iloc[:-1, :].copy()
        LOGGER.debug(f"Dropped last row: {microbes_df.shape}")

        # Drop rows containing the cycles
        microbes_df = microbes_df[
            ~microbes_df["Pathway"].str.endswith(" cycle")
        ]
        LOGGER.debug(f"Dropped cycles rows: {microbes_df.shape}")

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

        LOGGER.info("Loaded emi16291-sup-0005-datasets4.xlsx")
        LOGGER.info(
            "Number of entries (wide) - Illumina: "  + \
            f"{abundances_illumina_df.shape}"
        )

        abundances_roche_df = pd.read_excel(
            os.path.join(
                self.data_dir,
                "emi16291-sup-0006-datasets5.xlsx"
            ),
            sheet_name="DW Filtered OTUs",
            skiprows=11
        )

        LOGGER.info("Loaded emi16291-sup-0006-datasets5.xlsx")
        LOGGER.info(
            "Number of entries (wide) - Roche: "  + \
            f"{abundances_roche_df.shape}"
        )

        # Filter out species also present in the drilling water (possible contamination)
        # abundances_illumina_df = abundances_illumina_df[
        #     (abundances_illumina_df["DW_RG"] == 0) & \
        #     (abundances_illumina_df["IC"] == 0)
        # ]
        # abundances_roche_df = abundances_roche_df[
        #     abundances_roche_df["DWδ"] == 0
        # ]

        # Get only required columns (genus and BH10 samples)
        bh_cols_illumina = abundances_illumina_df\
            .filter(regex=r"BH10-\d+")\
            .columns\
            .tolist()
        abundances_illumina_df = abundances_illumina_df[
            ["Genus"] + bh_cols_illumina
        ].copy()

        LOGGER.debug(
            "Filtered genus and sample columns - Illumina: " + \
            f"{abundances_illumina_df.columns}"
        )

        bh_cols_roche = abundances_roche_df\
            .filter(regex=r"BH10-\d+")\
            .columns\
            .tolist()
        abundances_roche_df = abundances_roche_df[
            ["Genus"] + bh_cols_roche
        ].copy()

        LOGGER.debug(
            "Filtered genus and sample columns - Roche: " + \
            f"{abundances_roche_df.columns}"
        )

        # Group samples by genus to avoid repeats
        abundances_illumina_df = abundances_illumina_df\
            .groupby("Genus", as_index=False)\
            .sum()

        LOGGER.debug("Grouped samples by genus - Illumina:")
        LOGGER.debug(
            "\t" + abundances_illumina_df.to_string().replace("\n", "\n\t")
        )

        abundances_roche_df = abundances_roche_df\
            .groupby("Genus", as_index=False)\
            .sum()

        LOGGER.debug("Grouped samples by genus - Roche:")
        LOGGER.debug(
            "\t" + abundances_roche_df.to_string().replace("\n", "\n\t")
        )

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

        abundances_illumina_df["sample_id"] = abundances_illumina_df\
            ["sample_id"]\
            .apply(lambda row: f"BH10-{str(row)}-Illumina")\
            .astype(str)

        LOGGER.info(
            "Number of entries (long) - Illumina: "  + \
            f"{abundances_illumina_df.shape}"
        )

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

        abundances_roche_df["sample_id"] = abundances_roche_df\
            ["sample_id"]\
            .apply(lambda row: f"BH10-{str(row)}-Roche")

        LOGGER.info(
            "Number of entries (long) - Roche: "  + \
            f"{abundances_roche_df.shape}"
        )

        # Concatenate
        abundances_df = pd.concat(
            [abundances_illumina_df, abundances_roche_df],
            axis=0,
            ignore_index=True
        )

        LOGGER.info(
            "Number of entries (long) - Illumina and Roche: "  + \
            f"{abundances_df.shape}"
        )

        return abundances_df
