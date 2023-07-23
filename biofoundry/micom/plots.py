import numpy as np
import pandas as pd

import plotly
import plotly.express as px


def plot_abundances_depth(
    taxonomy_df: pd.DataFrame,
    config: dict
) -> plotly.graph_objects.Figure:
    """
    Plot the abundances of each species at different depths.

    Parameters
    ----------
    taxonomy_df : pandas.DataFrame
        Dataframe generated by MICOMPreloader containing both the taxonomy and
        the abundances of each species at different depths.
    config : dict
        The configuration dictionary.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The generated figure.

    Examples
    --------
    None

    """

    taxonomy_df["Code"] = taxonomy_df["file"]\
        .apply(lambda row: row.split("/")[-1].split("_")[0])

    taxonomy_df["Depth"] = taxonomy_df["sample_id"]\
        .str.split("-").str[1].astype(int, errors="ignore")

    taxonomy_df["log_abundance"] = np.log1p(taxonomy_df["abundance"])

    # Convert to long
    taxonomy_matrix = taxonomy_df.pivot(
        columns="Code",
        index="Depth",
        values="log_abundance"
    )

    fig = px.imshow(
        img=taxonomy_matrix,
        x=taxonomy_matrix.columns,
        y=taxonomy_matrix.index,
        labels=dict(
            x="Species",
            y="Depth",
            color="log(Abundance)"
        ),
        aspect="auto",
        template=config["figures"]["template"],
        title="Distribution of microbial functions across the vertical column"
    )

    return fig
