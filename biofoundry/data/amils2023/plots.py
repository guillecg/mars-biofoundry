import os

import pandas as pd

import plotly
import plotly.express as px


def plot_concentrations(
    data_df_long: pd.DataFrame,
    data_type: str,
    config: dict
) -> None:
    """
    Plot the different chemical species in Amils et al. 2023.

    Parameters
    ----------
    data_df_long : pandas.DataFrame
        Dataframe containing the chemical species in long format.
    config : dict
        The configuration dictionary.

    Returns
    -------
    None

    Examples
    --------
    None

    """

    # Get sorted elements by their maximum concentration
    data_sorted = data_df_long\
        .groupby("Species")\
        .max()\
        .sort_values("Concentration (ppm)", ascending=False)\
        .index\
        .to_list()

    fig = px.scatter(
        data_frame=data_df_long,
        x="Concentration (ppm)",
        y="Depth",
        color="Species",
        color_discrete_sequence=px.colors.qualitative.Pastel,
        category_orders={"Species": data_sorted},
        template=config["figures"]["template"],
        title=f"Concentration of {data_type}s across the vertical column"
    )
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig.show("png")

    fig = px.scatter(
        data_frame=data_df_long,
        x="Concentration (ppm)",
        log_x=True,
        y="Depth",
        color="Species",
        color_discrete_sequence=px.colors.qualitative.Pastel,
        category_orders={"Species": data_sorted},
        template=config["figures"]["template"],
        title=f"Concentration of {data_type}s across the vertical column"
    )
    fig.update_layout(xaxis_title="Concentration log(ppm)")
    fig['layout']['yaxis']['autorange'] = "reversed"
    fig.show("png")

    fig = px.violin(
        data_frame=data_df_long,
        y="Concentration (ppm)",
        log_y=True,
        color="Species",
        color_discrete_sequence=px.colors.qualitative.Pastel,
        category_orders={"Species": data_sorted},
        template=config["figures"]["template"],
        title=f"Distribution of concentrations per {data_type}"
    )
    fig.update_layout(
        xaxis_title="Species",
        yaxis_title="Concentration log(ppm)"
    )
    fig.show("png")


def plot_microbial_data(
    microbes_df: pd.DataFrame,
    config: dict
) -> plotly.graph_objects.Figure:
    """
    Plot the cycles present in the community by depth.

    Parameters
    ----------
    microbes_df : pandas.DataFrame
        Dataframe containing the microbial data from Amils et al. 2023.
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

    numeric_cols = [
        col for col in microbes_df.columns
        if col not in ["Pathway"]
    ]

    fig = px.imshow(
        img=microbes_df[numeric_cols].T.to_numpy(),
        x=microbes_df["Pathway"],
        y=numeric_cols,
        labels=dict(
            x="Pathway",
            y="Depth",
            color="Count"
        ),
        aspect="equal",
        template=config["figures"]["template"]
    )

    return fig
