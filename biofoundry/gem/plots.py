import os
import glob

import json

import pandas as pd

import plotly
import plotly.express as px


def plot_metabolic_models(config: dict) -> plotly.graph_objects.Figure:
    """
    Plot the number of genes, reactions and metabolites per model.

    Parameters
    ----------
    config : dict
        The configuration dictionary.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The generated figure.

    Examples
    --------
    >>> plot_metabolic_models(config)

    """

    plot_df = pd.DataFrame()

    glob_pattern = os.path.join(
        config["paths"]["models"],
        "*_formatted.json"
    )

    for filepath in glob.glob(glob_pattern):
        with open(filepath, "r") as fh:
            model_dict = json.load(fh)

            model_df = pd.Series({
                "Organism": model_dict["id"],
                "Genes": len(model_dict["genes"]),
                "Reactions": len(model_dict["reactions"]),
                "Metabolites": len(model_dict["metabolites"])
            })
            model_df = model_df.to_frame().T

            plot_df = pd.concat(
                [plot_df, model_df],
                axis=0,
                ignore_index=True
            )

    # Convert to long format for Plotly
    plot_df_long = pd.melt(
        frame=plot_df,
        id_vars=["Organism"],
        value_vars=[
            "Genes",
            "Reactions",
            "Metabolites"
        ],
        var_name="Element",
        value_name="Counts",
        ignore_index=True,
    )

    # Sort by ID
    plot_df_long = plot_df_long.sort_values("Organism")

    fig = px.bar(
        data_frame=plot_df_long,
        x="Organism",
        y="Counts",
        color="Element",
        barmode="group",
        category_orders={
            "Element": ["Genes", "Reactions", "Metabolites"]
        },
        color_discrete_sequence=px.colors.qualitative.Pastel,
        template=config["figures"]["template"],
        title="Distribution of metabolic model elements across species"
    )

    return fig
