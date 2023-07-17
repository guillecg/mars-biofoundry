import os
import glob

import json

import pandas as pd

import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"


DATA_DIR = "../data/modelseedpy/"


plot_df = pd.DataFrame()

for filepath in glob.glob(os.path.join(DATA_DIR, "*_formatted.json")):
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
    width=750,
    height=500,
    title="Distribution of metabolic model elements across species"
)

pio.write_image(
    fig=fig,
    file="../data/figures/metabolic-models.jpg",
    scale=6
)
