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
    # title="Distribution of metabolic model elements across species"
)

pio.write_image(
    fig=fig,
    file="../data/figures/metabolic-models.jpg",
    scale=6
)


# ---------------------------------------------------------------------------- #


import numpy as np
import pandas as pd

import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"


taxonomy = pd.read_csv("../data/micom/rio_tinto/amils_2023/taxonomy.csv")

taxonomy["Code"] = taxonomy["file"]\
    .apply(lambda row: row.split("/")[-1].split("_")[0])

taxonomy["Depth"] = taxonomy["sample_id"]\
    .str.split("-").str[1].astype(int, errors="ignore")

taxonomy["log_abundance"] = np.log1p(taxonomy["abundance"])


taxonomy_matrix = taxonomy.pivot(
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
    width=800,
    height=600,
    # title="Distribution of microbial functions across the vertical column"
)
fig.show()

pio.write_image(
    fig=fig,
    file="../data/figures/amils2023-abundances-depth.jpg",
    scale=6
)


medium_df = pd.read_csv(
    "../data/medium.csv",
    sep=";"
)

fig = px.scatter(
    data_frame=medium_df,
    x="Concentration (ppm)",
    log_x=True,
    y="Depth",
    color="Species",
    color_discrete_sequence=px.colors.qualitative.Pastel,
    #category_orders={"Species": elements_sorted},
    title="Concentration of elements across the vertical column"
)
fig.update_layout(xaxis_title="Concentration log(ppm)")
fig['layout']['yaxis']['autorange'] = "reversed"
fig.show()
