import os

import numpy as np
import pandas as pd

import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"

from micom import Community
from micom.media import minimal_medium


DATA_DIR = "../data/micom/rio_tinto/amils_2023/"
DEPTH = 487
TOP_N = 20


taxonomy = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "taxonomy.csv"
    )
)

# Drop species for which there are no abundance data
taxonomy = taxonomy.dropna(axis=0, subset="id")

# Filter by depth
taxonomy["Depth"] = taxonomy["sample_id"]\
    .str.split("-").str[1].astype(int)
taxonomy = taxonomy[taxonomy["Depth"] == DEPTH].copy()

# Build community
com = Community(taxonomy)

sol = com.cooperative_tradeoff()

# Extracellular medium has no growth rate
rates = sol.members.growth_rate.drop("medium")

# Get minimal medium
min_media = minimal_medium(
    community=com,
    community_growth=0.95*sol.growth_rate,
    min_growth=0.95*rates,
    minimize_components=True
)

min_media = min_media\
    .sort_values()\
    .to_frame("flux")\
    .reset_index(names="reaction")

fig = px.bar(
    data_frame=min_media\
        .tail(TOP_N)\
        .rename(columns={"reaction": "Reaction", "flux": "Flux"}),
    x="Reaction",
    y="Flux",
    color="Flux"
)
fig.update_coloraxes(showscale=False)

pio.write_image(
    fig=fig,
    file="../data/figures/medium-limiting-metabolite.jpg",
    scale=6
)


# ---------------------------------------------------------------------------- #


medium_df = pd.read_csv(
    "../data/medium.csv",
    sep=";"
)

# Select depth
medium_df_depth = medium_df[medium_df["Depth"] == DEPTH].copy()

# Remove pH
medium_df_depth = medium_df_depth[medium_df_depth["Species"] != "pH"]

# Get ratio for mapping concentrations to fluxes
nitrate_conc = medium_df_depth.loc[
    (medium_df_depth["Species"] == "Nitrate"), "Concentration (ppm)"
].values[0]

nitrate_flux = min_media.loc[
    (min_media["reaction"] == "EX_cpd00209=no3_m"), "flux"
].values[0]

medium_df_depth["flux"] = \
    medium_df_depth["Concentration (ppm)"] * nitrate_flux / nitrate_conc


# ---------------------------------------------------------------------------- #


import string


modelseed_cpd = pd.read_table(
    "../data/modelseed/compounds.tsv",
    dtype=object # Avoid warnings
)

# Try to get the maximum of chemical species
modelseed_map = modelseed_cpd[
    (modelseed_cpd["name"].isin(medium_df_depth["Species"])) |
    (modelseed_cpd["formula"].isin(medium_df_depth["Species"])) |
    (modelseed_cpd["abbreviation"].isin(medium_df_depth["Species"]))
]
modelseed_map = modelseed_map[
    modelseed_map["source"] == "Primary Database"
].copy()

# Avoid errors in COBRApy due to punctuation and whitespaces in names
modelseed_map["abbreviation"] = modelseed_map["abbreviation"]\
    .str.replace(
        pat=r"[{}|\s]".format(string.punctuation),
        repl="-",
        regex=True
    )

# Create new column to avoid abbreviation duplicates
modelseed_map["abbreviation_id"] = \
    "EX_" + modelseed_map["id"] + "=" + modelseed_map["abbreviation"] + "_m"

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
modelseed_map = modelseed_map[["abbreviation_id", "Species"]].drop_duplicates()

medium_mapped = pd.merge(
    left=medium_df_depth,
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

medium_mapped_min = pd.concat(
    [medium_mapped, min_media],
    axis=0,
    ignore_index=True
).drop_duplicates(keep="first")

medium_mapped_min = dict(zip(
    medium_mapped_min["reaction"],
    medium_mapped_min["flux"]
))


# ---------------------------------------------------------------------------- #


# Add medium
com.medium = medium_mapped_min
com.cooperative_tradeoff()


# ---------------------------------------------------------------------------- #


# from micom.workflows import fix_medium

# candidate_medium = pd.DataFrame.from_records([
#     {"reaction": "EX_glc__D_m", "flux": 10}
# ])
# candidate_medium

# medium = fix_medium(
#     manifest,
#     model_folder=MODEL_DIR,
#     medium=candidate_medium,
#     community_growth=0.1,
#     min_growth=0.01,
#     max_import=10,
#     threads=N_THREADS
# )
# medium