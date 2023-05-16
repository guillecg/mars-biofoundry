import os

import pandas as pd

import plotly.express as px

from tabula import read_pdf


# NOTE: Most isolates are present at depth 487 (see Table S2 in
# emi16291-sup-0001-supinfo.docx for Amils et al. 2023)


PAPER_DIR = "../data/amils2023/"


# ---------------------------------------------------------------------------- #
# Dataset S2 - ICP-MS elemental analysis of core samples (ppm)

elements_df = read_pdf(
    input_path=os.path.join(
        PAPER_DIR,
        "emi16291-sup-0003-datasets2.pdf"
    ),
    pages=1
)

# Extract the only dataframe since there's only one page
elements_df = elements_df[0]

# Convert to long for plotting
elements_df_long = pd.melt(
    elements_df,
    id_vars=["Depth"],
    value_vars=[
        col for col in elements_df.columns
        if col not in ["Depth"]
    ],
    var_name="Species",
    value_name="Concentration (ppm)"
)

# Convert numeric columns to float
elements_df_long[["Depth", "Concentration (ppm)"]] = \
    elements_df_long[["Depth", "Concentration (ppm)"]]\
    .apply(lambda row: row.str.replace(",", "."))\
    .astype(float)

# Round to fit Illumina and Roche datasets (Datasets S4 and S5)
elements_df_long["Depth"] = elements_df_long["Depth"].astype(int)

# Get sorted elements by their maximum concentration
elements_sorted = elements_df_long\
    .groupby("Species")\
    .max()\
    .sort_values("Concentration (ppm)", ascending=False)\
    .index\
    .to_list()

fig = px.scatter(
    data_frame=elements_df_long,
    x="Concentration (ppm)",
    y="Depth",
    color="Species",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Species": elements_sorted},
    title="Concentration of elements across the vertical column"
)
fig['layout']['yaxis']['autorange'] = "reversed"
fig.show()

fig = px.scatter(
    data_frame=elements_df_long,
    x="Concentration (ppm)",
    log_x=True,
    y="Depth",
    color="Species",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Species": elements_sorted},
    title="Concentration of elements across the vertical column"
)
fig.update_layout(xaxis_title="Concentration log(ppm)")
fig['layout']['yaxis']['autorange'] = "reversed"
fig.show()

fig = px.violin(
    data_frame=elements_df_long,
    y="Concentration (ppm)",
    log_y=True,
    color="Species",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Species": elements_sorted},
    title="Distribution of concentrations per element"
)
fig.update_layout(
    xaxis_title="Species",
    yaxis_title="Concentration log(ppm)"
)
fig.show()


# ---------------------------------------------------------------------------- #
# Dataset S3 - Ionic chromatography of BH10 soluble organic and iniorganic anions (ppm)

compounds_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
        "emi16291-sup-0004-datasets3.xls"
    ),
    sheet_name="BH10_CI",
    skiprows=1
)

# Create depth column
compounds_df["Depth"] = compounds_df["SAMPLE"]\
    .str.split("_").str[-1]\
    .str.split("-").str[-1]

compounds_df = compounds_df.drop("SAMPLE", axis=1)

# Remove letters from numbers (e.g. depth W90)
compounds_df["Depth"] = compounds_df["Depth"]\
    .str.strip(r"\W")

# Convert to numeric
compounds_df["Depth"] = compounds_df["Depth"]\
    .str.replace(",", ".")\
    .astype(float)

# Round to fit Illumina and Roche datasets (Datasets S4 and S5)
compounds_df["Depth"] = compounds_df["Depth"].astype(int)

compounds_df_long = pd.melt(
    compounds_df,
    id_vars=["Depth"],
    value_vars=[
        col for col in compounds_df.columns
        if col not in ["Depth"]
    ],
    var_name="Species",
    value_name="Concentration (ppm)"
)

# Capitalize to fit format
compounds_df_long["Species"] = compounds_df_long["Species"]\
    .str.capitalize()

# Replace capitalized pH
compounds_df_long["Species"] = compounds_df_long["Species"]\
    .str.replace("Ph", "pH")

# ---------------------------------------------------------------------------- #
# Table S1 - Soluble cations (ppm)

cations_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
        "emi16291-sup-0001-supinfo-tables1.ods"
    ),
    sheet_name="Sheet1",
    skiprows=1
)

# Round to fit Illumina and Roche datasets (Datasets S4 and S5)
cations_df["Depth"] = cations_df["Depth"].astype(int)

cations_df_long = pd.melt(
    cations_df,
    id_vars=["Depth"],
    value_vars=[
        col for col in cations_df.columns
        if col not in ["Depth"]
    ],
    var_name="Species",
    value_name="Concentration (ppm)"
)


# ---------------------------------------------------------------------------- #
# Table S7 - Occluded gases and natural activities at different depths (10/8/22)

gases_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
        "emi16291-sup-0001-supinfo-tables7.ods"
    ),
    sheet_name="Sheet1",
    skiprows=1
)

# Drop last row containing the explanation
gases_df = gases_df.iloc[:-1, :].copy()

# Drop last three columns since they correspond to activate metabolism
gases_df = gases_df.iloc[:, :-3].copy()

# Rename depth column
gases_df = gases_df.rename(columns={"Depth mbs": "Depth"})

# Round to fit Illumina and Roche datasets (Datasets S4 and S5)
gases_df["Depth"] = gases_df["Depth"].astype(int)

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

gases_df["H2"] = gases_df["H2"].replace(h2_co2_symbol_map).astype(float)
gases_df["CO2"] = gases_df["CO2"].replace(h2_co2_symbol_map).astype(float)
gases_df["CH4"] = gases_df["CH4"].replace(ch4_symbol_map).astype(float)

gases_df_long = pd.melt(
    gases_df,
    id_vars=["Depth"],
    value_vars=[
        col for col in gases_df.columns
        if col not in ["Depth"]
    ],
    var_name="Species",
    value_name="Concentration (ppm)"
)

# ---------------------------------------------------------------------------- #
# Table S8-2 - Number of microbial species detected at different depths which
# have the potential of carry out key metabolic pathways of the C, H, N, S and
# Fe cycles. 

microbes_df = pd.read_excel(
    os.path.join(
        PAPER_DIR,
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
    width=600,
    height=600,
    title="Distribution of microbial functions across the vertical column"
)
fig.show()

complete_depths = microbes_df[microbes_df["Pathway"] == "nº compl cyc"]
complete_depths = complete_depths[complete_depths == 5]\
    .dropna(axis=1)\
    .columns

print(
    "[INFO] Depths with all microbial functions analyzed:",
    complete_depths
)

fig = px.imshow(
    img=microbes_df[complete_depths].T.to_numpy(),
    x=microbes_df["Pathway"],
    y=complete_depths,
    labels=dict(
        x="Pathway",
        y="Depth",
        color="Count"
    ),
    aspect="equal",
    width=600,
    height=600,
    title="Distribution of microbial functions (only complete functions)"
)
fig.add_shape(
    type="rect",
    x0=0.0,
    x1=1.0,
    y0=481,
    y1=493,
    xref="paper",
    yref="y",
    line_color="red"
)
fig.add_shape(
    type="rect",
    x0=0.0,
    x1=1.0,
    y0=587,
    y1=625,
    xref="paper",
    yref="y",
    line_color="red"
)
fig.show()


# ---------------------------------------------------------------------------- #
# Merge all data
medium_df = pd.concat(
    [
        elements_df_long,
        compounds_df_long,
        cations_df_long,
        gases_df_long
    ],
    axis=0,
    ignore_index=True
)


# ---------------------------------------------------------------------------- #
# Taxonomy abundances

taxonomy_abundances = pd.read_csv(
    os.path.join(
        "../data/micom/rio_tinto/amils_2023/",
        "taxonomy.csv"
    ),
    sep=","
)

# Create depth column
taxonomy_abundances["Depth"] = taxonomy_abundances["sample_id"]\
    .str.split("-").str[1]
