
import pandas as pd

import plotly.express as px

from tabula import read_pdf


elements_df = read_pdf(
    input_path="../data/amils2023/emi16291-sup-0003-datasets2.pdf",
    pages=1
)

# Extract the only dataframe since there's only one page
elements_df = elements_df[0]

# Convert to long for plotting
elements_df_long = pd.melt(
    elements_df,
    id_vars=["Depth"],
    value_vars=[col for col in elements_df.columns if col != "Depth"],
    var_name="Element",
    value_name="Concentration (ppm)"
)

# Convert numeric columns to float
elements_df_long[["Depth", "Concentration (ppm)"]] = \
    elements_df_long[["Depth", "Concentration (ppm)"]]\
    .apply(lambda row: row.str.replace(",", "."))\
    .astype(float)

# Get sorted elements by their maximum concentration
elements_sorted = elements_df_long\
    .groupby("Element")\
    .max()\
    .sort_values("Concentration (ppm)", ascending=False)\
    .index\
    .to_list()

fig = px.scatter(
    data_frame=elements_df_long,
    x="Concentration (ppm)",
    y="Depth",
    color="Element",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Element": elements_sorted},
    title="Concentration of elements across the vertical column"
)
fig['layout']['yaxis']['autorange'] = "reversed"
fig.show()

fig = px.scatter(
    data_frame=elements_df_long,
    x="Concentration (ppm)",
    log_x=True,
    y="Depth",
    color="Element",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Element": elements_sorted},
    title="Concentration of elements across the vertical column"
)
fig.update_layout(xaxis_title="Concentration log(ppm)")
fig['layout']['yaxis']['autorange'] = "reversed"
fig.show()

fig = px.violin(
    data_frame=elements_df_long,
    y="Concentration (ppm)",
    log_y=True,
    color="Element",
    color_discrete_sequence=px.colors.qualitative.Alphabet,
    category_orders={"Element": elements_sorted},
    title="Distribution of concentrations per element"
)
fig.update_layout(
    xaxis_title="Element",
    yaxis_title="Concentration log(ppm)"
)
fig.show()
