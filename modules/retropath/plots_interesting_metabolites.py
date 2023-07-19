import os

import pandas as pd

from retropath2_wrapper import retropath2


DATA_DIR = "../data/retropath/interesting_metabolites/"


results_df = pd.DataFrame()

for source in os.listdir(os.path.join(DATA_DIR, "experiments")):

    if source.endswith(".csv"):
        continue

    source_files = os.listdir(os.path.join(DATA_DIR, "experiments", source))

    source_df = pd.DataFrame({
        "Source": [source],
        "Status": [None]
    })

    if "_scope.csv" in "".join(source_files):
        source_df["Status"] = "Scope"

    elif sorted(source_files) == sorted(["results.csv", "source-in-sink.csv"]):
        source_df["Status"] = "Results but no scope"

    elif source_files == ["source-in-sink.csv"]:
        source_sink_df = pd.read_csv(
            os.path.join(DATA_DIR, "experiments", source, "source-in-sink.csv")
        )

        if len(source_sink_df):
            source_df["Status"] = "Source in sink"
        else:
            source_df["Status"] = "Source in sink (empty)"

    else:
        source_df["Status"] = "Error"

    results_df = pd.concat(
        [results_df, source_df],
        axis=0,
        ignore_index=True
    )

# Get not produced metabolites
all_sources = os.listdir(os.path.join(DATA_DIR, "sources"))
all_sources = set([item.replace(".csv", "") for item in all_sources])
processed_sources = set(results_df["Source"].values)
not_processed_sources = all_sources - processed_sources
len(not_processed_sources), not_processed_sources

not_produced_df = pd.DataFrame({
    "Source": list(not_processed_sources),
    "Status": "Error"
})

results_df = pd.concat(
    [results_df, not_produced_df],
    axis=0,
    ignore_index=True
)

results_df.to_csv(
    os.path.join(DATA_DIR, "experiments", "interesting_metabolites.csv"),
    index=False,
    header=True
)


import plotly.io as pio
import plotly.express as px

pio.templates.default = "plotly_white"

status_counts = results_df["Status"]\
    .value_counts()\
    .reset_index(drop=False)\
    .rename(columns={"Status": "Counts", "index": "Status"})

fig = px.bar(
    data_frame=status_counts,
    x="Status",
    y="Counts",
    color="Status",
    color_discrete_sequence=px.colors.qualitative.Pastel,
    text_auto=True
)
fig.update_layout(showlegend=False)

pio.write_image(
    fig=fig,
    file="../data/figures/retropath-metabolites.jpg",
    scale=6
)
