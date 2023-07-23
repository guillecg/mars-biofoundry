import os

import pandas as pd

import plotly
import plotly.express as px


def get_retropath_results(config: dict) -> pd.DataFrame:
    """
    Get the results of the RetroPath2.0 analysis by inspecting each source
    folder.

    Parameters
    ----------
    config : dict
        The configuration dictionary.

    Returns
    -------
    results_df : pandas.DataFrame
        Dataframe containing the results for each source.

    Examples
    --------
    None

    """

    results_df = pd.DataFrame()

    experiments_dir = os.path.join(
        config["paths"]["retropath"],
        "interesting_metabolites/",
        "experiments/"
    )

    for source in os.listdir(experiments_dir):

        # Get all files for the given source
        source_files = os.listdir(os.path.join(experiments_dir, source))

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
                os.path.join(experiments_dir, source, "source-in-sink.csv")
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

    # Get metabolites not produced
    all_sources = os.listdir(
        os.path.join(
            config["paths"]["retropath"],
            "interesting_metabolites/",
            "sources/"
        )
    )
    all_sources = set([item.replace(".csv", "") for item in all_sources])

    processed_sources = set(results_df["Source"].values)
    not_processed_sources = all_sources - processed_sources

    # TODO: log len(not_processed_sources), not_processed_sources

    not_produced_df = pd.DataFrame({
        "Source": list(not_processed_sources),
        "Status": "Error"
    })

    results_df = pd.concat(
        [results_df, not_produced_df],
        axis=0,
        ignore_index=True
    )

    return results_df


def plot_retropath_results(
    results_df: pd.DataFrame,
    config: dict
) -> plotly.graph_objects.Figure:
    """
    Plot the results of the RetroPath2.0 analysis.

    Parameters
    ----------
    results_df : pandas.DataFrame
        Dataframe containing the results for each source.
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
        text_auto=True,
        template=config["figures"]["template"]
    )
    fig.update_layout(showlegend=False)

    return fig
