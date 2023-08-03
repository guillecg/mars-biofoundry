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
        .value_counts(normalize=True)\
        .apply(lambda row: 100 * row)\
        .reset_index(drop=False)\
        .rename(columns={"Status": "% Counts", "index": "Status"})

    fig = px.bar(
        data_frame=status_counts,
        x="Status",
        y="% Counts",
        color="Status",
        color_discrete_sequence=px.colors.qualitative.Pastel,
        text=[f"{i} %" for i in status_counts["% Counts"].round(2)],
        template=config["figures"]["template"]
    )
    fig.update_layout(showlegend=False)

    return fig


def get_classes_counts(
    results_df: pd.DataFrame,
    config: dict
) -> pd.DataFrame:
    """
    Get the compound classes to further analyse RetroPath results.

    Parameters
    ----------
    results_df : pandas.DataFrame
        Dataframe containing the results for each source obtained with function
        get_retropath_results.
    config : dict
        The configuration dictionary.

    Returns
    -------
    _ : pandas.DataFrame
        Dataframe containing the counts per status and class.

    Examples
    --------
    None

    """

    # Load compound SMILES and names
    detectable_df = pd.read_excel(
        os.path.join(
            config["paths"]["retropath_classes"],
            "Detectables_list.xlsx"
        ),
        header=None,
        names=["ID", "SMILES"] + list(range(11)) + ["Source"]
    )
    producible_df = pd.read_excel(
        os.path.join(
            config["paths"]["retropath_classes"],
            "Producibles_list.xlsx"
        ),
        header=None,
        names=["ID", "SMILES"] + list(range(7)) + ["Source"]
    )

    # Drop separator columns
    detectable_df = detectable_df.dropna(how="all", axis=1)
    producible_df = producible_df.dropna(how="all", axis=1)

    smiles_df = pd.concat(
        [detectable_df, producible_df],
        axis=0,
        ignore_index=True
    )

    # Load compound classes
    class_df = pd.concat(
        [
            pd.read_table(
                os.path.join(
                    config["paths"]["retropath_classes"],
                    "AnotacionDetectables.txt"
                ),
                names=["Index", "ID", "Class"]
            ),
            pd.read_table(
                os.path.join(
                    config["paths"]["retropath_classes"],
                    "AnotacionDetectablesEnvipath.txt"
                ),
                names=["Index", "ID", "Class"]
            ),
            pd.read_table(
                os.path.join(
                    config["paths"]["retropath_classes"],
                    "AnotacionProducibles.txt"
                ),
                names=["Index", "ID", "Class"]
            ),
            pd.read_table(
                os.path.join(
                    config["paths"]["retropath_classes"],
                    "AnotacionProduciblesEnvipath.txt"
                ),
                names=["Index", "ID", "Class"]
            )
        ],
        axis=0,
        ignore_index=True
    )

    # Drop potential duplicates
    smiles_df = smiles_df.drop_duplicates()
    class_df = class_df.drop_duplicates()

    merged_df = pd.merge(
        left=smiles_df,
        right=class_df,
        on="ID",
        how="left"
    )
    merged_df = merged_df[["Source", "Class"]]

    # Lower case to match the other dataframe
    merged_df["Source"] = merged_df["Source"].str.lower()
    results_df["Source"] = results_df["Source"].str.lower()

    results_class_df = pd.merge(
        left=results_df,
        right=merged_df,
        on="Source",
        how="left"
    )

    return results_class_df\
        .groupby(["Status", "Class"], as_index=False)\
        .count()\
        .rename(columns={"Source": "Counts"})\
        .sort_values("Counts", ascending=False)


def plot_classes_counts(
    counts_df: pd.DataFrame,
    config: dict,
    status: str = "Scope",
    counts_thr: int = 2
) -> plotly.graph_objects.Figure:
    """
    Plot the counts per class.

    Parameters
    ----------
    counts_df : pandas.DataFrame
        Dataframe containing the counts per status and class obtained with
        function get_classes_counts.
    config : dict
        The configuration dictionary.
    status : str
        The status to be plotted. Must be in counts_df["Status"]
    counts_thr: int
        Thresholt to avoid plotting too many points.

    Returns
    -------
    fig : plotly.graph_objects.Figure
        The generated figure.

    Examples
    --------
    None

    """

    fig = px.bar(
        data_frame=counts_df[
            (counts_df["Counts"] >= counts_thr) &
            (counts_df["Status"] == status)
        ],
        x="Class",
        y="Counts",
        color="Class",
        category_orders={
            "Status": [
                "Results but no scope",
                "Source in sink",
                "Source in sink (empty)",
                "Error",
                "Scope"
            ]
        },
        color_discrete_sequence=px.colors.qualitative.Pastel,
        template=config["figures"]["template"],
        height=600,
        width=900
    )
    fig.update_layout(
        showlegend=False,
        xaxis=dict(
            tickangle=45,
            tickfont=dict(size=11)
        )
    )

    return fig
