import os

import plotly
import plotly.io as pio


def save_fig(
    fig: plotly.graph_objects.Figure,
    filename: str,
    config: dict
) -> None:
    """
    Save the passed figure to the figures directory specified in the config.

    Parameters
    ----------
    fig : plotly.graph_objects.Figure
        The figure to save.
    filename: str
        The name of the output file.
    config : dict
        The configuration dictionary.

    Returns
    -------
    None

    Examples
    --------
    >>> import plotly.express as px
    >>> save_fig(
    >>>     fig=px.bar(),
    >>>     filename="empty_bar.jpg,
    >>>     config={"paths": {"figures": "."}}
    >>> )

    """

    pio.write_image(
        fig=fig,
        file=os.path.join(
            config["paths"]["figures"],
            filename
        ),
        scale=6
    )
