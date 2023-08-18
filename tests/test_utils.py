import os

import plotly.express as px

from biofoundry.utils import save_fig


def test_save_fig(config: dict) -> None:

    filename = "test-fig.png"

    save_fig(
        fig=px.bar(),
        filename=filename,
        config=config
    )

    fig_path = os.path.join(
        config["paths"]["figures"],
        filename
    )

    assert os.path.exists(fig_path), "Figure not saved correctly!"

    # Clean temporal data
    os.remove(fig_path)
