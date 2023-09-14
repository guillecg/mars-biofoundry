import os

import pytest

import json

import pandas as pd
from pandas.testing import assert_frame_equal

from biofoundry.retropath.preloader import RetroPathPreloader


@pytest.fixture(scope="module")
def preloader(config: dict) -> RetroPathPreloader:
    return RetroPathPreloader(config)


@pytest.fixture(scope="module")
def model_path(config: dict) -> str:
    return os.path.join(
        config["paths"]["models"],
        "test_model.json"
    )


def test_get_ec_from_model(
    config: dict,
    model_path: str,
    preloader: RetroPathPreloader
) -> None:

    modelseed_path = os.path.join(
        config["paths"]["modelseed"],
        "reactions.tsv"
    )
    modelseed_reactions = pd.read_table(modelseed_path)

    # Load model
    with open(model_path, mode="r") as fh:
        model_dict = json.loads(fh.read())

    ec_numbers = preloader.get_ec_from_model(
        model_dict=model_dict,
        modelseed_reactions=modelseed_reactions
    )

    assert len(ec_numbers) == 3, "EC numbers were not correctly extracted!"


def test_get_ec_numbers(
    config: dict,
    metadata_df: pd.DataFrame,
    preloader: RetroPathPreloader
) -> None:

    preloader.get_ec_numbers(metadata=metadata_df)

    ec_path = os.path.join(
        config["paths"]["retropath"],
        config["retropath"]["files"]["ec_numbers"]
    )
    ec_numbers_df = pd.read_csv(ec_path)

    # Clean temporal data
    os.remove(ec_path)

    ec_path_expected = os.path.join(
        config["paths"]["retropath"],
        "expected",
        config["retropath"]["files"]["ec_numbers"]
    )
    ec_numbers_df_expected = pd.read_csv(ec_path_expected)

    assert_frame_equal(
        left=ec_numbers_df,
        right=ec_numbers_df_expected
    )
