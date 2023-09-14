import os

import copy

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

    _ = preloader.get_ec_numbers(metadata=metadata_df)

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


def test_get_rules(
    config: dict,
    preloader: RetroPathPreloader
) -> None:

    # Change input files by the expected ones
    config_modified = copy.deepcopy(config)
    config_modified["retropath"]["files"]["ec_numbers"] = os.path.join(
        os.path.dirname(config_modified["retropath"]["files"]["ec_numbers"]),
        "expected",
        os.path.basename(config_modified["retropath"]["files"]["ec_numbers"])
    )

    preloader_modified = copy.deepcopy(preloader)
    preloader_modified.config = config_modified

    _ = preloader_modified.get_rules()

    rules_path = os.path.join(
        config["paths"]["retropath"],
        config["retropath"]["files"]["rules"]
    )
    rules_df = pd.read_csv(rules_path)

    # Clean temporal data
    os.remove(rules_path)

    rules_path_expected = os.path.join(
        config["paths"]["retropath"],
        "expected",
        config["retropath"]["files"]["rules"]
    )
    rules_df_expected = pd.read_csv(rules_path_expected)

    assert_frame_equal(
        left=rules_df,
        right=rules_df_expected
    )


def test_get_sink(
    config: dict,
    preloader: RetroPathPreloader
) -> None:

    # Change input files by the expected ones
    config_modified = copy.deepcopy(config)
    config_modified["retropath"]["files"]["ec_numbers"] = os.path.join(
        os.path.dirname(config_modified["retropath"]["files"]["ec_numbers"]),
        "expected",
        os.path.basename(config_modified["retropath"]["files"]["ec_numbers"])
    )
    config_modified["retropath"]["files"]["rules"] = os.path.join(
        os.path.dirname(config_modified["retropath"]["files"]["rules"]),
        "expected",
        os.path.basename(config_modified["retropath"]["files"]["rules"])
    )

    preloader_modified = copy.deepcopy(preloader)
    preloader_modified.config = config_modified

    _ = preloader_modified.get_sink()

    sink_path = os.path.join(
        config["paths"]["retropath"],
        config["retropath"]["files"]["sink"]
    )
    sink_df = pd.read_csv(sink_path)

    # Clean temporal data
    os.remove(sink_path)

    sink_path_expected = os.path.join(
        config["paths"]["retropath"],
        "expected",
        config["retropath"]["files"]["sink"]
    )
    sink_df_expected = pd.read_csv(sink_path_expected)

    assert_frame_equal(
        left=sink_df,
        right=sink_df_expected
    )


def test_get_sources(
    config: dict,
    preloader: RetroPathPreloader
) -> None:

    sources_df = preloader.get_sources()

    sources_dir = os.path.join(
        config["paths"]["retropath"],
        "interesting_metabolites/sources"
    )
    sources_list = os.listdir(sources_dir)

    # Clean temporal data
    for filename in sources_list:
        os.remove(
            os.path.join(
                sources_dir,
                filename
            )
        )

    assert len(sources_list) == len(sources_df), \
        "Sources files were not correctly created!"
