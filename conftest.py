import pytest

import yaml

import pandas as pd


@pytest.fixture(scope="session")
def config(config_path: str = "tests/config-test.yml") -> dict:

    with open(config_path) as config_file:
        config = yaml.safe_load(config_file)

    return config


@pytest.fixture(scope="session")
def metadata_df() -> pd.DataFrame:
    return pd.DataFrame.from_records([
        {
            "Species": "NA",
            "Code": "test_model",
            "Database": "NA",
            "ID": "NA",
            "Record": "NA",
            "Download": "NA",
            "Protein annotation file": "NA"
        }
    ])
