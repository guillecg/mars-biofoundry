import pytest

import yaml


@pytest.fixture(scope="session")
def config(config_path: str = "tests/config-test.yml") -> dict:

    with open(config_path) as config_file:
        config = yaml.safe_load(config_file)

    return config
