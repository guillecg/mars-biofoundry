import os

import pytest

from biofoundry.gem import (
    get_gene_counts,
    ModelValidator,
    ModelBuilder
)


@pytest.fixture(scope="module")
def model_validator() -> ModelValidator:
    return ModelValidator()


@pytest.fixture(scope="module")
def model_builder() -> ModelValidator:
    return ModelBuilder()


@pytest.fixture(scope="module")
def model_path(config: dict) -> str:
    return os.path.join(
        config["paths"]["models"],
        "test_model.json"
    )


def test_get_gene_counts(config: dict) -> None:
    genome_path = os.path.join(
        config["paths"]["genomes"],
        "test_genome.fsa_aa"
    )

    gene_counts = get_gene_counts(genome_path)

    assert gene_counts == 2, "Incorrect gene count!"


def test_validate_loading(
    model_validator: ModelValidator,
    model_path: str
) -> None:
    assert model_validator.validate_loading(model_path), \
        "Could not load model!"


def test_validate(
    model_validator: ModelValidator,
    model_path: str
) -> None:
    assert model_validator.validate(model_path), \
        "Model does not pass validation checks!"
