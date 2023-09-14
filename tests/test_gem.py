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
def model_builder(
    config: dict,
    model_validator: ModelValidator
) -> ModelBuilder:
    return ModelBuilder(
        config=config,
        model_validator=model_validator
    )


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


def test_format_model(
    model_builder: ModelBuilder,
    model_path: str
) -> None:
    model_path_formatted = model_builder.format_model(model_path)
    model_path_expected = os.path.join(
        os.path.dirname(model_path),
        "expected",
        os.path.basename(model_path)
    )
    model_path_expected = model_path_expected.replace(
        ".json",
        "_formatted.json"
    )

    with open(model_path_formatted, "r") as fh:
        model_formatted = fh.read()

    with open(model_path_expected, "r") as fh:
        model_expected = fh.read()

    # Clean temporal data
    os.remove(model_path_formatted)

    assert model_expected == model_formatted, \
        "Model is not correctly formatted!"
