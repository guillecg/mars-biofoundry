import os

import pytest

import pandas as pd

import plotly

from biofoundry.gem import (
    get_gene_counts,
    ModelValidator,
    ModelBuilder,
    plot_metabolic_models
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
        model = fh.read()

    with open(model_path_expected, "r") as fh:
        model_expected = fh.read()

    # Clean temporal data
    os.remove(model_path_formatted)

    assert model == model_expected, \
        "Model is not correctly formatted!"


def test_rename_compartments(
    model_builder: ModelBuilder,
    model_path: str
) -> None:

    model_path_expected = os.path.join(
        os.path.dirname(model_path),
        "expected",
        os.path.basename(model_path)
    )
    model_path_expected = model_path_expected.replace(
        ".json",
        "_renamed_compartments.json"
    )

    with open(model_path, "r") as fh:
        model_text = fh.read()
        model_text = model_builder.rename_compartments(model_text)

    with open(model_path_expected, "r") as fh:
        model_text_expected = fh.read()

    assert model_text == model_text_expected, \
        "Compartments are not correctly renamed!"


def test_rename_metabolites(
    model_builder: ModelBuilder,
    model_path: str
) -> None:

    model_path_expected = os.path.join(
        os.path.dirname(model_path),
        "expected",
        os.path.basename(model_path)
    )
    model_path_expected = model_path_expected.replace(
        ".json",
        "_renamed_metabolites.json"
    )

    with open(model_path, "r") as fh:
        model_text = fh.read()
        model_text = model_builder.rename_metabolites(model_text)

    with open(model_path_expected, "r") as fh:
        model_text_expected = fh.read()

    assert model_text == model_text_expected, \
        "Metabolites are not correctly renamed!"


def test_plot_metabolic_models(config: dict) -> None:

    # Use expected model JSON file
    config["paths"]["models"] = os.path.join(
        config["paths"]["models"],
        "expected"
    )

    metadata_df = pd.DataFrame({
        "Organism": ["organism"],
        "Genes (annotation)": [2]
    })

    fig = plot_metabolic_models(
        metadata_df=metadata_df,
        config=config
    )

    assert type(fig) == plotly.graph_objects.Figure, \
        "Plot was not correctly generated"
