import os

from biofoundry.gem import get_gene_counts


def test_get_gene_counts(config: dict) -> None:
    genome_path = os.path.join(
        config["paths"]["genomes"],
        "test_genome.fsa_aa"
    )

    gene_counts = get_gene_counts(genome_path)

    assert gene_counts == 2, "Incorrect gene count!"
