import os
import time

import pandas as pd

from modelseedpy import MSBuilder, MSGenome

from cobra.io import load_json_model, save_json_model


# (DEPRECATED) Download protein annotations for a given organims
# Check: https://astrobiomike.github.io/unix/ncbi_eutils
# Check: https://bioinformatics.stackexchange.com/a/16421

# Manually download protein annotations from GenBank's website
# For example, for tez: https://www.ncbi.nlm.nih.gov/nuccore/CP019229.1

DATA_DIR = "../data/genomes/"

metadata_df = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "genomes-metadata.csv"
    )
)

# Drop rercords with missing data
metadata_df = metadata_df.dropna(axis=0)

for _, row in metadata_df.iterrows():
    organism = row["Code"]

    genome_path = os.path.join(
        DATA_DIR,
        row["Code"],
        row["Protein annotation file"]
    )
    genome = MSGenome.from_fasta(
        genome_path,
        split=" "
    )
    print('Number of features:', len(genome.features))

    model = MSBuilder.build_metabolic_model(
        model_id=organism,
        genome=genome,
        gapfill_media=None,
        # template=template,
        allow_all_non_grp_reactions=True,
        annotate_with_rast=True
    )

    # Save model
    model_path = os.path.join(
        "../data/modelseedpy/",
        f"{organism}.json"
    )
    save_json_model(
        model=model,
        filename=model_path
    )

    # Test saved model
    assert load_json_model(model_path), "[ERROR] Could not read model!"

    # Wait between calls
    time.sleep(0.5)
