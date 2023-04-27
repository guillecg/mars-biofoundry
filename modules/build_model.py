import os
import time

from modelseedpy import MSBuilder, MSGenome

from cobra.io import load_json_model, save_json_model


# (DEPRECATED) Download protein annotations for a given organims
# Check: https://astrobiomike.github.io/unix/ncbi_eutils
# Check: https://bioinformatics.stackexchange.com/a/16421

# Manually download protein annotations from GenBank's website
# For example, for tez: https://www.ncbi.nlm.nih.gov/nuccore/CP019229.1

DATA_DIR = "../data/genomes/"
GENOMES_DICT = {
    "aci": "QOZT01.1.fsa_aa",
    "bme": "UXHF01P.1.fsa_aa",
    "dmi": "CP003629.1.faa",
    "pse": "CAJFAG01.1.fsa_aa",
    "rhi": "UEYP01.1.fsa_aa",
    "rho": "UWOC01.1.fsa_aa",
    "shw": "CACVBT03.1.fsa_aa",
    "tez": "CP019229.1.faa"
}

for organism, filename in GENOMES_DICT.items():
    genome_path = os.path.join(
        DATA_DIR,
        organism,
        filename
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
