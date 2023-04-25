import os

from modelseedpy import MSBuilder, MSGenome

from cobra.io import load_json_model, save_json_model


# Download protein annotations for a given organims
# Check: https://astrobiomike.github.io/unix/ncbi_eutils
# Check: https://bioinformatics.stackexchange.com/a/16421

ORGANISM = "tez"
DATA_DIR = "../data/genomes/"

# Manually download protein annotations from GenBank's website
# E.g. https://www.ncbi.nlm.nih.gov/nuccore/CP019229.1?report=genbank
genome_path = os.path.join(
    DATA_DIR,
    ORGANISM,
    f"{ORGANISM}.faa"
)
genome = MSGenome.from_fasta(
    genome_path,
    split=" "
)
print('Number of features:', len(genome.features))

model = MSBuilder.build_metabolic_model(
    model_id=ORGANISM,
    genome=genome,
    gapfill_media=None,
    # template=template,
    allow_all_non_grp_reactions=True,
    annotate_with_rast=True
)

# Format for COBRApy
model.compartments = {
    "c": "cytosol",
    "e": "extracellular space"
}

# Save model
model_path = genome_path.replace(".faa", ".json")
save_json_model(
    model=model,
    filename=model_path
)

# Test saved model
load_json_model(model_path)
