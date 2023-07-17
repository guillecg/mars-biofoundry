import os

import pandas as pd

from rdkit import Chem

from retropath2_wrapper import retropath2


DATA_DIR = "../data/retropath/interesting_metabolites/"

metabolites_df = pd.read_csv(
    os.path.join(
        DATA_DIR,
        "metabolitos_producibles.csv"
    )
)

# Get InChIs
metabolites_df["InChI"] = \
    metabolites_df["Smile"]\
    .dropna(how="all")\
    .apply(lambda row: Chem.MolToInchi(Chem.MolFromSmiles(row)))

metabolites_df.to_csv(
    os.path.join(
        DATA_DIR,
        "metabolitos_producibles_inchi.csv"
    ),
    header=True,
    index=False
)

# Drop metabolites withou InChI
metabolites_df = metabolites_df.dropna(subset="InChI")

# Rename to fit RetroPath2.0 format
metabolites_df = metabolites_df\
    .rename(columns={"Molecule": "Name"})

metabolites_df["Name"] = metabolites_df["Name"]\
    .str.lower()\
    .str.replace(" ", "_")

# Create a source file for each compound
for _, row in metabolites_df.iterrows():
    compound_name = row["Name"]
    row = row.to_frame().T[["Name", "InChI"]]

    row.to_csv(
        os.path.join(
            DATA_DIR,
            "sources",
            f"{compound_name}.csv"
        ),
        header=True,
        index=False
    )

for source_file in os.listdir(os.path.join(DATA_DIR, "sources")):

    r_code = retropath2(
        rules_file=os.path.abspath(
            "../data/retropath/ec_numbers_all_rules.csv"
        ),
        sink_file=os.path.abspath(
            "../data/retropath/ec_numbers_all_sink.csv"
        ),
        source_file=os.path.abspath(
            os.path.join(DATA_DIR, "sources", source_file)
        ),
        outdir=os.path.abspath(
            os.path.join(
                DATA_DIR,
                "experiments",
                os.path.splitext(source_file)[0]
            )
        ),
        dmax=16,
        dmin=6,
        max_steps=10,
        topx=100,
        mwmax_source=1000
    )
