import os

import pandas as pd

from micom import Community


DATA_DIR = "../data/genomes/"

# WARNING: Not working for more than one thread, caution is advised!
# N_THREADS = 1


taxonomy = [
    {
        "id": "tez",
        "genus": "Tessaracoccus",
        "species": "Tessaracoccus sp. T2.5-30",
        "reactions": 305,
        "metabolites": 236,
        "file": os.path.join(
            DATA_DIR, "tez/tez.json")
    }
]
taxonomy = pd.DataFrame.from_records(taxonomy * 4)
taxonomy["id"] = taxonomy["id"] + "_" + taxonomy.index.astype(str)

com = Community(taxonomy)
