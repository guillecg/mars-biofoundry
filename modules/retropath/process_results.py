import os

import pandas as pd


EXPERIMENT_NAME = "2023-04-19"
ORGANISM = "tez"


results_filepath = os.path.join(
    "../data/retropath/experiments/",
    f"{EXPERIMENT_NAME}/kegg_{ORGANISM}_result.csv"
)

results_df = pd.read_csv(results_filepath)

# Filter by first iterations
results_df = results_df[results_df["Iteration"].isin((0, 1))]

# Filter by "in sink"
results_df = results_df[results_df["In Sink"] == 1]

# Write to file
final_filepath = results_filepath.replace(".csv", "_filtered.csv")
results_df.to_csv(
    final_filepath,
    header=True,
    index=False,
    sep=","
)

rp2paths_dir = os.path.join(
    os.path.dirname(final_filepath),
    "rp2paths_results"
)

print(
    "[+] Run the following command to enumerate the paths:\n" + \
    f"python -m rp2paths all {final_filepath} " + \
    f"--outdir {rp2paths_dir}"
)
