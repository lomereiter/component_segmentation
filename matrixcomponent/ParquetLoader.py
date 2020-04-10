import json
import pathlib

import pandas as pd
from sortedcontainers import SortedDict

# pyarrow package must be installed in order to load Parquet;
#
# the code below can be easily adopted to work with CSV:
# it's literally just s/parquet/csv/g

import matrixcomponent.matrix as matrix

def load_parquet(directory):
    directory = pathlib.Path(directory)

    metadata = json.loads((directory / "metadata.json").read_text())
    pangenome_length = metadata['pangenome_length']
    bin_width = metadata['bin_width']

    paths_df = pd.read_parquet(directory / "paths.parquet")
    links_df = pd.read_parquet(directory / "links.parquet")
    ranges_df = pd.read_parquet(directory / "path_bin_ranges.parquet")
    bins_df = pd.read_parquet(directory / "path_bins.parquet")

    paths = []
    for name in paths_df["name"]:
        paths.append(matrix.Path(name))

    for item in bins_df.itertuples():
        bin = matrix.Path.Bin(item.bin_id, item.mean_cov, item.mean_inv, 0, 0)
        paths[item.path_id].bins[item.bin_id] = bin

    for item in ranges_df.itertuples():
        bin = paths[item.path_id].bins[item.bin_id]
        bin.first_nucleotide = item.start
        bin.last_nucleotide = item.end

    for path_id, g in links_df.groupby("path_id"):
        paths[path_id].links = g.to_numpy()[:, 1:3]

    return (paths, pangenome_length, bin_width)
