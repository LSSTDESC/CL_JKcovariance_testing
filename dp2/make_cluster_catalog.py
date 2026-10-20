#!/usr/bin/env python
"""
Build the DP2 cluster catalog for TXPipe from the DESC data registry.

Reads registry dataset 133 (dp2_known_clusters 0.1.0), keeps the DES Y3
redMaPPer clusters inside a DP2 tract, and writes the flat "clusters" HDF5
group that CLClusterBinningRedshiftRichness and CLClusterShearCatalogs expect.

redMaPPer only: the registry mixes richness definitions (redMaPPer lambda,
HSC N_mem) and its catalog_id is only unique within one catalogue.
The previous catalog (1" DP2 object match of all catalogues) is built by
make_cluster_catalog_external_0723.py.

Run with the TXPipe env (needs dataregistry):
    module load python; conda activate /global/cfs/cdirs/m1727/kabelo/TXPipe/conda
"""

from datetime import date

import h5py
import numpy as np
import pandas as pd
from dataregistry import DataRegistry

DATASET_ID = 133
CATALOG_SOURCE = "DES Y3 RedMaPPer"
OUTFILE = "/pscratch/sd/k/kabelo/dp2/dp2_known_clusters_0.1.0_redmapper.hdf5"

COLUMN_MAP = {
    "catalog_id": "cluster_id",
    "RA": "ra",
    "DEC": "dec",
    "redshift": "redshift",
    "redshift_unc": "redshift_err",
    "richness": "richness",
    "richness_unc": "richness_err",
    "dp2_tract": "dp2_tract",
    "bcg_specz": "bcg_specz",  # -1 where there is no BCG spec-z
}

query = DataRegistry().Query
infile = query.get_dataset_absolute_path(DATASET_ID)
meta = query.find_datasets(
    ["dataset.name", "dataset.version_string"],
    [query.gen_filter("dataset.dataset_id", "==", DATASET_ID)],
    return_format="dataframe",
).iloc[0]

df = pd.read_csv(infile)
print(f"{len(df)} rows read from {infile}")
df = df[(df["catalog_source"] == CATALOG_SOURCE) & (df["dp2_tract"] >= 0)].copy()
print(f"{len(df)} {CATALOG_SOURCE} clusters in a DP2 tract")
# catalog_id is read as strings (other catalogues use names); DES ids are integers
df["catalog_id"] = df["catalog_id"].astype(np.int64)

# the registry uses -1 for missing values; none are expected for redMaPPer
for col in ["redshift", "richness", "redshift_unc", "richness_unc"]:
    n_bad = int((df[col] < 0).sum())
    assert n_bad == 0, f"{n_bad} rows with negative {col}"
# TXPipe looks clusters up by cluster_id, so it must be unique
assert df["catalog_id"].is_unique, "catalog_id is not unique"

with h5py.File(OUTFILE, "w") as f:
    g = f.create_group("clusters")
    for src_col, out_col in COLUMN_MAP.items():
        g.create_dataset(out_col, data=df[src_col].to_numpy())
    # not used anywhere downstream, only copied by bin_cluster.py
    g.create_dataset("scaleval", data=np.ones(len(df)))
    g.attrs["registry_dataset_id"] = DATASET_ID
    g.attrs["registry_dataset_name"] = meta["dataset.name"]
    g.attrs["registry_dataset_version"] = meta["dataset.version_string"]
    g.attrs["registry_path"] = infile
    g.attrs["selection"] = f"catalog_source == '{CATALOG_SOURCE}' & dp2_tract >= 0"
    g.attrs["build_date"] = date.today().isoformat()

print(f"Wrote {len(df)} clusters to {OUTFILE}")
