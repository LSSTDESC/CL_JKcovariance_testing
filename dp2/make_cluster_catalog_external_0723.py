#!/opt/lsst/software/stack/conda/envs/lsst-scipipe-13.0.0/bin/python
"""
Convert the raw external x DP2 cluster cross-match CSV into the flat
"clusters" HDF5 group TXPipe's CLClusterBinningRedshiftRichness expects,
so it can do its own z/richness binning.
"""

import h5py
import pandas as pd

INFILE = "/pscratch/sd/k/kabelo/dp2/external_0723_x_dp2.csv"
OUTFILE = "/pscratch/sd/k/kabelo/dp2/external_0723_x_dp2.hdf5"

COLUMN_MAP = {
    "catalog_id_clusters": "cluster_id",
    "ra": "ra",
    "dec": "dec",
    "redshift": "redshift",
    "redshift_unc_clusters": "redshift_err",
    "richness": "richness",
    "richness_unc_clusters": "richness_err",
}

df = pd.read_csv(INFILE)
print(f"{len(df)} clusters read from {INFILE}")

with h5py.File(OUTFILE, "w") as f:
    g = f.create_group("clusters")
    for src_col, out_col in COLUMN_MAP.items():
        g.create_dataset(out_col, data=df[src_col].to_numpy())
    # not used anywhere downstream (grepped the full txpipe package: only
    # ever read-and-copied by bin_cluster.py, never used in a calculation)
    g.create_dataset("scaleval", data=[1.0] * len(df))

print(f"Wrote {OUTFILE}")
