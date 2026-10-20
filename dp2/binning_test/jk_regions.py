"""Count occupied jackknife regions (unique healpix pixels of clusters that have
background sources) per (z, richness) bin, for candidate binnings and nsides."""
import sys
import h5py
import healpy
import numpy as np

OUT = "/pscratch/sd/k/kabelo/dp2/outputs"
with h5py.File("/pscratch/sd/k/kabelo/dp2/external_0723_x_dp2.hdf5") as f:
    c = {k: f["clusters"][k][:] for k in ["cluster_id", "ra", "dec", "redshift", "richness"]}

has_src = {}
for cat in ["rubin", "sheldon"]:
    for prof in ["delta_sigma", "gt"]:
        with h5py.File(f"{OUT}/{cat}/{prof}/cluster_shear_catalogs.hdf5") as f:
            ids = f["catalog/cluster_id"][:]
            cnt = f["catalog/cluster_sample_count"][:]
        ok = set(ids[cnt > 0].tolist())
        has_src[(cat, prof)] = np.array([i in ok for i in c["cluster_id"]])
# worst case over catalogs/profiles: a cluster counts only if it has sources in all four
src_all = np.logical_and.reduce(list(has_src.values()))

STRATS = {
    "current": ([0.2, 0.4, 0.6, 0.8], [5, 20, 30, 50, 70, 200]),
    "S1": ([0.2, 0.4, 0.6, 0.8], [5, 20, 30, 50, 200]),
    "S2/S3": ([0.2, 0.4, 0.6, 0.8], [5, 20, 30, 200]),
}
extra = sys.argv[1:]  # optional "zedges;richedges" strings
for i, s in enumerate(extra):
    z, r = s.split(";")
    STRATS[f"X{i}"] = ([float(v) for v in z.split(",")], [float(v) for v in r.split(",")])

for name, (zedge, redge) in STRATS.items():
    print(f"\n== {name}  zedge={zedge}  richedge={redge}")
    print(f"{'bin':>10} {'N':>5} {'Nsrc':>5}  " + "  ".join(f"ns{n:<3}" for n in (8, 16, 32)))
    for iz in range(len(zedge) - 1):
        for ir in range(len(redge) - 1):
            m = ((c["redshift"] >= zedge[iz]) & (c["redshift"] < zedge[iz + 1])
                 & (c["richness"] >= redge[ir]) & (c["richness"] < redge[ir + 1]))
            ms = m & src_all
            regs = [np.unique(healpy.ang2pix(n, c["ra"][ms], c["dec"][ms], nest=True, lonlat=True)).size
                    for n in (8, 16, 32)]
            print(f"   z{iz}_r{ir:<3} {m.sum():5d} {ms.sum():5d}  " + "  ".join(f"{r:5d}" for r in regs))
