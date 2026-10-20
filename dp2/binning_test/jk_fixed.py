"""Occupied regions per bin when 30 k-means regions are defined once, on all
binned clusters that have background sources, and shared by every bin."""
import random, numpy as np
exec(open("/global/u2/k/kabelo/CL_JKcovariance_testing/dp2/binning_test/jk_regions.py").read().split("STRATS = ")[0])
from kmeans_radec import kmeans_sample, find_nearest
random.seed(11); np.random.seed(11)
sel = src_all & (c["redshift"] >= 0.2) & (c["redshift"] < 0.8) & (c["richness"] >= 5)
X = np.vstack((c["ra"][sel], c["dec"][sel])).T
km = kmeans_sample(X, ncen=30, nsample=len(X), maxiter=100, verbose=0)
lab = np.full(len(c["ra"]), -1); lab[sel] = find_nearest(X, km.centers)
print("clusters per region:", np.bincount(lab[sel]).min(), "-", np.bincount(lab[sel]).max())
Z = [0.2, 0.4, 0.6, 0.8]
for name, R in {"A": [5,20,30,200], "B": [5,19,25,200], "C": [5,18,22,30,200], "D": [5,20,30,70,200]}.items():
    out = []
    for iz in range(3):
        for ir in range(len(R) - 1):
            m = sel & (c["redshift"] >= Z[iz]) & (c["redshift"] < Z[iz+1]) & (c["richness"] >= R[ir]) & (c["richness"] < R[ir+1])
            out.append(f"{m.sum()}:{np.unique(lab[m]).size}")
    print(name, " ".join(out))
