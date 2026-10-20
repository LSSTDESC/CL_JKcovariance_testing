"""Occupied regions per bin when 30 k-means regions are defined once, on the full
cluster catalog (option b), and shared by every bin. Bins count only clusters
with background sources, as in jk_fixed.py (option a)."""
import random, numpy as np
exec(open("/global/u2/k/kabelo/CL_JKcovariance_testing/dp2/binning_test/jk_regions.py").read().split("STRATS = ")[0])
from kmeans_radec import kmeans_sample, find_nearest
random.seed(11); np.random.seed(11)
X = np.vstack((c["ra"], c["dec"])).T
# same first-pass sample size as jk_kmeans.jk_regions
nsample = min(len(X), max(2 * np.sqrt(len(X)), 10 * 30))
km = kmeans_sample(X, ncen=30, nsample=nsample, maxiter=100, verbose=0)
lab = find_nearest(X, km.centers)
sel = src_all & (c["redshift"] >= 0.2) & (c["redshift"] < 0.8) & (c["richness"] >= 5)
print("N all:", len(X), " N with sources in z/richness range:", sel.sum())
print("clusters per region (all):", np.bincount(lab, minlength=30).min(), "-", np.bincount(lab, minlength=30).max())
n_sel = np.bincount(lab[sel], minlength=30)
print("clusters per region (with sources):", n_sel.min(), "-", n_sel.max(), " empty regions:", (n_sel == 0).sum())
Z = [0.2, 0.4, 0.6, 0.8]
for name, R in {"A": [5,20,30,200], "B": [5,19,25,200], "C": [5,18,22,30,200], "D": [5,20,30,70,200]}.items():
    out = []
    for iz in range(3):
        for ir in range(len(R) - 1):
            m = sel & (c["redshift"] >= Z[iz]) & (c["redshift"] < Z[iz+1]) & (c["richness"] >= R[ir]) & (c["richness"] < R[ir+1])
            out.append(f"{m.sum()}:{np.unique(lab[m]).size}")
    print(name, " ".join(out))
