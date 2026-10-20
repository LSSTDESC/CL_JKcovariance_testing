"""
Build a ~100 deg^2 synthetic cluster-shear dataset from the Takahashi et al.
2017 (T17) full-sky ray-tracing simulation, in the same schema
CLClusterShearCatalogs/CLClusterEnsembleProfiles expect (mirrors
make_dc2_cutout_100deg2.py's conventions closely).

Box: ra in [10, 20] deg, dec in [-5, 5] deg (10x10 = 100 deg^2), placed near
the equator and away from the phi=0/2pi seam so a plain box mask is safe.
Realization: 000.

- Raw binary halo catalog and full-sky healpix (RING, nside=4096) kappa/
gamma1/gamma2/omega map, parsed via the existing, unmodified
read_functions.read_halo_catalog / read_map (see
~/CL_JKcovariance_testing/cljkcov/read_functions.py).
- phi = RA directly (no rotation); dec = 90 - colatitude_deg.
- Shear sign convention CONFIRMED EMPIRICALLY (validate_t17_shear_
convention.py): CLMM's coordinate_system="celestial" frame requires
g1 = -gamma1_healpix, g2 = -gamma2_healpix (a single-axis parity flip on
BOTH components -- the "no flip" and partial-flip candidates were all
ruled out by the validation script's clean positive-g_t test).
- Source plane is a single fixed redshift: zs18 -> z_source = 1.218
(confirmed 3 independent ways: public T17 zs-tag table, first-principles
recomputation from T17's 150 h^-1 Mpc comoving source-plane spacing, and
the cljkcov config's own self-flagged stale-placeholder comment).
- Every source trivially passes the background cut (z_source=1.218 well
behind any cluster z<=0.8 + delta_z=0.2) -- expected, not a bug: isolates
the radius/angle/cosmology selection effect from photo-z-driven source-
count noise, since there IS no photo-z scatter in this dataset.
- richness proxy = log10(M200b) (compresses the ~2-dex mass range into
richness-like bin edges; raw mass wouldn't map sensibly onto richedge
cuts).

Run via sbatch -- read_halo_catalog/read_map load their full files into memory (~10GB+ for halos, map is way bigger)
"""
import os
import sys
import numpy as np
import h5py

sys.path.insert(0, "/global/homes/k/kabelo/CL_JKcovariance_testing/cljkcov")
import read_functions

RA_MIN, RA_MAX = 10.0, 20.0
DEC_MIN, DEC_MAX = -5.0, 5.0

Z_MIN, Z_MAX = 0.2, 0.8  # for zsource=1.218 the lensing efficiency peaks aroud z=0.5 so this is fine
LOG_MASS_MIN = 13.8  # log10(M200b / Msun)

Z_SOURCE = 1.218

HALO_FILE = "/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/catalogs/skyhalo_nres12r000.halo"
MAP_FILE = "/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/allskymaps/allskymap_nres12r000.zs18.mag.dat"

OUT_DIR = "/pscratch/sd/k/kabelo/txp_cluster_compare/t17_cutout_100deg2"


def in_box(ra, dec):
    return (ra >= RA_MIN) & (ra <= RA_MAX) & (dec >= DEC_MIN) & (dec <= DEC_MAX)


def make_cluster_cutout():
    print("Reading halo catalog (full ~10GB read) ...")
    (n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo,
     Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list) = \
        read_functions.read_halo_catalog(HALO_FILE)

    ra = np.degrees(phi_i.astype(np.float64))
    dec = 90.0 - np.degrees(theta_i.astype(np.float64))

    box_mask = in_box(ra, dec)
    z_mask = (z_halo >= Z_MIN) & (z_halo <= Z_MAX)
    log_mass = np.log10(M200b.astype(np.float64), out=np.full(len(M200b), -np.inf, dtype=np.float64), where=M200b > 0)
    mass_mask = log_mass >= LOG_MASS_MIN
    mask = box_mask & z_mask & mass_mask
    n = int(mask.sum())
    print(f"  {n} halos in box (ra {RA_MIN}-{RA_MAX}, dec {DEC_MIN}-{DEC_MAX}, "
          f"z {Z_MIN}-{Z_MAX}, log10(M200b)>={LOG_MASS_MIN})")
    if n > 0:
        lm = log_mass[mask]
        print(f"  log10(M200b) percentiles [10,50,90]: {np.percentile(lm, [10, 50, 90])}")
        print(f"  z_halo percentiles [10,50,90]: {np.percentile(z_halo[mask], [10, 50, 90])}")

    with h5py.File(f"{OUT_DIR}/cluster_catalog.hdf5", "w") as fout:
        g = fout.create_group("clusters")
        g.create_dataset("cluster_id", data=ID[mask].astype(np.int64))
        g.create_dataset("ra", data=ra[mask])
        g.create_dataset("dec", data=dec[mask])
        g.create_dataset("redshift", data=z_halo[mask].astype(np.float32))
        g.create_dataset("redshift_err", data=np.zeros(n, dtype=np.float32))
        g.create_dataset("richness", data=log_mass[mask].astype(np.float32))
        g.create_dataset("richness_err", data=np.zeros(n, dtype=np.float32))
        g.create_dataset("scaleval", data=np.ones(n, dtype=np.float32))
    print(f"  wrote {n} clusters")
    return n


def make_shear_cutout():
    print("Reading full-sky map (full ~3.2GB read) ...")
    theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(MAP_FILE)

    ra = np.degrees(phi)
    dec = 90.0 - np.degrees(theta)
    mask = in_box(ra, dec)
    n = int(mask.sum())
    print(f"  {n} map pixels in box")

    g1 = -gamma1[mask]
    g2 = -gamma2[mask]

    with h5py.File(f"{OUT_DIR}/shear_catalog.hdf5", "w") as fout:
        g_shear = fout.create_group("shear")
        g_shear.attrs["catalog_type"] = "metadetect"
        g_ns = g_shear.create_group("ns")
        g_ns.create_dataset("ra", data=ra[mask])
        g_ns.create_dataset("dec", data=dec[mask])
        g_ns.create_dataset("g1", data=g1)
        g_ns.create_dataset("g2", data=g2)
        g_ns.create_dataset("weight", data=np.ones(n, dtype=np.float64))
        g_ns.create_dataset("redshift_true", data=np.full(n, Z_SOURCE, dtype=np.float64))
    print(f"  wrote {n} rows under shear/ns")

    with h5py.File(f"{OUT_DIR}/shear_tomography_catalog.hdf5", "w") as fout:
        g_tomo = fout.create_group("tomography")
        g_tomo.attrs["catalog_type"] = "metadetect"
        g_tomo.attrs["nbin_source"] = 1
        g_tomo.attrs["source_zmin_0"] = Z_SOURCE - 1e-6
        g_tomo.attrs["source_zmax_0"] = Z_SOURCE + 1e-6
        g_tomo.create_dataset("bin", data=np.zeros(n, dtype=np.int32))
        g_tomo.create_dataset("source_bin", data=np.zeros(n, dtype=np.int32))
        for name in ("mean_e1", "mean_e1_2d", "mean_e2", "mean_e2_2d", "sigma_e", "sigma_e_2d"):
            g_tomo.create_dataset(name, data=np.zeros(1, dtype=np.float32))
        for name in ("N_eff", "N_eff_2d"):
            g_tomo.create_dataset(name, data=np.array([float(n)], dtype=np.float32))
        for name in ("source_counts", "source_counts_2d"):
            g_tomo.create_dataset(name, data=np.array([n], dtype=np.int32))

        g_resp = fout.create_group("response")
        g_resp.create_dataset("R", data=np.eye(2, dtype=np.float32)[None, :, :])
        g_resp.create_dataset("R_2d", data=np.eye(2, dtype=np.float32))

        g_counts = fout.create_group("counts")  # required by MetaDetectCalibrator.load
        for name in ("mean_e1", "mean_e2", "mean_e1_2d", "mean_e2_2d"):
            g_counts.create_dataset(name, data=np.zeros(1, dtype=np.float32))
    print(f"  wrote tomography/{{bin,source_bin}} for {n} rows")

    with h5py.File(f"{OUT_DIR}/source_photoz_pdfs.hdf5", "w") as fout:
        g_meta = fout.create_group("meta")
        g_meta.create_dataset("xvals", data=np.array([[Z_SOURCE]]))
        g_meta.create_dataset("pdf_name", data=[b"delta"])
        g_meta.create_dataset("pdf_version", data=[0])

        g_ancil = fout.create_group("ancil")
        g_ancil.create_dataset("zmode", data=np.full(n, Z_SOURCE, dtype=np.float64))
        g_ancil.create_dataset("zmean", data=np.full(n, Z_SOURCE, dtype=np.float64))

        g_data = fout.create_group("data")
        g_data.create_dataset("yvals", data=np.ones((n, 1), dtype=np.float64))
    print(f"  wrote ancil/{{zmode,zmean}} and data/yvals for {n} rows")

    return n


def verify():
    print("Verifying outputs ...")
    with h5py.File(f"{OUT_DIR}/cluster_catalog.hdf5") as f:
        print("  cluster_catalog:", f["clusters/ra"].shape)
    with h5py.File(f"{OUT_DIR}/shear_catalog.hdf5") as f:
        print("  shear_catalog:", f["shear/ns/ra"].shape, "catalog_type =", f["shear"].attrs["catalog_type"])
    with h5py.File(f"{OUT_DIR}/shear_tomography_catalog.hdf5") as f:
        print("  tomography:", f["tomography/bin"].shape, "counts/mean_e1 =", f["counts/mean_e1"][:])
    with h5py.File(f"{OUT_DIR}/source_photoz_pdfs.hdf5") as f:
        print("  photoz:", f["data/yvals"].shape)


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    make_cluster_cutout()
    make_shear_cutout()
    verify()
    print("Done.")
