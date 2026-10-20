"""
Validate the sign convention needed to convert Takahashi (T17) HEALPix-native
(gamma1, gamma2) shear into the standard North/East frame CLMM expects
(coordinate_system="celestial"), before building the full cutout.

Picks the single most massive halo in realization 000's pre-made FITS cut,
pulls its nearby map pixels, and tests all 4 sign candidates through
clmm.compute_tangential_and_cross_components. The physically-correct
candidate must give clearly positive mean tangential shear (a real halo
lenses coherently) and mean cross shear consistent with ~0.
"""
import sys
import numpy as np
import healpy as hp
import astropy.io.fits as fits

sys.path.insert(0, "/global/homes/k/kabelo/CL_JKcovariance_testing/cljkcov")
import read_functions

import clmm

HALO_FITS = "/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/catalogs/skyhalo_nres12r000.halo_z_0.508_0.574_m_13.8.fits"
MAP_FILE = "/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/allskymaps/allskymap_nres12r000.zs18.mag.dat"
NSIDE = 4096
DISC_RADIUS_ARCMIN = 10.0

with fits.open(HALO_FITS) as f:
    d = f[1].data
i = np.argmax(d["M200b"])
ra_cl = np.degrees(d["phi_i"][i])
dec_cl = 90.0 - np.degrees(d["theta_i"][i])
z_cl = float(d["z_halo"][i])
m_cl = float(d["M200b"][i])
print(f"Most massive halo: M200b={m_cl:.3e} Msun, z={z_cl:.3f}, ra={ra_cl:.4f}, dec={dec_cl:.4f}")

print("Reading full-sky map (this reads the whole ~3.2GB file) ...")
theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(MAP_FILE)

vec = hp.ang2vec(np.radians(90 - dec_cl), np.radians(ra_cl))
pix = hp.query_disc(NSIDE, vec, np.radians(DISC_RADIUS_ARCMIN / 60))
ra_p = np.degrees(phi[pix])
dec_p = 90.0 - np.degrees(theta[pix])
print(f"n_pixels within {DISC_RADIUS_ARCMIN} arcmin: {len(pix)}")
print(f"mean kappa near halo (sanity check, should be > 0): {kappa[pix].mean():.4e}")

candidates = [
    ("no flip", gamma1[pix], gamma2[pix]),
    ("flip g2", gamma1[pix], -gamma2[pix]),
    ("flip g1", -gamma1[pix], gamma2[pix]),
    ("flip both", -gamma1[pix], -gamma2[pix]),
]

print()
for label, g1c, g2c in candidates:
    _, gt, gx = clmm.compute_tangential_and_cross_components(
        ra_cl, dec_cl, ra_p, dec_p, g1c, g2c,
        coordinate_system="celestial", geometry="curve",
        is_deltasigma=False, sigma_c=None, validate_input=True,
    )
    print(f"{label:12s}  mean g_t = {gt.mean():+.4e}   mean g_x = {gx.mean():+.4e}")

print()
print("The correct convention is whichever of {no flip, flip both} gives POSITIVE mean g_t.")
