import os
import sys
import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
import scipy
import gc

import healpy as hp
import clmm
from clmm import GalaxyCluster, ClusterEnsemble, GCData
from clmm import Cosmology
from clmm.support import mock_data as mock
from clmm.jk_kmeans import jk_cov
import treecorr
from astropy.table import Table, vstack


jktest_dir = '/global/homes/k/kabelo/CL_JKcovariance_testing/cljkcov'
sys.path.insert(0, jktest_dir)
from use_config import load_yaml_config
import read_functions

import warnings
warnings.simplefilter("ignore")

PLOTS_DIR = "./plots/m14_4_15/"
os.makedirs(PLOTS_DIR, exist_ok=True)


np.random.seed(11)

cosmo = ccl.Cosmology(
    Omega_c=0.265,
    Omega_b=0.0448,
    h=0.71,
    sigma8=0.8,
    n_s=0.96,
    Neff=3.04,
    m_nu=1.0e-05,
    mass_split="single",
)
cosmo_clmm = Cosmology(H0=71.0, Omega_dm0=0.265 - 0.0448, Omega_b0=0.0448, Omega_k0=0.0)
cosmo_clmm.set_be_cosmo(cosmo)
hmd_200c = ccl.halos.MassDef(200, "critical")

# loading the yaml file that was used to generate measurements.
yaml_file = os.path.join(jktest_dir, 'configs/m14_15_z0_508_0_574.yml')
yd = load_yaml_config(yaml_file)
# the simulation names and area of skies of the measurements
nskys = yd['nsky']
output_path = yd['outputpath']
sim_name = yd['inputfile']['filename_sets'][0]
print('sim name:', sim_name)
shear_filename = yd['inputfile']['shear_filename_base']+sim_name+yd['inputfile']['shear_filename_suffix']
theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(shear_filename)

halo_filename = yd['inputfile']['halocat_filename_base']+sim_name+yd['inputfile']['halocat_filename_suffix']
n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list = read_functions.read_halo_catalog(halo_filename)


min_ang=yd['tree_corr_angle']['min']
max_ang=yd['tree_corr_angle']['max']
nang_bins = yd['tree_corr_angle']['nbins']
zmin=yd['source_redshifts']['min']
zmax=yd['source_redshifts']['max']
mmin=yd['source_mass']['min']
mmax=yd['source_mass']['max']

mask = (M200c > 10**mmin) & (M200c < 10**mmax) & (z_halo > zmin) & (z_halo < zmax)
theta_cl = theta_i[mask]
phi_cl = phi_i[mask]
z_cl = z_halo[mask]
halo_id = ID[mask]

cluster_m = M200c[mask]
cluster_z = z_halo[mask]

# Concentration CCL object to compute the theoretical concentration
conc_obj = ccl.halos.ConcentrationDuffy08(mass_def=hmd_200c)
conc_list = []
for number in range(0, len(cluster_m)):
    a = 1.0 / (1.0 + (cluster_z[number]))
    # mean value of the concentration for that cluster
    lnc_mean = np.log(conc_obj(cosmo, M=(cluster_m[number]), a=a))
    conc_list.append(np.exp(lnc_mean))

conc_list = np.array(conc_list)

# Build the cluster ensemble object
gclist = []
tables = []

radius = np.radians((max_ang+10)/60)
for i, coords in enumerate(zip(theta_cl, phi_cl, z_cl)):
    vec = hp.ang2vec(coords[0], coords[1])
    # da = cosmo_clmm.eval_da(coords[2])
    pix_in_disc = hp.query_disc(4096, vec, radius)
    t, p, gam1, gam2, k = theta[pix_in_disc], phi[pix_in_disc], gamma1[pix_in_disc], gamma2[pix_in_disc], kappa[pix_in_disc]
    z = np.full_like(t, fill_value=1.218)
    t, p = np.degrees(t), np.degrees(p)
    t = 90 - t
    g1 = gam1 / (1-k)
    # above i changed coords from HEALPix --> IAU
    # according to https://healpix.sourceforge.io/html/intro_HEALPix_conventions.htm
    # flip signs of the 'U' parameter which is analogous to gamma 2
    g2 = -gam2 / (1-k)
    galcat = GCData(
        {
            "id": pix_in_disc,
            "ra": p,
            "dec": t,
            "e1": g1,
            "e2": g2,
            "z": z,
        }
    )
    
    tcl = np.degrees(coords[0])
    tcl = 90 - tcl
    pcl = np.degrees(coords[1])
    
    cl = clmm.GalaxyCluster(
        halo_id[i],
        pcl,
        tcl,
        coords[2],
        galcat=galcat
        )
    
    cl.compute_tangential_and_cross_components(
        shape_component1="e1",
        shape_component2="e2",
        tan_component="g_t",
        cross_component="g_x",
        cosmo=cosmo_clmm,
        is_deltasigma=False,
    )

    cl.compute_galaxy_weights(
        use_pdz=True,
        use_shape_noise=False,
        shape_component1="e1",
        shape_component2="e2",
        use_shape_error=False,
        weight_name="w_ls",
        cosmo=cosmo_clmm,
        is_deltasigma=False,
        add=True,
    )
    gclist.append(cl)
    tables.append(cl.galcat)
    del t, p, gam1, gam2, k, z, g1, g2

galcat_tot = vstack([Table(t) for t in tables])

bins = np.linspace(min_ang, max_ang, nang_bins+1)
ensemble_id = 1
clusterensemble = ClusterEnsemble(ensemble_id)
for cluster in gclist:
    clusterensemble.make_individual_radial_profile(
        galaxycluster=cluster,
        tan_component_in="g_t",
        cross_component_in="g_x",
        tan_component_out="g_t",
        cross_component_out="g_x",
        weights_in="w_ls",
        weights_out="W_l",
        bins=bins,
        bin_units="arcmin",
        cosmo=cosmo_clmm,
    )
    
del n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, \
    theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list, omega

clusterensemble.make_stacked_radial_profile(tan_component="g_t", cross_component="g_x", weights="W_l")
# Jackknife covariance of the stack between radial bins
jk_n_side = 16 
clusterensemble.compute_jackknife_covariance(
    n_side=jk_n_side, tan_component="g_t", cross_component="g_x"
)
njk = 20 
kmeans_gt_cov, kmeans_gx_cov, km_centers = jk_cov(clusterensemble, njk=njk)

moo = clmm.Modeling(massdef="critical", delta_mdef=200, halo_profile_model="nfw")
moo.set_cosmo(cosmo_clmm)
# Average values of mass and concentration of the ensemble to be used below
# to overplot the model on the stacked profile
moo.set_concentration(float(conc_list.mean()))
moo.set_mass(float(cluster_m.mean()))

r_stack, gt_stack, gx_stack = (clusterensemble.stacked_data[c] for c in ("radius", "g_t", "g_x"))
plt.rcParams["axes.linewidth"] = 2
fig, axs = plt.subplots(1, 2, figsize=(17, 6))
err_gt = clusterensemble.cov["tan_jk"].diagonal() ** 0.5     
err_gx = clusterensemble.cov["cross_jk"].diagonal() ** 0.5
km_err_gt = np.sqrt(np.diag(kmeans_gt_cov))
km_err_gx = np.sqrt(np.diag(kmeans_gx_cov))

axs[0].errorbar(
    r_stack,
    gt_stack,
    err_gt,
    markersize=5,
    c="r",
    fmt="o",
    capsize=10,
    elinewidth=1,
    zorder=1000,
    alpha=1,
    label="Stack",
)
axs[1].errorbar(
    r_stack,
    gx_stack,
    err_gx,
    markersize=5,
    c="r",
    fmt="o",
    capsize=10,
    elinewidth=1,
    zorder=1000,
    alpha=1,
    label="Stack",
)

axs[0].set_xscale("log")
axs[0].set_yscale("log")
axs[1].set_xscale("log")


for i in range(len(theta_cl)):
    axs[0].plot(
        clusterensemble.data["radius"][i],
        clusterensemble.data["g_t"][i],
        color="cyan",
        label="Individual",
        alpha=1,
        linewidth=1,
    )
    axs[1].plot(
        clusterensemble.data["radius"][i],
        clusterensemble.data["g_x"][i],
        color="cyan",
        label="Individual",
        alpha=1,
        linewidth=1,
    )
    if i == 0:
        axs[0].legend(frameon=False, fontsize=15)
        axs[1].legend(frameon=False, fontsize=15)

axs[0].set_xlabel("separation [arcmin]", fontsize=20)
axs[1].set_xlabel("separation [arcmin]", fontsize=20)
axs[0].tick_params(axis="both", which="major", labelsize=18)
axs[1].tick_params(axis="both", which="major", labelsize=18)
axs[0].set_ylabel(r"$g_+$", fontsize=20)
axs[1].set_ylabel(r"$g_\times$", fontsize=20)
axs[0].set_title(r"Tangential", fontsize=20)
axs[1].set_title(r"Cross", fontsize=20)

for ax in axs:
    ax.minorticks_on()
    ax.grid(lw=0.5)
    ax.grid(which="minor", lw=0.1)
fig.tight_layout()
plt.savefig(os.path.join(PLOTS_DIR, "individual_and_stacked_profiles.png"))

clusterensemble.save("ce.pkl")

# Cluster catalog
cat_cluster = treecorr.Catalog(
    ra = phi_cl,
    dec = np.pi / 2 - theta_cl,
    ra_units = "radians",
    dec_units = "radians",
    patch_centers=km_centers
)
g1 = gamma1 / (1-kappa)
# again i do the sign flip in prep for coordinate change
g2 = -gamma2/ (1-kappa)
# Source catalog
cat_source = treecorr.Catalog(
    ra = phi,
    dec = np.pi / 2 - theta,
    g1 = g1,
    g2 = g2,
    ra_units = "radians",
    dec_units = "radians",
    patch_centers=cat_cluster.patch_centers
)
del gamma1, gamma2, kappa, g1, g2, phi

ng = treecorr.NGCorrelation(nbins=nang_bins, min_sep=min_ang, max_sep=max_ang, sep_units="arcmin",bin_type = "Log")
ng.process(cat_cluster,cat_source)
theta = ng.meanr # average angular separation in each bin
cov_jk = ng.estimate_cov('jackknife')
var_jk = np.diag(cov_jk)  

da = cosmo.angular_diameter_distance(1./(1.+float(np.mean(cluster_z))))
da_clmm = cosmo_clmm.eval_da(float(np.mean(cluster_z)))
print('da:', da)
print('da_clmm:', da_clmm)

arcmin_to_Mpc =  np.pi / (60. * 180) * da_clmm

print("gt_stack:", gt_stack)
print("ng.xi:", ng.xi)

z_src = np.mean([yd['source_redshifts']['min'], yd['source_redshifts']['max']])
moo = clmm.Modeling(massdef="critical", delta_mdef=200, halo_profile_model="nfw")
moo.set_cosmo(cosmo_clmm)
# Average values of mass and concentration of the ensemble to be used below
# to overplot the model on the stacked profile
moo.set_concentration(float(conc_list.mean()))
moo.set_mass(float(cluster_m.mean()))

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 5))
ax.errorbar(
    theta,
    ng.xi,
    np.sqrt(var_jk),
    color="blue",
    marker="x",
    markersize=15,
    linestyle="",
    label="Treecorr - cluster positions x source shears",
)

ax.errorbar(
    r_stack,
    gt_stack,
    err_gt,
    markersize=5,
    c="r",
    fmt=".",
    capsize=10,
    elinewidth=1,
    zorder=1000,
    alpha=1,
    label="CLMM - stacked reduced rangential shear",
)

ax.plot(
    theta,
    moo.eval_reduced_tangential_shear(theta*arcmin_to_Mpc, float(cluster_z.mean()),
                                    z_src=z_src, z_src_info='discrete'),
    "--k",
    linewidth=1,
    label="CLMM prediction for stack 'mean' cluster",
    zorder=100,
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.axhline(0, color="k", ls=":")
ax.set_xlabel("separation [arcmin]")
ax.set_ylabel("$\\widehat{g_t}$")
ax.legend()
plt.savefig(os.path.join(PLOTS_DIR, "treecorr_vs_clmm.png"))
plt.close(fig)

fig, ax = plt.subplots()
ax.scatter(err_gt, np.sqrt(var_jk), label='healpy')
ax.scatter(km_err_gt, np.sqrt(var_jk), label='kmeans')
ax.plot([0, 0.007], [0, 0.007], "k:")
ax.set_xlabel("Error from CLMM")
ax.set_ylabel("Error from TreeCorr")
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.savefig(os.path.join(PLOTS_DIR, "error_comparison.png"))
plt.close(fig)

fig, ax = plt.subplots()
ax.scatter(km_err_gt, err_gt)
ax.plot([0, 0.007], [0, 0.007], "k:")
ax.set_xlabel("Error from kmeans")
ax.set_ylabel("Error from healpy")
ax.set_xscale('log')
ax.set_yscale('log')
plt.legend()
plt.savefig(os.path.join(PLOTS_DIR, "kmeans_healpy_error_comparison.png"))
plt.close(fig)
