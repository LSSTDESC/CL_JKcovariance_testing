import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
import os

from use_config import load_yaml_config


import sys

def select_survey_pixels(area_deg2, nside_survey=64, center_theta=np.pi/2, center_phi=np.pi):
  '''
  Select an approximately circular survey footprint of area_deg2 (in deg^2)
  around (center_theta, center_phi), as a list of healpixels at nside_survey.
  Returns (pix_list, radius) where radius is the spherical-cap radius in radians.
  '''
  area_sr = area_deg2 * (np.pi / 180.0)**2
  radius = np.arccos(1.0 - area_sr / (2.0 * np.pi))
  vec = hp.ang2vec(center_theta, center_phi)
  pix_list = hp.query_disc(nside_survey, vec, radius)
  return pix_list, radius


NSIDE_RANDOM = 2**17  # fine healpix grid for random sampling, resolution ~1.6 arcsec


def sample_survey_randoms(pix_list, n_random, nside_survey=64, nside_fine=NSIDE_RANDOM):
  '''
  Sample n_random points uniformly over the survey footprint (pix_list at
  nside_survey, RING ordering) by drawing random healpixel centers at
  nside_fine, using the NESTED parent-child hierarchy.
  Returns (theta_randoms, phi_randoms).
  '''
  coarse_nest = hp.ring2nest(nside_survey, pix_list)
  n_child = (nside_fine // nside_survey)**2
  i_coarse = np.random.randint(0, len(coarse_nest), n_random)
  i_child = np.random.randint(0, n_child, n_random)
  fine_nest = coarse_nest[i_coarse].astype(np.int64) * n_child + i_child
  return hp.pix2ang(nside_fine, fine_nest, nest=True)


def compute_treecorr_cov(yd) :
  '''
  Use info from the yaml file to compute jackknife covariances
  '''
  print(yd)

  min_ang=yd['tree_corr_angle']['min']
  max_ang=yd['tree_corr_angle']['max']
  nang_bins = yd['tree_corr_angle']['nbins']
  zmin=yd['source_redshifts']['min']
  zmax=yd['source_redshifts']['max']
  mmin=yd['source_mass']['min']
  mmax=yd['source_mass']['max']
  low_memory_flag = yd.get('low_memory_flag', True)
  
  for sim_name in yd['inputfile']['filename_sets']:
      print(sim_name)
      halo_filename = yd['inputfile']['halocat_filename_base']+sim_name+yd['inputfile']['halocat_filename_suffix']
      n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list = read_functions.read_halo_catalog(halo_filename)

      shear_filename = yd['inputfile']['shear_filename_base']+sim_name+yd['inputfile']['shear_filename_suffix']
      theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(shear_filename)

      for area in yd['survey_area']:
              ### select an ~circular survey footprint around theta, phi = pi/2, pi
              NSIDE_SURVEY = 64  # coarse healpix grid defining the survey footprint (pixel area ~0.84 deg2)
              pix_list, radius = select_survey_pixels(area)
              in_survey = np.zeros(hp.nside2npix(NSIDE_SURVEY), dtype=bool)
              in_survey[pix_list] = True
              actual_area = len(pix_list) * hp.nside2pixarea(NSIDE_SURVEY, degrees=True)

              ind_cl, =np.where(in_survey[hp.ang2pix(NSIDE_SURVEY, theta_i, phi_i)] & \
                        (M200b>10.0**mmin) & (M200b<10.0**mmax) & \
                        (z_halo > zmin) & (z_halo < zmax)  )
              ind_map, =np.where(in_survey[hp.ang2pix(NSIDE_SURVEY, theta, phi)])
              print("Running on sims %s, survey area %g deg2 (%i healpixels, %.1f deg2 selected)."%(sim_name, area, len(pix_list), actual_area))
              print("Number of clusters in this survey %i"%(len(ind_cl)) )

              ## set up tree corr measurements
              NJK_patch = yd['njk_patch'][area]
              cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
              cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers, save_patch_dir='temp')

              ## perform tree corr measurements
              if yd['run_treecorr']:
                try:
                    os.system('rm -f temp/*')
                    ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                    ng.process(cat1, cat2, low_mem=low_memory_flag)

                    if yd['use_randoms']:
                        ## draw randoms uniformly within the survey footprint from
                        ## high-resolution healpixel centers
                        n_random = len(ind_cl) * 10
                        theta_randoms, phi_randoms = sample_survey_randoms(pix_list, n_random)
                        cat_random = treecorr.Catalog(x=theta_randoms, y=phi_randoms, patch_centers=cat1.patch_centers, save_patch_dir='temp')
                        rg = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                        rg.process(cat_random, cat2, low_mem=low_memory_flag)
                        ng.calculateXi(rg=rg)

                    cov_jk = ng.estimate_cov('jackknife')

                    ng.write(yd['outputpath']+'datavector_'+sim_name+'_area%gdeg2.txt'%(area))
                    np.savetxt(yd['outputpath']+'cov_'+sim_name+'_area%gdeg2.txt'%(area), cov_jk)
                    del ng
                except:
                    print('TreeCorr not successful!')
                del cat1
                del cat2
              else:
                print('TreeCorr not run.')  


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(f"Script name: {sys.argv[0]}")
        print("Arguments received:")
        for i, arg in enumerate(sys.argv[1:]):
            print(f"  Argument {i+1}: {arg}")
            yaml_data = load_yaml_config(arg)
            compute_treecorr_cov(yaml_data)
    else:
        print("No arguments provided, using default config in configs/t19lite.yml")
        yaml_data = load_yaml_config('configs/t19lite.yml')
        compute_treecorr_cov(yaml_data)

