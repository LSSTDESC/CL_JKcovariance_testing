import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
import os

from use_config import load_yaml_config


import sys


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

      nsky=yd['nsky']
      for nn in nsky:
          #for jj in range(nn):
              ### dividing up the sky according to the phi coordinate value
              jj = 0
              phi_lo = np.pi*2.0/nn*jj
              phi_hi = np.pi*2.0/nn*(jj+1)
              ind_cl, =np.where((phi_i > phi_lo) & (phi_i < phi_hi) & (M200b>10.0**mmin) & (M200b<10.0**mmax) & \
                        (z_halo > zmin) & (z_halo < zmax)  )
              ind_map, =np.where((phi > phi_lo) & (phi < phi_hi))
              print("Running on sims %s, sky division %i out of %i."%(sim_name, jj, nn))
              print("Number of clusters in this division %i"%(len(ind_cl)) )

              ## set up tree corr measurements
              NJK_patch = yd['njk_patch'][nn]
              cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
              cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers, save_patch_dir='temp')

              ## perform tree corr measurements
              if yd['run_treecorr']:
                try:
                    os.system('rm -f temp/*')
                    ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                    ng.process(cat1, cat2, low_mem=low_memory_flag)

                    if yd['use_randoms']:
                        theta_randoms = np.random.rand(len(ind_cl) * 5) * (np.pi - 0)
                        phi_randoms = np.random.rand(len(ind_cl) * 5) * (phi_hi - phi_lo) + phi_lo
                        cat_random = treecorr.Catalog(x=theta_randoms, y=phi_randoms, patch_centers=cat1.patch_centers, save_patch_dir='temp')
                        rg = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                        rg.process(cat_random, cat2, low_mem=low_memory_flag)
                        ng.calculateXi(rg=rg)

                    cov_jk = ng.estimate_cov('jackknife')

                    ng.write(yd['outputpath']+'datavector_'+sim_name+'_patchN%i_%i.txt'%(nn, jj))
                    np.savetxt(yd['outputpath']+'cov_'+sim_name+'_patchN%i_%i.txt'%(nn, jj), cov_jk)
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

