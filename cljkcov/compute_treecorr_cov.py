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
  
  for sim_name in yd['inputfile']['filename_sets']:
      print(sim_name)
      halo_filename = yd['inputfile']['halocat_filename_base']+sim_name+yd['inputfile']['halocat_filename_suffix']
      n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list = read_functions.read_halo_catalog(halo_filename)

      shear_filename = yd['inputfile']['shear_filename_base']+sim_name+yd['inputfile']['shear_filename_suffix']
      theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(shear_filename)

      nsky=yd['nsky']
      for nn in nsky:
          for jj in range(nn):
              ### dividing up the sky according to the phi coordinate value
              phi_lo = np.pi*2.0/nn*jj
              phi_hi = np.pi*2.0/nn*(jj+1)
              ind_cl, =np.where((phi_i > phi_lo) & (phi_i < phi_hi) & (M200b>10.0**14.0) & \
                        (z_halo > zmin) & (z_halo < zmax)  )
              ind_map, =np.where((phi > phi_lo) & (phi < phi_hi))
              print("Running on sims %s, sky division %i out of %i."%(sim_name, jj, nn))
              print("Number of clusters in this division %i"%(len(ind_cl)) )

              ## set up tree corr measurements
              NJK_patch=50
              cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
              cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers, save_patch_dir='temp')

              ## perform tree corr measurements
              if yd['run_treecorr']:
                try:
                    os.system('rm -f temp/*')
                    ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                    ng.process(cat1, cat2, low_mem=True)
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



  
# # this code uses ~9.6GB memory to store the data

# # input file

# sim_names = ['000', '064', '065', '066', '067', '068']
# #sim_names = ['067', '068']
# #sim_names = ['068']
# nsky = [2]#, 4, 8, 16]
# len_sim = len(sim_names)

# min_ang=0.1
# max_ang=10.0
# nang_bins = 10
# zmin=0.508
# zmax=0.574
# for ii in range(len_sim):
#   sim_name = sim_names[ii]
  
#   '''
#   filename = 'data/catalogs/skyhalo_nres12r'+sim_name+'.halo_z_0.508_0.574_m_13.8.fits'

#   data = fits.open(filename) 
#   z_halo = data[1].data['z_halo']
#   M200b = data[1].data['M200b']
#   theta_i = data[1].data['theta_i']
#   phi_i = data[1].data['phi_i']
#   '''
#   filename = 'data/catalogs/skyhalo_nres12r'+sim_name+'.halo'

#   with open(filename, 'rb') as f:
#     n_halo = np.fromfile(f, dtype='int32', count=1)[0]
#     ID = np.fromfile(f, dtype='int32', count=n_halo)
#     PID = np.fromfile(f, dtype='int32', count=n_halo)
#     Mvir = np.fromfile(f, dtype='float32', count=n_halo)
#     M200b = np.fromfile(f, dtype='float32', count=n_halo)
#     M200c = np.fromfile(f, dtype='float32', count=n_halo)
#     M500c = np.fromfile(f, dtype='float32', count=n_halo)
#     M2500c = np.fromfile(f, dtype='float32', count=n_halo)
#     Rvir = np.fromfile(f, dtype='float32', count=n_halo)
#     Rs = np.fromfile(f, dtype='float32', count=n_halo)
#     z_halo = np.fromfile(f, dtype='float32', count=n_halo)
#     r_halo = np.fromfile(f, dtype='float32', count=n_halo)
#     Vr = np.fromfile(f, dtype='float32', count=n_halo)
#     theta_i = np.fromfile(f, dtype='float32', count=n_halo)
#     phi_i = np.fromfile(f, dtype='float32', count=n_halo)
#     theta_s = np.fromfile(f, dtype='float32', count=n_halo)
#     phi_s = np.fromfile(f, dtype='float32', count=n_halo)
#     ipix = np.fromfile(f, dtype='int64', count=n_halo)
#     multi = np.fromfile(f, dtype='int16', count=n_halo)
#     lplane = np.fromfile(f, dtype='int16', count=n_halo)
#     hc_list = np.fromfile(f, dtype='int16', count=n_halo)


#   #plt.plot(theta_i, phi_i, 'k.')
#   #plt.show()

#   ### read in shear catalog
#   filename = 'data/allskymaps/allskymap_nres12r'+sim_name+'.zs18.mag.dat'
#   theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(filename)
#   #print(theta, phi, gamma1, gamma2, kappa, omega)
  
#   for nn in nsky:
#    for jj in range(nn):
#       print(ii, sim_name, jj)
#       ### dividing up the sky according to the phi coordinate value
#       phi_lo = np.pi*2.0/nn*jj
#       phi_hi = np.pi*2.0/nn*(jj+1)
#       ind_cl, =np.where((phi_i > phi_lo) & (phi_i < phi_hi) & (M200b>10.0**14.5) & \
#                         (z_halo > zmin) & (z_halo < zmax)  )
#       ind_map, =np.where((phi > phi_lo) & (phi < phi_hi))
#       print(len(ind_cl), len(ind_map))

#       ## set up tree corr measurements
#       NJK_patch=50
#       cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
#       cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers, save_patch_dir='temp_dir')

#       ## perform tree corr measurements
#       #try:
#       if 1:
#          os.system('rm -f temp_dir/*')
#          ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
#          ng.process(cat1, cat2, low_mem=True)
#          #ng.process(cat1, cat2)
#          #cov_jk = ng.estimate_cov('jackknife')
        
#          ng.write('data/measurements/datavector_'+sim_name+'_patchN%i_%i.txt'%(nn, jj))
#          #np.savetxt('data/measurements/cov_'+sim_name+'_patchN%i_%i.txt'%(nn, jj), cov_jk)
#       #except:
#       else:
#          print('TreeCorr not successful:', ii, sim_name, nn, jj)
#       del ng
#       del cat1
#       del cat2
