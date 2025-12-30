import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
import os

# this code uses ~9.6GB memory to store the data

# input file

sim_names = ['000', '021', '023', '025', '064', '065', '066', '067', '068']
#sim_names = ['067', '068']
#sim_names = ['068']
nsky = [4, 8, 16]
#nsky = [80]#, 4, 8, 16]
len_sim = len(sim_names)

min_ang=0.1
max_ang=100.0
nang_bins = 50
zmin=0.508
zmax=0.574
for ii in range(len_sim):
  sim_name = sim_names[ii]
  
  '''
  filename = 'data/catalogs/skyhalo_nres12r'+sim_name+'.halo_z_0.508_0.574_m_13.8.fits'

  data = fits.open(filename) 
  z_halo = data[1].data['z_halo']
  M200b = data[1].data['M200b']
  theta_i = data[1].data['theta_i']
  phi_i = data[1].data['phi_i']
  '''
  filename = '/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/catalogs/skyhalo_nres12r'+sim_name+'.halo'

  with open(filename, 'rb') as f:
    n_halo = np.fromfile(f, dtype='int32', count=1)[0]
    ID = np.fromfile(f, dtype='int32', count=n_halo)
    PID = np.fromfile(f, dtype='int32', count=n_halo)
    Mvir = np.fromfile(f, dtype='float32', count=n_halo)
    M200b = np.fromfile(f, dtype='float32', count=n_halo)
    M200c = np.fromfile(f, dtype='float32', count=n_halo)
    M500c = np.fromfile(f, dtype='float32', count=n_halo)
    M2500c = np.fromfile(f, dtype='float32', count=n_halo)
    Rvir = np.fromfile(f, dtype='float32', count=n_halo)
    Rs = np.fromfile(f, dtype='float32', count=n_halo)
    z_halo = np.fromfile(f, dtype='float32', count=n_halo)
    r_halo = np.fromfile(f, dtype='float32', count=n_halo)
    Vr = np.fromfile(f, dtype='float32', count=n_halo)
    theta_i = np.fromfile(f, dtype='float32', count=n_halo)
    phi_i = np.fromfile(f, dtype='float32', count=n_halo)
    theta_s = np.fromfile(f, dtype='float32', count=n_halo)
    phi_s = np.fromfile(f, dtype='float32', count=n_halo)
    ipix = np.fromfile(f, dtype='int64', count=n_halo)
    multi = np.fromfile(f, dtype='int16', count=n_halo)
    lplane = np.fromfile(f, dtype='int16', count=n_halo)
    hc_list = np.fromfile(f, dtype='int16', count=n_halo)

  #print(np.min(theta_i), np.max(theta_i), np.min(phi_i), np.max(phi_i))
  #plt.plot(theta_i, phi_i, 'k.')
  #plt.show()

  ### read in shear catalog
  filename = '/global/cfs/cdirs/lsst/groups/CL/takahashi_sims/allskymaps/allskymap_nres12r'+sim_name+'.zs18.mag.dat'
  theta, phi, gamma1, gamma2, kappa, omega = read_functions.read_map(filename)
  #print(theta, phi, gamma1, gamma2, kappa, omega)
  
  for nn in nsky:
   #for jj in range(nn):
   for jj in np.arange(0, nn, 5):
      print(ii, sim_name, jj)
      ### dividing up the sky according to the phi coordinate value
      phi_lo = np.pi*2.0/nn*jj
      phi_hi = np.pi*2.0/nn*(jj+1)
      ind_cl, =np.where((phi_i > phi_lo) & (phi_i < phi_hi) & (M200b>10.0**14.0) & \
                        (z_halo > zmin) & (z_halo < zmax)  )
      ind_map, =np.where((phi > phi_lo) & (phi < phi_hi))
      print(len(ind_cl), len(ind_map))

      ## set up tree corr measurements
      NJK_patch=50
      cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
      cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers, save_patch_dir='temp_dir')

      ## perform tree corr measurements
      #try:
      if 1:
         os.system('rm -f temp_dir/*')
         ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
         ng.process(cat1, cat2, low_mem=True)
         #ng.process(cat1, cat2)
         cov_jk = ng.estimate_cov('jackknife')
        
         ng.write('measurements/datavector_'+sim_name+'_patchN%i_%i.txt'%(nn, jj))
         np.savetxt('measurements/cov_'+sim_name+'_patchN%i_%i.txt'%(nn, jj), cov_jk)
      #except:
      else:
         print('TreeCorr not successful:', ii, sim_name, nn, jj)
      del ng
      del cat1
      del cat2
