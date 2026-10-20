import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
import os
import gc
import time

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
              ind_cl, =np.where((phi_i > phi_lo) & (phi_i < phi_hi) & (M200b>10.0**mmin) & (M200b<10.0**mmax) & \
                        (z_halo > zmin) & (z_halo < zmax)  )
              ind_map, =np.where((phi > phi_lo) & (phi < phi_hi))
              print("Running on sims %s, sky division %i out of %i."%(sim_name, jj, nn))
              print("Number of clusters in this division %i"%(len(ind_cl)) )

              ## set up tree corr measurements
              NJK_patch=50
              cat1 = treecorr.Catalog(x=theta_i[ind_cl], y=phi_i[ind_cl], npatch=NJK_patch)
              cat2 = treecorr.Catalog(x=theta[ind_map], y=phi[ind_map], g1 = gamma1[ind_map], g2=gamma2[ind_map], patch_centers=cat1.patch_centers)#, save_patch_dir='temp')

              ## perform tree corr measurements
              if yd['run_treecorr']:
                ng = treecorr.NGCorrelation(min_sep=min_ang, max_sep=max_ang, nbins=nang_bins, sep_units='arcmin', bin_type='Log')
                start = time.time()
                ng.process(cat1, cat2, low_mem=False)
                end = time.time()
                print(f'ng procees took {end - start:.1f}s')
                cov_jk = ng.estimate_cov('jackknife')

                ng.write(yd['outputpath']+'datavector_'+sim_name+'_patchN%i_%i.txt'%(nn, jj))
                np.savetxt(yd['outputpath']+'cov_'+sim_name+'_patchN%i_%i.txt'%(nn, jj), cov_jk)
                del ng
                del cat1
                del cat2
              else:
                print('TreeCorr not run.')
                  
      del n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list
      del theta, phi, gamma1, gamma2, kappa, omega
      gc.collect()


def compute_treecorr_txpipe_bins(yd) :
  '''
  Stacked reduced shear + cluster-based healpix jackknife on a TXPipe-format cutout,
  using the same (zbin, richbin) selection and physical-Mpc radial bins as
  TXPipe's CLClusterEnsembleProfiles, so the results can be compared bin by bin.
  '''
  import pickle
  from astropy.cosmology import FlatLambdaCDM

  print(yd)
  clusters, sources = read_functions.read_txpipe_cutout(yd['inputfile']['cluster_catalog'],
                                                        yd['inputfile']['shear_catalog'])
  cosmo = FlatLambdaCDM(H0=yd['cosmology']['H0'], Om0=yd['cosmology']['Om0'], Ob0=yd['cosmology']['Ob0'],
                        Tcmb0=2.725, Neff=3.046, m_nu=0.0)
  rbins = yd['radial_bins']
  nside = yd['jackknife_nside']
  zedge = np.array(yd['zedge'])
  massedge = np.array(yd['massedge'])
  outpath = yd['outputpath']
  os.makedirs(outpath, exist_ok=True)

  ## halo selection: z < halo_zmax, then TXPipe's redshift/richness(=log10 M200b) bins
  sel = clusters['redshift'] < yd['halo_zmax']
  print("Halos with z < %.2f: %i of %i"%(yd['halo_zmax'], sel.sum(), len(sel)))
  izbin = np.digitize(clusters['redshift'], zedge) - 1
  imbin = np.digitize(clusters['richness'], massedge) - 1

  ## TXPipe source cuts, trivially satisfied for these catalogs; check rather than apply
  zs_min = sources['redshift_true'].min()
  assert zs_min > zedge[-1] + yd['delta_z'], "Sources do not all pass TXPipe's delta_z cut"
  theta_max = np.degrees(rbins['max'] / cosmo.angular_diameter_distance(max(zedge[0], clusters['redshift'][sel].min())).value)*60
  if theta_max > yd['max_angle_arcmin']:
      print("WARNING: r_max subtends %.1f arcmin > TXPipe max_angle %.1f arcmin for the lowest-z halos"%(theta_max, yd['max_angle_arcmin']))

  ## the source r is not used by the Rlens metric, but 3D positions are required
  src_r = cosmo.angular_diameter_distance(np.median(sources['redshift_true'])).value * np.ones(len(sources['ra']))
  cat_src = treecorr.Catalog(ra=sources['ra'], dec=sources['dec'], r=src_r, ra_units='deg', dec_units='deg',
                             g1=sources['g1'], g2=sources['g2'], w=sources['weight'], flip_g2=yd['flip_g2'])

  def run_ng(ind):
      ng = treecorr.NGCorrelation(min_sep=rbins['min'], max_sep=rbins['max'], nbins=rbins['nbins'],
                                  bin_type='Log', metric='Rlens', min_rpar=0)
      cat_lens = treecorr.Catalog(ra=clusters['ra'][ind], dec=clusters['dec'][ind],
                                  r=cosmo.angular_diameter_distance(clusters['redshift'][ind]).value,
                                  ra_units='deg', dec_units='deg')
      ng.process(cat_lens, cat_src)
      return ng

  results = {}
  for iz in range(len(zedge)-1):
      for im in range(len(massedge)-1):
          name = 'bin_zbin_%i_richbin_%i'%(iz, im)
          ind, = np.where(sel & (izbin == iz) & (imbin == im))
          print("%s: %i clusters"%(name, len(ind)))
          if len(ind) < 2 or not yd['run_treecorr']:
              continue

          start = time.time()
          ng = run_ng(ind)
          ng.write(outpath+'datavector_'+name+'.txt')

          ## cluster-based jackknife over occupied healpix pixels (same as CLMM compute_jackknife_covariance)
          pix = hp.ang2pix(nside, clusters['ra'][ind], clusters['dec'][ind], nest=True, lonlat=True)
          upix = np.unique(pix)
          njk = len(upix)
          ng_jk = [run_ng(ind[pix != p]) for p in upix]
          gt_jk = np.array([n.xi for n in ng_jk])
          gx_jk = np.array([n.xi_im for n in ng_jk])
          del ng_jk
          dgt = gt_jk - gt_jk.mean(axis=0)
          dgx = gx_jk - gx_jk.mean(axis=0)
          ## CLMM normalisation: (N-1)^2/N * np.cov(ddof=0) = (N-1)^2/N^2 * sum,
          ## i.e. (N-1)/N times the standard (N-1)/N * sum jackknife estimate (kept as cov_gt_std)
          cov_gt = (njk-1)**2/njk**2 * dgt.T @ dgt
          cov_gx = (njk-1)**2/njk**2 * dgx.T @ dgx
          cov_gt_std = njk/(njk-1) * cov_gt
          end = time.time()
          print('%s: njk = %i, took %.1fs'%(name, njk, end - start))

          np.savetxt(outpath+'cov_'+name+'.txt', cov_gt)
          results[name] = {'n_cl': len(ind), 'cluster_id': clusters['cluster_id'][ind],
                           'meanr': ng.meanr, 'rnom': ng.rnom, 'npairs': ng.npairs, 'weight': ng.weight,
                           'gt': ng.xi, 'gx': ng.xi_im, 'varxi': ng.varxi,
                           'njk': njk, 'cov_gt': cov_gt, 'cov_gx': cov_gx, 'cov_gt_std': cov_gt_std}
          del ng
          gc.collect()

  with open(outpath+'treecorr_profiles.pkl', 'wb') as f:
      pickle.dump(results, f)
  print('Wrote '+outpath+'treecorr_profiles.pkl')
  return results


if __name__ == "__main__":
    if len(sys.argv) > 1:
        print(f"Script name: {sys.argv[0]}")
        print("Arguments received:")
        for i, arg in enumerate(sys.argv[1:]):
            print(f"  Argument {i+1}: {arg}")
            yaml_data = load_yaml_config(arg)
            if yaml_data.get('input_format') == 'txpipe_hdf5':
                compute_treecorr_txpipe_bins(yaml_data)
            else:
                compute_treecorr_cov(yaml_data)
    else:
        print("No arguments provided, using default config in configs/t19lite.yml")
        yaml_data = load_yaml_config('configs/t19lite.yml')
        compute_treecorr_cov(yaml_data)

