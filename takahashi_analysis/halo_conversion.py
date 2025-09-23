import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions
import astropy.io.fits as fits
from astropy.table import Table

# this code uses ~9.6GB memory to store the data

# input file

sim_names = ['000', '064', '065', '066', '067', '068']
#sim_names = ['068']
len_sim = len(sim_names)
for ii in range(len_sim):
  sim_name = sim_names[ii]
  filename = 'skyhalo_nres12r'+sim_name+'.halo'
  
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
    
    zmin=0.508
    zmax=0.574
    Mmin=13.8
    ind, =np.where( (z_halo > zmin) & (z_halo < zmax) & \
                (M200b > 10.0**Mmin))
    print(n_halo, len(ind))
    my_table = Table({'ID':ID[ind], 'PID':PID[ind], \
                      'Mvir':Mvir[ind], 'M200b':M200b[ind], 'M200c':M200c[ind], 'M500c':M500c[ind], 'M2500c':M2500c[ind], \
                      'Rvir':Rvir[ind], 'Rs':Rs[ind], 'z_halo':z_halo[ind], 'r_halo':r_halo[ind], 'Vr':Vr[ind], \
                      'theta_i':theta_i[ind], 'phi_i':phi_i[ind], 'theta_s':theta_s[ind], 'phi_s':phi_s[ind], \
                      'ipix':ipix[ind], 'multi':multi[ind], 'lplane':lplane[ind], 'hc_list':hc_list[ind] })

    table_hdu = fits.BinTableHDU(my_table)
    hdul_table = fits.HDUList([fits.PrimaryHDU(), table_hdu]) # PrimaryHDU can be empty if no image data
    hdul_table.writeto('data/catalogs/'+filename+'_z_%.3f_%.3f_m_%.1f.fits'%(zmin, zmax, Mmin), overwrite=True)


