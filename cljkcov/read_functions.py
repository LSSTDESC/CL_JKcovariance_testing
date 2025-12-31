import numpy as np
import healpy as hp

# input file
#filename = 'allskymap_nres12r001.zs16.mag.dat'
def read_halo_catalog(filename = 'skyhalo_nres12r067.halo'):
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

    print('loading completed:', filename)
    return n_halo, ID, PID, Mvir, M200b, M200c, M500c, M2500c, Rvir, Rs, z_halo, r_halo, Vr, theta_i, phi_i, theta_s, phi_s, ipix, multi, lplane, hc_list



def read_map(filename = 'allskymap_nres12r000.zs18.mag.dat'):

 skip = [0, 536870908, 1073741818, 1610612728, 2147483638, 2684354547, 3221225457]
 load_blocks = [skip[i+1]-skip[i] for i in range(0, 6)]

 with open(filename, 'rb') as f:
    rec = np.fromfile(f, dtype='uint32', count=1)[0]
    nside = np.fromfile(f, dtype='int32', count=1)[0]
    npix = np.fromfile(f, dtype='int64', count=1)[0]
    rec = np.fromfile(f, dtype='uint32', count=1)[0]

    rec = np.fromfile(f, dtype='uint32', count=1)[0]

    kappa = np.array([])
    r = npix
    for i, l in enumerate(load_blocks):
        blocks = min(l, r)
        load = np.fromfile(f, dtype='float32', count=blocks)
        np.fromfile(f, dtype='uint32', count=2)
        kappa = np.append(kappa, load)
        r = r-blocks
        if r == 0:
            break
        elif r > 0 and i == len(load_blocks)-1:
            load = np.fromfile(f, dtype='float32', count=r)
            np.fromfile(f, dtype='uint32', count=2)
            kappa = np.append(kappa, load)

    gamma1 = np.array([])
    r = npix
    for i, l in enumerate(load_blocks):
        blocks = min(l, r)
        load = np.fromfile(f, dtype='float32', count=blocks)
        np.fromfile(f, dtype='uint32', count=2)
        gamma1 = np.append(gamma1, load)
        r = r-blocks
        if r == 0:
            break
        elif r > 0 and i == len(load_blocks)-1:
            load = np.fromfile(f, dtype='float32', count=r)
            np.fromfile(f, dtype='uint32', count=2)
            gamma1 = np.append(gamma1, load)

    gamma2 = np.array([])
    r = npix
    for i, l in enumerate(load_blocks):
        blocks = min(l, r)
        load = np.fromfile(f, dtype='float32', count=blocks)
        np.fromfile(f, dtype='uint32', count=2)
        gamma2 = np.append(gamma2, load)
        r = r-blocks
        if r == 0:
            break
        elif r > 0 and i == len(load_blocks)-1:
            load = np.fromfile(f, dtype='float32', count=r)
            np.fromfile(f, dtype='uint32', count=2)
            gamma2 = np.append(gamma2, load)

    omega = np.array([])
    r = npix
    for i, l in enumerate(load_blocks):
        blocks = min(l, r)
        load = np.fromfile(f, dtype='float32', count=blocks)
        np.fromfile(f, dtype='uint32', count=2)
        omega = np.append(omega, load)
        r = r-blocks
        if r == 0:
            break
        elif r > 0 and i == len(load_blocks)-1:
            load = np.fromfile(f, dtype='float32', count=r)
            np.fromfile(f, dtype='uint32', count=2)
            omega = np.append(omega, load)

 theta, phi = hp.pix2ang(nside, range(npix))
 print('loading completed:', filename)
 return theta, phi, gamma1, gamma2, kappa, omega 

 # print data
 #for i in range (npix):  
 #    print(i, kappa[i], gamma1[i], gamma2[i], omega[i])

 # example of saving data as a fits file
 #hp.fitsfunc.write_map('output.fits', kappa)
