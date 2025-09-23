import numpy as np

# this code uses ~9.6GB memory to store the data

# input file
filename = 'skyhalo_nres12r000.halo'

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

print(n_halo)

for k in range (n_halo):  
    hc_catalog='hc_'+str(450*(hc_list[k]//100+1))+'Mpc_r'+format((hc_list[k]%100)//10,'03')+'_'+format(hc_list[k]%10,'03')+'.list' # <- name of original rockstar halo catalog 
    print(k, ID[k], PID[k], Mvir[k], M200b[k], M200c[k], M500c[k], M2500c[k], Rvir[k], Rs[k], z_halo[k], r_halo[k], Vr[k], theta_i[k], phi_i[k], theta_s[k], phi_s[k], ipix[k], multi[k], lplane[k], hc_catalog)


