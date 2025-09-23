import numpy as np
import treecorr
import healpy as hp
import matplotlib.pyplot as plt
import read_functions

fig, axs = plt.subplots(2)

sim_names = ['000', '064', '065', '067', '068']
#sim_names = ['068']
len_sim = len(sim_names)
nsky = 2
for ii in range(len_sim):
  for jj in range(nsky):
      sim_name = sim_names[ii]
      dat= np.genfromtxt('data/measurements/datavector_'+sim_name+'_patchN%i_%i.txt'%(nsky, jj))
      axs[0].plot(dat[:, 0], dat[:, 3], '--', label = sim_name)

      if ii == 0 and jj == 0:
         meas_rec = dat[:, 3]
      else:
         meas_rec = np.vstack((meas_rec, dat[:, 3]))
      print(ii, jj, np.max(dat[:, 3]))
print(np.shape(meas_rec))
cov_jk = np.genfromtxt('data/measurements/cov_000_patchN2_0.txt')
x_val = dat[:, 0]
y_val=np.mean(meas_rec, axis=0)
y_sim_err=np.std(meas_rec, axis=0)
y_jk_err=np.sqrt(np.diag(cov_jk))
#axs[0].errorbar(x_val, y_val, yerr=y_sim_err, fmt='bs', label='Simulation Variation')
#axs[0].errorbar(x_val, y_val, yerr=y_jk_err, fmt='ro', label = 'Jackknife uncertainty')

axs[1].errorbar(x_val, y_sim_err, fmt='bs', label='Simulation Variation')
axs[1].errorbar(x_val, y_jk_err, fmt='ro', label = 'Jackknife uncertainty')


print(np.shape(meas_rec))
axs[0].legend(loc=0)
axs[0].set_xscale('log')
axs[0].set_yscale('log')
axs[0].set_xlim(0.7, 10)
axs[0].set_xlabel(r'$r$ [arcmin]')
axs[0].set_ylabel(r'$\gamma_T$')


axs[1].set_xscale('log')
axs[1].set_xlim(0.7, 10)
axs[1].set_xlabel(r'$r$ [arcmin]')
axs[1].set_ylabel(r'Uncertainty on $\gamma_T$')
plt.show()

