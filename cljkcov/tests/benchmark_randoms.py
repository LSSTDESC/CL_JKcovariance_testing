'''
Benchmark of two ways to draw randoms uniformly within the survey footprint:
 - rejection: uniform points in a bounding box around the cap, kept if their
   nside-64 healpixel is in the footprint (the previous implementation)
 - healpix: draw a random coarse footprint pixel and a random NESTED child
   pixel at high resolution, use the child pixel center coordinates

Run from the cljkcov/ directory:  python tests/benchmark_randoms.py
'''
import os
import sys
import time

import numpy as np
import healpy as hp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compute_treecorr_cov import select_survey_pixels

NSIDE_SURVEY = 64
NSIDE_RANDOM = 2**17


def rejection_randoms(in_survey, radius, n_random, nside_survey=NSIDE_SURVEY):
    margin = radius + 2.0 * hp.nside2resol(nside_survey)
    costh_lo = np.cos(np.pi/2 + margin)
    costh_hi = np.cos(np.pi/2 - margin)
    theta_randoms = np.empty(0)
    phi_randoms = np.empty(0)
    while len(theta_randoms) < n_random:
        nbatch = 2 * (n_random - len(theta_randoms)) + 100
        th = np.arccos(np.random.rand(nbatch) * (costh_hi - costh_lo) + costh_lo)
        ph = np.pi + (np.random.rand(nbatch) * 2.0 - 1.0) * margin
        keep = in_survey[hp.ang2pix(nside_survey, th, ph)]
        theta_randoms = np.append(theta_randoms, th[keep])
        phi_randoms = np.append(phi_randoms, ph[keep])
    return theta_randoms[:n_random], phi_randoms[:n_random]


def healpix_randoms(pix_list, n_random, nside_survey=NSIDE_SURVEY, nside_fine=NSIDE_RANDOM):
    coarse_nest = hp.ring2nest(nside_survey, pix_list)
    n_child = (nside_fine // nside_survey)**2
    i_coarse = np.random.randint(0, len(coarse_nest), n_random)
    i_child = np.random.randint(0, n_child, n_random)
    fine_nest = coarse_nest[i_coarse].astype(np.int64) * n_child + i_child
    return hp.pix2ang(nside_fine, fine_nest, nest=True)


def best_of(func, nrepeat=5):
    times = []
    for _ in range(nrepeat):
        t0 = time.perf_counter()
        func()
        times.append(time.perf_counter() - t0)
    return min(times)


if __name__ == '__main__':
    area = 320.0
    pix_list, radius = select_survey_pixels(area)
    in_survey = np.zeros(hp.nside2npix(NSIDE_SURVEY), dtype=bool)
    in_survey[pix_list] = True

    ## containment check for both samplers
    for name, sampler in [('rejection', lambda n: rejection_randoms(in_survey, radius, n)),
                          ('healpix', lambda n: healpix_randoms(pix_list, n))]:
        th, ph = sampler(50000)
        assert len(th) == 50000
        assert in_survey[hp.ang2pix(NSIDE_SURVEY, th, ph)].all(), name + ': points outside footprint'
    print('Both samplers keep all points inside the footprint.\n')

    print('%10s %14s %14s' % ('n_random', 'rejection [s]', 'healpix [s]'))
    for n_random in [3000, 100000, 1000000]:
        t_rej = best_of(lambda: rejection_randoms(in_survey, radius, n_random))
        t_hpx = best_of(lambda: healpix_randoms(pix_list, n_random))
        print('%10i %14.5f %14.5f' % (n_random, t_rej, t_hpx))
