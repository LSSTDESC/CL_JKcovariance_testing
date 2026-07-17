'''
Tests for sample_survey_randoms in compute_treecorr_cov.py:
 - the requested number of randoms is returned
 - all randoms fall inside the survey footprint
 - randoms are uniform across the footprint pixels (Poisson-consistent counts)

Run from the cljkcov/ directory:  python tests/test_survey_randoms.py
'''
import os
import sys

import numpy as np
import healpy as hp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compute_treecorr_cov import select_survey_pixels, sample_survey_randoms

NSIDE_SURVEY = 64  # matches the default nside_survey in select_survey_pixels


def test_survey_randoms():
    n_random = 50000
    pix_list, radius = select_survey_pixels(320.0)
    in_survey = np.zeros(hp.nside2npix(NSIDE_SURVEY), dtype=bool)
    in_survey[pix_list] = True

    theta_r, phi_r = sample_survey_randoms(pix_list, n_random)
    assert len(theta_r) == n_random

    pix_r = hp.ang2pix(NSIDE_SURVEY, theta_r, phi_r)
    assert in_survey[pix_r].all(), 'randoms fall outside the survey footprint'

    ## uniformity: counts per footprint pixel consistent with Poisson
    counts = np.bincount(pix_r, minlength=hp.nside2npix(NSIDE_SURVEY))[pix_list]
    mean = n_random / float(len(pix_list))
    assert counts.min() > 0, 'empty footprint pixel'
    assert np.abs(counts - mean).max() < 5.0 * np.sqrt(mean), \
        'random counts per pixel deviate more than 5 sigma from uniform'

    print('randoms: %i points, all inside footprint; counts per pixel %i-%i (mean %.1f)'
          % (n_random, counts.min(), counts.max(), mean))


if __name__ == '__main__':
    test_survey_randoms()
    print('All checks passed.')
