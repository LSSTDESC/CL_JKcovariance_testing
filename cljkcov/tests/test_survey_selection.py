'''
Tests for select_survey_pixels in compute_treecorr_cov.py:
 - selected healpixel area matches the requested survey area to within 10 deg^2
 - all selected pixel centers lie within the spherical-cap radius of theta, phi = pi/2, pi
 - the footprint is centered on theta, phi = pi/2, pi

Run from the cljkcov/ directory:  python tests/test_survey_selection.py
'''
import os
import sys

import numpy as np
import healpy as hp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from compute_treecorr_cov import select_survey_pixels

NSIDE_SURVEY = 64  # matches the default nside_survey in select_survey_pixels


def test_survey_selection():
    center_vec = hp.ang2vec(np.pi/2, np.pi)
    for area in [320.0, 650.0, 1300.0]:
        pix_list, radius = select_survey_pixels(area)
        actual_area = len(pix_list) * hp.nside2pixarea(NSIDE_SURVEY, degrees=True)

        assert abs(actual_area - area) < 10.0, \
            'area %g deg2: selected %.2f deg2, off by more than 10 deg2' % (area, actual_area)

        vecs = np.array(hp.pix2vec(NSIDE_SURVEY, pix_list)).T
        dist = np.arccos(np.clip(vecs @ center_vec, -1.0, 1.0))
        assert dist.max() <= radius, \
            'area %g deg2: pixel center outside the cap radius' % (area)

        mean_vec = vecs.mean(axis=0)
        mean_vec = mean_vec / np.linalg.norm(mean_vec)
        offset = np.degrees(np.arccos(np.clip(mean_vec @ center_vec, -1.0, 1.0)))
        assert offset < 0.5, \
            'area %g deg2: footprint center offset from ra, dec = 0, 0 by %.2f deg' % (area, offset)

        print('area %6g deg2: %4i pixels, selected %7.2f deg2, cap radius %5.2f deg, center offset %.3f deg'
              % (area, len(pix_list), actual_area, np.degrees(radius), offset))


if __name__ == '__main__':
    test_survey_selection()
    print('All checks passed.')
