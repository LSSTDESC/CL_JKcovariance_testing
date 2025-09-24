# CL_JKcovariance_testing

This code base provides analysis tools and a pipeline for testing Jackknife covariance matrix for cluster lensing.

Galaxy clusters span scales that require high resolution in the cluster inner regions and large boxes to capture large scale structure contributions to the lensing signature.

Here, we utilize the package TreeCorr and simulation data to calculate a jackknife covariance matrix.  In a default example, we use simulations from [Takahashi et al. 2019](https://arxiv.org/pdf/1805.11629), which includes ray-traced full-sky lensing maps and halos in lightcones.  These simulations include several realizations of the universe, from which we can repeatedly measure the correlation function between cluster centers and source galaxies, and quantify the variance across all measurements.  This approach is similar to what was done in [Wu et al. 2020](https://arxiv.org/abs/1907.06611).  The variance should be comparable to the jackknife covariance.  However, with decreasing area, we expect that the jackknife covariance will no longer match the simulation variance due to extreme fluctuations in the large scale structure in any given small area of a survey.  This procedure will allow us to stress test the covariance calculation in the regimes of small area survey during the early LSST data releases.


## Getting Started

We include a default YAML file to define the simulation of interest, redshift snapshots, and area of sky that you wish to perform jackknife error estimation on.  The default file is `t19lite.yml`, which corresponds to utilizing a minimal amount of data from Takahashi et al. 2019 to test the code.

t19lite.yml
