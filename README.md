# Multivariate Space-Time Dynamic Linear Model (MVSTDM)

## Authors: Robert Garret, Lyndsay Shand, J. Gabriel Huerta

## This code repository includes code to implement the multivariate space-time dynamic linear model and to recreate the results highlighted in our manuscript. 

### Funding Statement: Sandia National Laboratories is a multimission laboratory managed and operated by National Technology & Engineering Solutions of Sandia, LLC, a wholly owned subsidiary of Honeywell International Inc., for the U.S. Department of Energy’s National Nuclear Security Administration under contract DE-NA0003525.

# dlm_code: repositry for model code
# manuscript_code: run scripts to recreate results from manuscript

## 01_data_processing

Data processing script to obtain the multivariate spatiotemporal field matrices for each pathway. Also provides codes for spatial plots and indices for the space and time locations of each holdout set scenario. Scenario 1, or s1, is the spatial block holdout, and Scenario 2, or s2, is the random location holdout.

## 02_simulation

This file contains the data generation, MCMC, and figures for the simulated data example.

## 03, 04, 05_mcmc...

These files contain the code to run the DLM sampler, contained in spatial_stats/dlm_code, and produce the results for the full data and both s1 and s2. Results are saved in cee drive, cldera/obs_thrust/dm_results.

## 06, 07, 08 plots

These files contain the figures for the results sections and a few additional explorations of the posterior samples.
