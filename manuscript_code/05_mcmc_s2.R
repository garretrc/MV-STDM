#setwd('')
source('manuscript_code/01_data_processing.R')

n_iter = 1500
mcmc_param <- list(n_mcmc=n_iter, n_burn=0, sigma_param=c(1,1), tau_param = c(1,1)) #should probably add more here
lk_options <- list(long=c(-180,180),lat=c(-90,90),nlevel=1, startingLevel=3, kappa=sqrt(4), alpha=c(1))

sourceCpp("dlm_code/sampling_eqs.cpp") 
source("dlm_code/mvdlm_mcmc.R")

#set holdout to NA
y_strat_missing = y_strat
y_strat_missing[t_s2, s2_uni] = NA

y_full_missing = y_full_strat
y_full_missing[t_s2, s2_multi] = NA


#Univariate
tic()
results_uni <- mvstdm_mcmc(data=y_strat_missing, coords=coords, mcmc_param = mcmc_param, lk_options = lk_options)
toc()

save(results_uni,file = 'data/results/results_s2_uni_84_95_1500.RData')
rm(results_uni)
gc()


#Univariate-RW
tic()
results_rw <- mvstdm_mcmc(data=y_strat_missing, coords=coords, mcmc_param = mcmc_param, lk_options = lk_options, fix_A = T)
toc()

save(results_rw,file = 'data/results/results_s2_rw_84_95_1500.RData')
rm(results_rw)
gc()


#Multivariate
tic()
results_multi <- mvstdm_mcmc(data=y_full_missing, coords=coords, mcmc_param = mcmc_param, lk_options = lk_options)
toc()

save(results_multi,file = 'data/results/results_s2_multi_84_95_1500.RData')
rm(results_multi)
gc()


