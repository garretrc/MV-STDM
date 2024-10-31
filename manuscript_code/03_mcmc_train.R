#setwd('')
source('manuscript_code/01_data_processing.R')

n_iter = 1500
mcmc_param <- list(n_mcmc=n_iter, n_burn=0, sigma_param=c(1,1), tau_param = c(1,1)) #should probably add more here
lk_options <- list(long=c(-180,180),lat=c(-90,90),nlevel=1, startingLevel=3, kappa=sqrt(4), alpha=c(1))

sourceCpp("dlm_code/sampling_eqs.cpp") 
source("dlm_code/mvdlm_mcmc.R")


#Multivariate full data run
tic()
results_multi <- mvstdm_mcmc(data=y_full_strat, coords=coords, mcmc_param = mcmc_param, lk_options = lk_options)
#results_multi_2 <- mvstdm_mcmc(data=y_full_strat, coords=coords, mcmc_param = mcmc_param, lk_options = lk_options)
toc()

save(results_multi,file = "data/results/results_strat_full_multi_84_95_1500_c1.RData")

