#' Bayesian Estimation Multivariate Spatial Dynamic Linear Model
#'
#' This function computes the posterior of the parameters of a multivariate space-time dynamic linear model (MVSTDM)
#' using Gibbs sampling. Relies on LatticeKrig package.
#'
#' Branching models specify gamma priors for mu, alpha and beta parameters.
#'
#' @param data - T x n x m array of training data with T time points, n spatial locations and m variables
#' @param coords - n x 2 matrix of location coordinates [lon, lat] or [x, y]
#' @param m - number of variables
#' @param lk_options - list of parameters for LKrigSetup()
#' @param param_init - list of parameters of initial guess (default = NULL, will start with generic values)
#' @param fix_A - boolean value to determine if A matrix is resampled at every mcmc iteration, TRUE = A is fixed and not resampled (default=FALSE), if param_init is also NULL, this specifies a random walk
#' @param lower_triangular_A - boolean value that specifies if A should be lower triangular (default=FALSE)
#' @param mcmc_param - list of mcmc parameters
#' @return A List containing the mcmc samples
#' @export
mvstdm_mcmc <- function(data, coords, lk_options=NULL, param_init=NULL, fix_A=FALSE, lower_triangular_A=FALSE, mcmc_param=NULL){
    n = dim(coords)[1]
    nt = dim(data)[1]

    #reformatting data array into T x [n x m] matrix y
    if(length(dim(data))>2){
        y = matrix(nrow=dim(data)[1],ncol=prod(dim(data)[-1]))
        m = dim(data)[3]
        for(i in 1:m){
            y[,((i-1)*n+1):((i-1)*n+n)] = data[,,i]
        }
    }else{
        m = ncol(data)/n
        y = data
    }

    #setting default mcmc_param
    if (is.null(mcmc_param)){
      mcmc_param <- list(n_mcmc=2000, n_burn=500, sigma_param=c(1,1), tau_param = c(1,1))
    }

    #setting up basis function construction
    if (is.null(lk_options)){
      lk_options <- list(long=c(-180,180),lat=c(-90,90),nlevel=1, startingLevel=3, kappa=sqrt(4), alpha=c(1))
    }

    LKinfo = LKrigSetup(x=data.frame(long=lk_options$long,lat=lk_options$lat),
                    nlevel=lk_options$nlevel,
                    LKGeometry = 'LKSphere',
                    a.wght = 1+lk_options$kappa^2,
                    startingLevel=lk_options$startingLevel,
                    alpha=lk_options$alpha)
    lattice_basis = LKrig.basis(x1 = coords,LKinfo = LKinfo)
    Phi = as.matrix(lattice_basis)
    B = LKrig.precision(LKinfo,return.B = T)

    if (is.null(param_init)){
      k = dim(Phi)[2]
      param_init <- list(alpha = matrix(rnorm((nt+1)*m*k,0,1),nt+1,k*m), 
                         A = diag(1,k*m), 
                         sigma = matrix(1,nrow=nt,ncol=m), 
                         tau = rep(1,m))
    }

    #MCMC sampler in C++
    post_samps <- mcmc_dlm(y, as.matrix(B), Phi, 
                           n_mcmc = mcmc_param$n_mcmc, n_burn = mcmc_param$n_burn,
                           fix_A = fix_A,
                           lower_triangular_A = lower_triangular_A,
                           alpha_init = param_init$alpha, A_init = param_init$A, 
                           sig_init = param_init$sigma, tau_init = param_init$tau,
                           a_sigma = mcmc_param$sigma_param[1], b_sigma = mcmc_param$sigma_param[2], 
                           a_tau = mcmc_param$tau_param[1], b_tau = mcmc_param$tau_param[2])
    
    #Burn in
    post_samps$sigma <- post_samps$sigma[,,(mcmc_param$n_burn + 1):mcmc_param$n_mcmc]
    post_samps$tau <- post_samps$tau[,(mcmc_param$n_burn + 1):mcmc_param$n_mcmc]
    post_samps$A <- post_samps$A[,(mcmc_param$n_burn + 1):mcmc_param$n_mcmc]
    post_samps$alpha <- post_samps$alpha[,,(mcmc_param$n_burn + 1):mcmc_param$n_mcmc]
    
    return(post_samps)
}
