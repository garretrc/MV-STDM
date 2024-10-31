#Customized Versions of DLM Filtering and Sampling

#library(dlm)

#### 0. AR1 Gibbs sampler from Gabriel Huerta ####

AR1_gibbs = function(alpha){ 
  # Code from Gabriel Huerta
  
  n = length(alpha)
  f1 = alpha[1:(n-1)]
  y = alpha[2:n]
  FF.t <- sum(f1^2)
  FF.inv <- 1/FF.t
  a_hat <- sum(f1*y)*FF.inv
  R <- (y - f1*a_hat)
  ssq = sum(R*R)
  
  # sampling innovation variance
  tau2 = 1/rgamma(1,shape=n-1,rate=ssq/2)
  
  # sampling AR1 coefficient
  a = rnorm(1, mean = a_hat, sd = sqrt(tau2*FF.inv))
  
  return(c(a,tau2))
}

#### 3. SAR model for V ####

dlmFilterSAR <- function(y, mod, simplify = FALSE)
{
  ## Note: V must be nonsingular. It will be forced to be so,
  ## with a warning, otherwise (time-invariant case only).
  eps <- .Machine$double.eps^.3
  mod1 <- mod
  yAttr <- attributes(y)
  ytsp <- tsp(y)
  y <- as.matrix(y)
  timeNames <- dimnames(y)[[1]]
  stateNames <- names(mod$m0)
  
  ## define flags for time-varying components
  if (is.null(mod$JFF)) {
    tvFF <- FALSE
  } else {
    tvFF <- TRUE
    nz <- mod$JFF != 0
    mod$JFF <- cbind(row(mod$JFF)[nz], col(mod$JFF)[nz], mod$JFF[nz])
  }
  if (is.null(mod$JV)) {
    tvV <- FALSE
  } else {
    tvV <- TRUE
    nz <- mod$JV != 0
    mod$JV <- cbind(row(mod$JV)[nz], col(mod$JV)[nz], mod$JV[nz])
  }
  if (is.null(mod$JGG)) {
    tvGG <- FALSE
  } else {
    tvGG <- TRUE
    nz <- mod$JGG != 0
    mod$JGG <- cbind(row(mod$JGG)[nz], col(mod$JGG)[nz], mod$JGG[nz])
  }
  if (is.null(mod$JW)) {
    tvW <- FALSE
  } else {
    tvW <- TRUE
    nz <- mod$JW != 0
    mod$JW <- cbind(row(mod$JW)[nz], col(mod$JW)[nz], mod$JW[nz])
  }
  tvFV <- tvFF || tvV
  m <- rbind(mod$m0,matrix(0,nrow=nrow(y),ncol=length(mod$m0))) # filtered values
  a <- matrix(0,nrow=nrow(y),ncol=length(mod$m0))
  f <- matrix(0,nrow=nrow(y),ncol=ncol(y))
  U.C <- vector(1+nrow(y),mode="list")
  D.C <- matrix(0,1+nrow(y),length(mod$m0))
  U.R <- vector(nrow(y),mode="list")
  D.R <- matrix(0,nrow(y),length(mod$m0))
  ## preliminary calculations, if possible (non time-varying case)
  if ( !tvV ) {
    tmp <- La.svd(mod$V,nu=0)
    Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
    if (any(Dv < eps)) {
      Dv <- pmax(Dv, eps)
      warning("a numerically singular 'V' has been slightly perturbed to make it nonsingular")
    }
    Dv.inv <- 1/Dv
    sqrtVinv <- Dv.inv * tmp$vt # t() %*% () = V^{-1}
    if ( !tvFF ) {
      tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
    }
  }
  if ( !tvW ) {
    temp = matrix(0,nrow=nrow(mod$W),ncol=ncol(mod$W))
    #ind = order(diag(mod$W),decreasing=T)
    #temp[cbind(1:nrow(temp),ind)] = 1
    #svdW = list(d=diag(mod$W)[ind], vt = temp)
    svdW <- La.svd(mod$W,nu=0)
    sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
  }
  tmp = list(d = diag(mod$C0), vt = diag(nrow(C0)))
  #tmp <- La.svd(mod$C0,nu=0)
  U.C[[1]] <- t(tmp$vt)
  D.C[1,] <- sqrt(tmp$d)
  for (i in seq(length=nrow(y))) {
    ## set time-varying matrices
    if ( tvFF ) {
      mod$FF[mod$JFF[,-3,drop=FALSE]] <- mod$X[i,mod$JFF[,3]]
    }
    if ( tvV ) {
      mod$V[mod$JV[,-3,drop=FALSE]] <- mod$X[i,mod$JV[,3]]
      tmp = list(d = diag(mod$V), vt = diag(nrow(mod$V)))
      #tmp <- La.svd(mod$V,nu=0)
      Uv <- t(tmp$vt); Dv <- sqrt(tmp$d)
      Dv.inv <- 1/Dv; Dv.inv[abs(Dv.inv)==Inf] <- 0
      sqrtVinv <- Dv.inv * tmp$vt
    }
    if ( tvGG ) {
      mod$GG[mod$JGG[,-3,drop=FALSE]] <- mod$X[i,mod$JGG[,3]]
    }
    if ( tvW ) {
      mod$W[mod$JW[,-3,drop=FALSE]] <- mod$X[i,mod$JW[,3]]
      svdW <- La.svd(mod$W,nu=0)
      sqrtW <- sqrt(svdW$d) * svdW$vt # t()%*%() = W
    }
    if ( tvFV ) {
      #tF.Vinv <- t(mod$FF) %*% crossprod(sqrtVinv)
      #tF.Vinv = (1/mod$V[1,1])*t(mod$FF)
      tF.Vinv = t(mod$FF) %*% diag(1/diag(mod$V))
    }
    
    ## prior
    a[i,] <- mod$GG %*% m[i,]
    tmp <- La.svd(rbind( D.C[i,]*t(mod$GG%*%U.C[[i]]), sqrtW), nu=0) #why transpose here?
    U.R[[i]] <- t(tmp$vt)
    D.R[i,] <- tmp$d
    ## one-step forecast
    f[i,] <- mod$FF %*% a[i,]
    ## posterior
    D.Rinv <- 1/D.R[i,]
    D.Rinv[abs(D.Rinv)==Inf] <- 0
    #tmp <- La.svd(rbind(sqrtVinv[1,1]*mod$FF %*% U.R[[i]],
    #                    diag(x=D.Rinv,nrow=length(D.Rinv))), nu=0)
    tmp <- La.svd(rbind(sqrtVinv %*% mod$FF %*% U.R[[i]],
                        diag(x=D.Rinv,nrow=length(D.Rinv))), nu=0)
    U.C[[i+1]] <- U.R[[i]] %*% t(tmp$vt)
    foo <- 1/tmp$d; foo[abs(foo)==Inf] <- 0
    D.C[i+1,] <- foo
    m[i+1,] <- a[i,] + crossprod(D.C[i+1,]*t(U.C[[i+1]])) %*% 
      tF.Vinv %*% as.matrix(y[i,]-f[i,]) #transpose necessary?
  }
  ans <- list(m=m,U.C=U.C,D.C=D.C,a=a,U.R=U.R,D.R=D.R,f=f)
  
  ans$m <- drop(ans$m); ans$a <- drop(ans$a); ans$f <- drop(ans$f)
  attributes(ans$f) <- yAttr
  if (!is.null(ytsp)) {
    tsp(ans$a) <- ytsp
    tsp(ans$m) <- c(ytsp[1] - 1/ytsp[3], ytsp[2:3])
    class(ans$a) <- class(ans$m) <- if (length(mod$m0) > 1) c("mts","ts") else "ts"
  }
  if (!(is.null(timeNames) && is.null(stateNames)))
    if (is.matrix(ans$a))
    {
      dimnames(ans$a) <- list(timeNames, stateNames)
      dimnames(ans$m) <- list(if(is.null(timeNames)) NULL else c("",timeNames),
                              stateNames)
    }
  else
    if (!is.null(timeNames))
    {
      names(ans$a) <- timeNames
      names(ans$m) <- c("", timeNames)
    }
  if (simplify)
    ans <- c(mod=list(mod1), ans)
  else
  {
    attributes(y) <- yAttr
    ans <- c(y=list(y), mod=list(mod1), ans)
  }
  class(ans) <- "dlmFiltered"
  return(ans)
}

#Start preprocessing for A
#record indices for the response vector and design matrix
ind_design = 1:(nt*k*m)
ind_response = (k*m+1):((nt+1)*k*m)

#create indices to order the alphas by basis instead of variable
order_by_basis = unlist(lapply(1:k,function(b){b + (0:(m-1))*k}))

#initialize design matrix X_alpha to compute vectors of indices we'll need later:
#1. ind_X_alpha, a list of indices of X_alpha which correspond to the position of nonzero elements in X_alpha
#2. ind_Y_alpha, a list of indices of Y_alpha telling us which elements of Y_alpha correspond to the nonzero elements in X_alpha
X_alpha = do.call(rbind,lapply(1:nt,function(t){
  Diagonal(m)%x%Matrix(rowsum(diag((t-1)*k*m + order_by_basis),rep(1:k,each=m)),sparse=TRUE)
}))
ind_X_alpha = which(as.matrix(X_alpha)!=0)
ind_Y_alpha = X_alpha[ind_X_alpha]

#initialize beta_A vector as indices 1:km^2
beta_A = 1:(k*m^2)

#initialize A matrix with the desired pattern
A = matrix(1,m,m)%x%diag(k)
ind_A = which(A!=0)
A[ind_A] = beta_A
A = t(A) #now, the A matrix contains the right indices for the elements of beta_A
ind_beta_A = A[ind_A]

# Some stuff we'll need for the innovation variance
BTB = t(B)%*%B
BTB.inv = solve(BTB)
L.BTB.inv = t(chol(BTB.inv))
L.inv.BTB.inv = solve(L.BTB.inv)

#### MCMC ####

for(i in 1:n_iter){
  #Gibbs Sampling for sigma2_t
  for(j in 1:m){
    j_ind = 1:n + (j-1)*n
    for(t in 1:nt){
      #sigma2
      shape_sigma = a_sigma_0 + n/2
      rate_sigma = b_sigma_0 + 0.5*sum((Y_tilde[t,j_ind] - Y_tilde_hat[t,j_ind])^2)
      sigma2[t,j] = 1/rgamma(1,shape_sigma,rate_sigma)
    }
  }
  
  #Update tau2
  for(j in 1:m){
    shape_tau = a_tau_0 + k*nt/2
    ind = 1:k + (j-1)*k
    rate_tau = b_tau_0 + 0.5*sum(sapply(1:nt,function(t){
      crossprod(L.inv.BTB.inv%*%(alpha[t+1,ind]-(A%*%alpha[t,])[ind,]))
    }))
    tau2[j] = 1/rgamma(1,shape_tau,rate_tau)
  }
  
  #Update A using Gibbs
  #Create vector by reordering the alphas by time, variable, basis
  Y_alpha = as.vector(t(alpha))
  #Place the correct elements of Y_alpha into the nonzero slots of X_alpha
  X_alpha[ind_X_alpha] = (Y_alpha[ind_design])[ind_Y_alpha]
  
  #Cholesky decomposition of linear model covariance inverse
  L.inv.Sigma_full = Diagonal(nt)%x%
    (Diagonal(x = 1/sqrt(tau2))%x%
       Matrix(L.inv.BTB.inv,sparse = TRUE))
  
  #Use cholesky to give diagonal variance
  L_inv_Y = L.inv.Sigma_full%*%Y_alpha[ind_response]
  L_inv_X = L.inv.Sigma_full%*%X_alpha
  
  #posterior precision and mean
  precision_beta_A = crossprod(L_inv_X) #this step takes the longest (big matrix multiplication)
  beta_A_hat = solve(precision_beta_A,
                     crossprod(L_inv_X,L_inv_Y))
  
  #one posterior sample
  beta_A = beta_A_hat + backsolve(chol(precision_beta_A),rnorm(k*m^2))
  
  #finally, place the beta_A terms in the correct spot
  A[ind_A] = beta_A[ind_beta_A]
  
  #Gibbs update for alpha via FFBS step
  #Time-varying sigma2 and tau2
  dlm_alpha = dlm(m0 = m0, C0 = C0,
                  FF = diag(m)%x%Phi, GG = A,
                  V = matrix(0,n*m,n*m), #JV will specify that the diagonal is time-varying
                  W = diag(tau2,m,m)%x%BTB.inv, #W is time-invariant
                  JV = diag(rep(1:m,each=n)),    #The index "i" means that at time t, replace the matching value in V with X[t,i]
                  X = sigma2)
  dlmfilt = dlmFilterSAR(Y_tilde, dlm_alpha, simplify=T) #W is SAR
  alpha = dlmBSample(dlmfilt) #backwards sampling
  
  #update current predicted values with new alpha value
  results[[i]] = list(sigma2 = sigma2, A = A, tau2 = tau2, alpha = alpha)
  Y_tilde_hat = alpha[-1,]%*%t(diag(m)%x%Phi)
}

# for(i in 1:n_iter){
# 
#   #Gibbs Sampling for sigma2_t
#   for(j in 1:m){
#     j_ind = 1:n + (j-1)*n
#     for(t in 1:nt){
#       #sigma2
#       shape_sigma = a_sigma_0 + n/2
#       rate_sigma = b_sigma_0 + 0.5*sum((Y_tilde[t,j_ind] - Y_tilde_hat[t,j_ind])^2)
#       sigma2[t,j] = 1/rgamma(1,shape_sigma,rate_sigma)
#     }
#   }
#   
#   #Update tau2 (and thus Sigma_alpha)
#   for(j in 1:m){
#     shape_tau = a_tau_0 + k*nt/2
#     ind = 1:k + (j-1)*k
#     rate_tau = b_tau_0 + 0.5*sum(sapply(1:nt,function(t){
#       #crossprod(L.inv.Sigma_alpha%*%(alpha[t+1,ind]-(A%*%alpha[t,])[ind,]))
#       nu_t = cbind(alpha[t+1,ind]-(A%*%alpha[t,])[ind,])
#       t(nu_t)%*%BTB%*%nu_t
#     }))
#     tau2[j] = 1/rgamma(1,shape_tau,rate_tau)
#   }
#   
#   #Update Innovation Variance
#   Sigma_alpha = tau2*BTB.inv
#   Sigma_alpha.inv = (1/tau2)*BTB
#   L.Sigma_alpha = chol(Sigma_alpha)
#   L.inv.Sigma_alpha = solve(L.Sigma_alpha)
#   
#   #Update A using Gibbs
#   if(m>1){
#     for(b in 1:k){
#       b_ind = b + (0:(m-1))*k
#       alpha_b = alpha[,b_ind] #alpha coef. for basis b with columns for each variable
#       A[b_ind,b_ind] = suppressWarnings(Acoef(VAR(alpha_b, p=1, type='none'))[[1]])
#     }
#   } else{
#   ar1_out = apply(L.inv.Sigma_alpha%*%t(alpha),1,AR1_gibbs) # After cholesky, the true errors are iid N(0,1). Well, at least if the SAR fits well...
#   A = diag(ar1_out[1,])
#   }
#   
#   #Gibbs update for alpha via FFBS step
#   #Time-varying sigma2 and tau2
#   dlm_alpha = dlm(m0 = m0, C0 = C0,
#                   FF = Phi, GG = A,
#                   V = matrix(0,n,n), #JV will specify that the diagonal is time-varying
#                   W = Sigma_alpha,    #W is now time-invariant
#                   JV = diag(1,n),    #The index "1" means that at time t, replace the diagonal of V with X[t,1]
#                   X = sigma2)
#    dlmfilt = dlmFilterSAR(Y_tilde, dlm_alpha, simplify=T) #W is SAR
#    alpha = dlmBSample(dlmfilt) #backwards sampling
# 
#   
#   #update current predicted values with new alpha value
#   results[[i]] = list(sigma2 = sigma2, A = A, tau2 = tau2, alpha = alpha)
#   Y_tilde_hat = alpha[-1,]%*%t(Phi)
# }
