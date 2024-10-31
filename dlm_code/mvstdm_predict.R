# need to extend to multivariate case
# will also need to take exp transform for AOD
mvstdm_predict <- function(mcmc_results, coords, t_idx, loc_idx, y_test){  
    #testdata is matrix of time index, long and lat
    #y_test is true values
    #n_var is number of variables
    #for now coding this for known grid points, observed in the full dataset
    #for truly unobserved locations, we need to code up the wendland basis fucntions, to "project" onto Phi
    
    Phi = mcmc_results$Phi
    alpha_samps = mcmc_results$alpha[-1,,] # remove alpha0
    sigma_samps = mcmc_results$sigma

    n_mcmc = dim(results$sigma)[3]
    n_pred = nrow(y_test)
    n_var = dim(results$sigma)[2]
    n_basis = dim(Phi)[2]
    metrics = matrix(nrow=n_var, ncol = 4)
    
    for(m in n_var){
      
      y_pred = matrix(ncol=n_mcmc, nrow=n_pred)
      for(i in 1:n_mcmc){
          y_pred[,i] = diag(Phi[loc_idx,] %*% t(alpha_samps[t_idx,((m-1)*n_basis+1):((m-1)*n_basis+n_basis),i]))
          sig_vec = sigma_samps[t_idx,m,i] 
          y_pred[,i] = y_pred[,i] + sapply(sig_vec, function(x) rnorm(1,0,sqrt(x)))
      }

      #MSPE
      pred<-apply(y_pred,1,mean)
      MSPE<-mean((pred-y_test[,m])^2)

      #PMCC
      PMCC<-mean((y_test[,m]-apply(y_pred,1,mean))^2/apply(y_pred,1,function(x) var(x)))+sum(apply(y_pred,1,function(x) log(var(x))))

      #CRPS
      meandiff<-function(obs){
        md<-0
        for(k in 1:dim(y_pred)[2]){
          md<-((k-1)*md+mean(abs(y_pred[obs,k]-y_pred[obs,])))/k
        }
        return(md)
      }
      
      crps1<-sapply(1:n_pred,function(x) meandiff(x))
      crps2<-rep(0,n_pred)
      for(k in 1:n_pred){
        crps2[k]<-mean(abs(y_pred[k,]-y_test[k,m]))
      }
      crps_star=-0.5*crps1+crps2
      CRPS<-mean(crps_star)

      #coverage
      pred_intervals<-t(apply(y_pred,1,function(x) quantile(x,probs=c(0.025,0.975))))
      pred_intervals<-cbind(y_test[,m],pred_intervals)
      count=0
      for(k in 1:n_pred){
          if(pred_intervals[k,1]<=pred_intervals[k,3] && pred_intervals[k,1]>=pred_intervals[k,2]){
          count=count+1
          }else{
            count=count
          }
      }
      coverage<-count/n_pred
    
      metrics[m,] <- c(MSPE, PMCC, CRPS, coverage)
    
    }
    colnames(metrics) <-  c("MSPE", "PMCC", "CRPS", "coverage")
}