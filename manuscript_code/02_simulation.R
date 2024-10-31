#### Simulated DLM Data ####
#setwd('') #set to location of repo on your computer

#Source MERRA-2 Data Processing script to load packages and obtain coordinates/missing data scenarios
source('manuscrip_code/01_data_processing.R')

#### Simulation Parameters ####

## data dimension
n = nrow(coords)  #locations
nt = 144          #time points
m = 3             #variables

## basis
lk_options <- list(long=c(-180,180),lat=c(-90,90), nlevel=1, startingLevel=2, kappa=sqrt(4), alpha=c(1))
LKinfo = LKrigSetup(x=data.frame(long=lk_options$long,lat=lk_options$lat),
                    nlevel=lk_options$nlevel,                 #Number of resolutions
                    LKGeometry = 'LKSphere',
                    a.wght = 1+lk_options$kappa^2,
                    startingLevel=lk_options$startingLevel,          #40 locations at coarsest res
                    alpha=lk_options$alpha)           #Relative variance between resolutions
lattice_basis = LKrig.basis(x1 = coords,LKinfo = LKinfo)
Phi = as.matrix(lattice_basis)
B = LKrig.precision(LKinfo,return.B = T)
k = dim(Phi)[2]

res = LKinfo$latticeInfo$startingLevel
centers = data.frame(toSphere(IcosahedronGrid(res)[[res]]))
colnames(centers) = c('lon','lat')

## sigma2
sigma2 = matrix(2,nt,m)

## tau2
tau2 = rep(5,m)
Sigma_alpha = diag(tau2)%x%solve(t(B)%*%B)

## A
A = matrix(0,k*m,k*m)
A_nonspatial = matrix(c(
   0.8,  0,    0,
  -0.2,  0.6,  0.3,
   0.4, -0.2,  0.6
),m,m,byrow=T)

for(i in 1:m){
  for(j in 1:m){
    #default trend is spatially constant
    A_spatial = rep(1,k)
    
    #For some cross variable dependencies, a linear trend starting at 0 in the south pole
    # if( ((j==1)&(i==2)) | ((j==2)&(i==3)) ){ 
    #   A_spatial = (centers$lat + 90)/180 #
    # }
    
    #For other cross-variable dependencies, a quadratic trend strongest near the equator
    if((j==3)&(i==2)){
      A_spatial = (centers$lat/90)
    }
    
    #For other cross-variable dependencies, a quadratic trend strongest near the equator
    if((j==1)&(i==3)){
      A_spatial = 1-abs(centers$lat/90)^0.5
    }
    
    #Assign A matrix based on nonspatial coefficient + spatial trend
    A[1:k + k*(i-1), 1:k + k*(j-1)] = diag(A_nonspatial[i,j]*A_spatial)
  }
}

## alpha
burn = 100 #time points to add to the initial run
m0 = 0
C0 = 1
alpha_sim = matrix(0,nt+1+burn,m*k)

set.seed(1000)
alpha_sim[1,] = rnorm(m*k,m0,C0)
for(t in 1:(nrow(alpha_sim)-1)){
  alpha_sim[t+1,] = A%*%cbind(alpha_sim[t,]) + mvrnorm(1,rep(0,m*k),Sigma_alpha)
}
image.plot(alpha_sim)
alpha_sim = alpha_sim[-(1:burn),]

#### Simulated Dataset ####

y_sim = matrix(0,nt,m*n)

for(t in 1:nt){
  y_sim[t,] = (diag(m)%x%Phi)%*%cbind(alpha_sim[t+1,]) + rnorm(n*m,rep(0,n*m),rep(sqrt(sigma2[t,]),each=n))   #alpha_sim[t+1,] to account for alpha[0,]
}

spatial_plot = function(t,t_name='',v = c('Sim AOD','Sim Long Wave','Sim Strat Temp')){
  plot_data = data.frame(coords)
  names(plot_data) = c("lon","lat")
  plot_data = rbind(plot_data,plot_data,plot_data)
  plot_data$var = c(rep(v[1],n),rep(v[2],n),rep(v[3],n))
  plot_data$y = c(y_sim[t,])
  
  ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y)) + 
    facet_grid(.~var)+
    scale_fill_gradient2(low='blue',mid='white',high='red')+
    theme_minimal() + 
    ggtitle(paste("Simulated data",t_name))+ 
    theme(aspect.ratio=1/1.6)
}

t=102
spatial_plot(t-2)/
spatial_plot(t-1)/
spatial_plot(t)


#### MCMC Run on Synthetic Dataset ####

n_iter = 2500
mcmc_param <- list(n_mcmc=n_iter, n_burn=0, sigma_param=c(1,1), tau_param = c(1,1)) #should probably add more here

# #MCMC code, commented out to avoid re-running
 sourceCpp("dlm_code/sampling_eqs.cpp")
 source("dlm_code/mvdlm_mcmc.R")
# 
# tic()
 results_sim <- mvstdm_mcmc(data = y_sim, coords = coords,
                            mcmc_param = mcmc_param, lk_options = lk_options)
# toc()
 save(results_sim, file = 'data/results/results_sim_multi_2500.RData')

#load('../data/results/results_sim_multi_2500.RData')

#### Results ####

# Convert vectorized A matrix into matrix form
A_vec_to_mat = function(A_vec){
  A_mat = matrix(1,m,m)%x%diag(k)
  A_index = which(A_mat==1)
  A_mat[A_index] = A_vec
  A_mat
}

# Plot A*Phi matrix 
p=1 # controls power of the color scale. p=1 is no transformation.
s_sqrt <- function(x){sign(x)*(abs(x)^(1/p))}
is_sqrt <- function(x){(abs(x)^p)*sign(x)}
s_sqrt_trans <- function(){trans_new("s_sqrt",s_sqrt,is_sqrt)}
A_Phi_plot = function(
    f = mean,
    title = 'A matrix posterior mean',
    lims=c(-1,1),
    v = c('= 1','= 2','= 3'),
    view_diag = TRUE)
{
  #setup variables we need
  k=dim(Phi)[2]
  labels = c()
  var_t = c()
  var_tm1 = c()
  for(vj in v){
    for(vi in v){
      var_t = c(var_t,paste('i',vi))
      var_tm1 = c(var_tm1,paste('j',vj))
      labels = c(labels,paste(vi,"~ lag 1",vj))
    }
  }

  #wrangle data for ggplot
  if(identical(f,'true')){
    A_post_mean = A
  }else if(identical(f,'diff')){
    A_post_mean = A_vec_to_mat(apply(results$A[,burn],1,mean))-
      A_vec_to_mat(apply(results2$A[,burn],1,mean))
  }else{
    A_post_mean = A_vec_to_mat(apply(results$A[,burn],1,f))
  }
  A_post_mean = t(sapply(1:k,function(j){A_post_mean[j-1+c(1,k+1,2*k+1),j-1+c(1,k+1,2*k+1)]}))
  A_alpha_df = data.frame(do.call("rbind",replicate(9,coords,simplify=F)))
  A_alpha_df$index = rep(c(1,4,7,2,5,8,3,6,9),each=n)
  A_alpha_df$var_t = rep(var_t,each=n)
  A_alpha_df$var_tm1 = rep(var_tm1,each=n)
  A_alpha_df$desc = rep(labels,each=n)
  A_alpha_df$desc = factor(A_alpha_df$desc,ordered = T)
  temp=c()
  for(i in 1:9){
    #set diagonal to 0 if requested
    if((view_diag==FALSE) & (i %in% c(1,5,9))){
      temp_i = rep(0,n)
    }else{
      temp_i = as.numeric(results$Phi%*%A_post_mean[,i])/as.numeric(results$Phi%*%rep(1,k))
    }
    temp=c(temp,temp_i)
  }
  A_alpha_df$A = temp#/11.56 #/max(abs(temp))

  lim = max(abs(A_alpha_df$A))
  pal = colorRampPalette(colors = c("blue", "white", "red"))

  #https://andrewpwheeler.com/2015/07/31/custom-square-root-scale-with-negative-values-in-ggplot2-r/
  p = ggplot(A_alpha_df)+
    geom_tile(aes(lon,lat,fill=A,color=A),linewidth=0)+
    #facet_wrap(~desc)+
    facet_grid(var_t~var_tm1,switch='y')+
    theme_bw()+
    scale_y_continuous(position = 'right')+
    scale_fill_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=c(-1,1),
      trans = 's_sqrt')+
    scale_color_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=c(-1,1),
      trans = 's_sqrt',
      guide='none')+
    labs(x=NULL,y=NULL,title = title,fill=NULL,color=NULL)+
    #coord_quickmap()+
    theme(legend.key.height = unit(1.5, "cm"))#+

  p
}

#Plot raw A matrix
A_plot = function(
    f = mean,
    title = 'A matrix posterior mean',
    lims=c(-1,1),
    v = c('= 1','= 2','= 3'),
    view_diag = TRUE)
{
  #setup variables we need
  k=dim(Phi)[2]
  labels = c()
  var_t = c()
  var_tm1 = c()
  for(vj in v){
    for(vi in v){
      var_t = c(var_t,paste('i',vi))
      var_tm1 = c(var_tm1,paste('j',vj))
      labels = c(labels,paste(vi,"~ lag 1",vj))
    }
  }
  
  #wrangle data for ggplot
  if(identical(f,'true')){
    A_post_mean = A
  }else if(identical(f,'diff')){
    A_post_mean = A_vec_to_mat(apply(results$A[,burn],1,mean))-
      A_vec_to_mat(apply(results2$A[,burn],1,mean))
  }else{
    A_post_mean = A_vec_to_mat(apply(results$A[,burn],1,f))
  }
  A_post_mean = t(sapply(1:k,function(j){A_post_mean[j-1+c(1,k+1,2*k+1),j-1+c(1,k+1,2*k+1)]}))
  A_alpha_df = data.frame(do.call("rbind",replicate(9,centers,simplify=F)))
  A_alpha_df$index = rep(c(1,4,7,2,5,8,3,6,9),each=k)
  A_alpha_df$var_t = rep(var_t,each=k)
  A_alpha_df$var_tm1 = rep(var_tm1,each=k)
  A_alpha_df$desc = rep(labels,each=k)
  A_alpha_df$desc = factor(A_alpha_df$desc,ordered = T)
  temp=c()
  for(i in 1:9){
    #set diagonal to 0 if requested
    if((view_diag==FALSE) & (i %in% c(1,5,9))){
      temp_i = rep(0,n)
    }else{
      temp_i = A_post_mean[,i]
    }
    temp=c(temp,temp_i)
  }
  A_alpha_df$A = temp
  
  lim = max(abs(A_alpha_df$A))
  pal = colorRampPalette(colors = c("blue", "white", "red"))
  
  #https://andrewpwheeler.com/2015/07/31/custom-square-root-scale-with-negative-values-in-ggplot2-r/
  p = ggplot(A_alpha_df)+
    geom_point(aes(lon,lat,fill=A),color='black',size=2.5,shape=21)+
    #facet_wrap(~desc)+
    facet_grid(var_t~var_tm1,switch='y')+
    theme_bw()+
    scale_y_continuous(position = 'right')+
    scale_fill_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=c(-1,1),
      trans = 's_sqrt')+
    scale_color_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=c(-1,1),
      trans = 's_sqrt',
      guide='none')+
    labs(x=NULL,y=NULL,title = title,fill=NULL,color=NULL)+
    #coord_quickmap()+
    theme(legend.key.height = unit(1.5, "cm"))+
    ylim(-95,95)+
    xlim(-160,195)+
    scale_x_continuous(breaks=c(-100,0,100),limits=c(-155,185))+
    scale_y_continuous(breaks=c(-50,0,50),limits=c(-98,98),position='right')#+
  #scale_x_continuous(sec.axis = sec_axis(~., name = "Variable at current time point", breaks=NULL, labels=NULL)) +
  #scale_y_continuous(sec.axis = sec_axis(~., name = "Variable at previous time point", breaks=NULL, labels=NULL))
  
  p
}

#extra ggplot stuff to remove axis labels where applicable
rm_x_axis = scale_x_continuous(breaks=NULL,labels=NULL)
rm_y_axis = scale_y_continuous(breaks=NULL,labels=NULL)


burn = c(501:2500) 
results=results_sim

#Figure 2 in manuscript
(A_plot('true',TeX('Simulation $A_{ij}$ True Values'))+guides(fill='none'))|
(A_Phi_plot('true',TeX('Simulation $\\Phi A_{ij}$ True Values'))+theme(legend.key.height = unit(1.25, "cm")))

#Figure 3 in manuscript
(A_Phi_plot('true',TeX('Simulation $\\Phi A_{ij}$ True Values'))+ guides(fill="none")+rm_y_axis)+
(A_Phi_plot(mean,TeX('$\\Phi A_{ij}$ Posterior Mean')))+
(A_Phi_plot(function(x){quantile(x,0.025)},TeX('$\\Phi A_{ij}$ Credible Interval (Low)'))+rm_y_axis)+
(A_Phi_plot(function(x){quantile(x,0.975)},TeX('$\\Phi A_{ij}$ Credible Interval (High)')))+
  plot_layout(guides='collect')

# Posterior mean only
(A_Phi_plot(mean,TeX('$\\Phi A_{ij}$ Posterior Mean (Low Res, 1986-1995)'),lims=c(-8.85,8.85))+rm_x_axis+rm_y_axis)


#Tau2 and sigma2 summary (Table 1 in manuscript)

#tau2 summary
digits = 2
#Posterior Means
round(apply(results$tau[,burn],1,mean),digits)
#Posterior q2.5
round(apply(results$tau[,burn],1,quantile,0.025),digits)
#Posterior q97.5
round(apply(results$tau[,burn],1,quantile,0.975),digits)

#Traces
burn=501:2500
(ggplot()+
  geom_line(aes(x=burn,y=results$tau[1,burn]))+
  geom_hline(aes(yintercept=tau2[1])))/
(ggplot()+
  geom_line(aes(x=burn,y=results$tau[2,burn]))+
  geom_hline(aes(yintercept=tau2[2])))/
(ggplot()+
  geom_line(aes(x=burn,y=results$tau[3,burn]))+
  geom_hline(aes(yintercept=tau2[3])))


#sigma2_t summary
t=1
#Posterior Means
round(apply(results$sigma[t,,burn],1,mean),digits)
#Posterior q2.5
round(apply(results$sigma[t,,burn],1,quantile,0.025),digits)
#Posterior q97.5
round(apply(results$sigma[t,,burn],1,quantile,0.975),digits)

#Traces
(ggplot()+
    geom_line(aes(x=burn,y=results$sigma[t,1,burn]))+
    geom_hline(aes(yintercept=sigma2[1])))/
(ggplot()+
    geom_line(aes(x=burn,y=results$sigma[t,2,burn]))+
    geom_hline(aes(yintercept=sigma2[2])))/
(ggplot()+
    geom_line(aes(x=burn,y=results$sigma[t,3,burn]))+
    geom_hline(aes(yintercept=sigma2[3])))

#sigma2 coverage over all time points
ci_low = apply(results$sigma[,,burn],1:2,quantile,probs=0.025)
ci_high = apply(results$sigma[,,burn],1:2,quantile,probs=0.975)
sig_coverage = (results$sigma[,,burn[1]]>ci_low)&(results$sigma[,,burn[1]]<ci_high)
for(i in 2:length(burn)){
  sig_coverage = sig_coverage + (results$sigma[,,burn[i]]>ci_low)&(results$sigma[,,burn[i]]<ci_high)
}

#average over all variables
mean(sig_coverage)

#average per each variable
apply(sig_coverage,2,mean)

