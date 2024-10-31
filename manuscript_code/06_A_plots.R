source('manuscript_code/01_data_processing.R')
library(abind)

#### Load basis details ####
lk_options <- list(long=c(-180,180),lat=c(-90,90), nlevel=1, startingLevel=3, kappa=sqrt(4), alpha=c(1))
LKinfo = LKrigSetup(x=data.frame(long=lk_options$long,lat=lk_options$lat),
                    nlevel=lk_options$nlevel,                 #Number of resolutions
                    LKGeometry = 'LKSphere',
                    a.wght = 1+lk_options$kappa^2,
                    startingLevel=lk_options$startingLevel,          #40 locations at coarsest res
                    alpha=lk_options$alpha)           #Relative variance between resolutions
lattice_basis = LKrig.basis(x1 = coords,LKinfo = LKinfo)
Phi = as.matrix(lattice_basis)
k = dim(Phi)[2]
res = LKinfo$latticeInfo$startingLevel
centers = data.frame(toSphere(IcosahedronGrid(res)[[res]]))
colnames(centers) = c('lon','lat')

#### Plotting functions ####

A_vec_to_mat = function(A_vec,m=3){
  A_mat = matrix(1,m,m)%x%diag(k)
  A_index = which(A_mat==1)
  A_mat[A_index] = A_vec
  A_mat
}
p=1
s_sqrt <- function(x){sign(x)*(abs(x)^(1/p))}
is_sqrt <- function(x){(abs(x)^p)*sign(x)}
s_sqrt_trans <- function(){trans_new("s_sqrt",s_sqrt,is_sqrt)}
A_Phi_plot = function(
    f = mean,
    title = 'A matrix posterior mean',
    lims=c(-1,1),
    v = c('AOD','LWR','T50'),
    view_diag = TRUE)
{
  #setup variables we need
  k=dim(Phi)[2]
  labels = c()
  var_t = c()
  var_tm1 = c()
  for(vj in v){
    for(vi in v){
      var_t = c(var_t,paste(vi))
      var_tm1 = c(var_tm1,paste(vj,'(t−1)'))
      labels = c(labels,paste(vi,"~ lag 1",vj))
    }
  }
  
  #wrangle data for ggplot
  if(identical(f,'true')){
    A_post_mean = A
  }else if(identical(f,'diff')){
    A_post_mean = A_vec_to_mat(apply(results$A,1,mean))-
      A_vec_to_mat(apply(results2$A,1,mean))
  }else{
    A_post_mean = A_vec_to_mat(apply(results$A,1,f))
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
  A_alpha_df$A = temp#/max(abs(temp))
  
  lim = max(abs(A_alpha_df$A))
  pal = colorRampPalette(colors = c("blue", "white", "red"))
  
  #https://andrewpwheeler.com/2015/07/31/custom-square-root-scale-with-negative-values-in-ggplot2-r/
  p = ggplot(A_alpha_df)+
    geom_tile(aes(lon,lat,fill=A,color=A),linewidth=0)+
    #facet_wrap(~desc)+
    facet_grid(var_t~var_tm1,switch='y')+
    scale_y_continuous(position = "right")+
    theme_bw()+
    scale_fill_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=lims,
      trans = 's_sqrt')+
    scale_color_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=lims,
      trans = 's_sqrt',
      guide='none')+
    labs(x=NULL,y=NULL,title = title,fill=NULL)+
    coord_map(projection = 'rectangular',lat0=0)+
    theme(legend.key.height = unit(1.25, "cm"),aspect.ratio=1/1.5)#+ 
  #scale_x_continuous(sec.axis = sec_axis(~., name = "Variable at current time point", breaks=NULL, labels=NULL)) +
  #scale_y_continuous(sec.axis = sec_axis(~., name = "Variable at previous time point", breaks=NULL, labels=NULL))
  
  p
}
A_Phi_plot_uni = function(
    f = mean,
    title = 'A matrix posterior mean',
    lims=c(-1,1),
    v = 'Strat Temp')
{
  #setup variables we need
  k=dim(Phi)[2]

  #wrangle data for ggplot
  if(identical(f,'true')){
    A_post_mean = A
  }else if(identical(f,'diff')){
    A_post_mean = A_vec_to_mat(apply(results$A,1,mean),1)-
      A_vec_to_mat(apply(results2$A,1,mean),1)
  }else{
    A_post_mean = A_vec_to_mat(apply(results$A,1,f),1)
  }
  A_alpha_df = data.frame(do.call("rbind",replicate(1,coords,simplify=F)))
  A_alpha_df$A = as.numeric(results$Phi%*%cbind(diag(A_post_mean)))/as.numeric(results$Phi%*%rep(1,k))
  
  lim = max(abs(A_alpha_df$A))
  pal = colorRampPalette(colors = c("blue", "white", "red"))
  
  #https://andrewpwheeler.com/2015/07/31/custom-square-root-scale-with-negative-values-in-ggplot2-r/
  p = ggplot(A_alpha_df)+
    geom_tile(aes(lon,lat,fill=A,color=A),linewidth=0)+
    theme_bw()+
    scale_fill_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=lims,
      trans = 's_sqrt')+
    scale_color_gradientn(
      colors=pals::ocean.balance(100),
      #colors = pal(100),
      limits=lims,
      trans = 's_sqrt',
      guide='none')+
    labs(x=NULL,y=NULL,title = title,fill=NULL)+
    #coord_quickmap()+
    theme(legend.key.height = unit(1, "cm"))#+ 
  #scale_x_continuous(sec.axis = sec_axis(~., name = "Variable at current time point", breaks=NULL, labels=NULL)) +
  #scale_y_continuous(sec.axis = sec_axis(~., name = "Variable at previous time point", breaks=NULL, labels=NULL))
  
  p
}
rm_x_axis = scale_x_continuous(breaks=NULL,labels=NULL)
rm_y_axis = scale_y_continuous(breaks=NULL,labels=NULL)


#### Phi %*% A Plots ####

# load full runs
load("data/results/results_strat_full_multi_84_95_1500_c1.RData")
load("data/results/results_strat_full_multi_84_95_1500_c2.RData")
burn = 1:500
results_multi$sigma = abind(results_multi$sigma[,,-burn],results_multi_2$sigma[,,-burn])
results_multi$tau = cbind(results_multi$tau[,-burn],results_multi_2$tau[,-burn])
results_multi$A = cbind(results_multi$A[,-burn],results_multi_2$A[,-burn])
results_multi$alpha = abind(results_multi$alpha[,,-burn],results_multi_2$alpha[,,-burn])
results = results_multi

#Posterior Mean
(A_Phi_plot(mean,TeX('MERRA2 $\\Phi A_{ij}$ Posterior Mean (1984-1995)'))+theme(aspect.ratio = 1/1.5))#+rm_x_axis+rm_y_axis)

#Credible Interval
((A_Phi_plot(function(x){quantile(x,0.025)},TeX('MERRA-2 $\\Phi A_{ij}$ Credible Interval (Low)'))+rm_y_axis)|
    (A_Phi_plot(function(x){quantile(x,0.975)},TeX('MERRA-2 $\\Phi A_{ij}$ Credible Interval (High)'))))+
  plot_layout(guides='collect')

