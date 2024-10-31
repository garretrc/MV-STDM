source('manuscript_code/01_data_processing.R')
library(abind)

#### Load s2 Results ####

load("data/results/results_s2_rw_84_95_1500.RData")
load("data/results/results_s2_rw_84_95_1500_2.RData")
load("data/results/results_s2_uni_84_95_1500.RData")
load("data/results/results_s2_uni_84_95_1500_2.RData")
load("data/results/results_s2_multi_84_95_1500.RData")
load("data/results/results_s2_multi_84_95_1500_2.RData")

#combine chains
burn = 1:500
results_multi$sigma = abind(results_multi$sigma[,,-burn],results_multi_2$sigma[,,-burn])
results_multi$tau = cbind(results_multi$tau[,-burn],results_multi_2$tau[,-burn])
results_multi$A = cbind(results_multi$A[,-burn],results_multi_2$A[,-burn])
results_multi$alpha = abind(results_multi$alpha[,,-burn],results_multi_2$alpha[,,-burn])

results_uni$sigma = cbind(results_uni$sigma[,-burn],results_uni_2$sigma[,-burn])
results_uni$tau = c(results_uni$tau[-burn],results_uni_2$tau[-burn])
results_uni$A = cbind(results_uni$A[,-burn],results_uni_2$A[,-burn])
results_uni$alpha = abind(results_uni$alpha[,,-burn],results_uni_2$alpha[,,-burn])

results_rw$sigma = cbind(results_rw$sigma[,-burn],results_rw_2$sigma[,-burn])
results_rw$tau = c(results_rw$tau[-burn],results_rw_2$tau[-burn])
results_rw$A = cbind(results_rw$A[,-burn],results_rw_2$A[,-burn])
results_rw$alpha = abind(results_rw$alpha[,,-burn],results_rw_2$alpha[,,-burn])

rm(results_multi_2,results_uni_2,results_rw_2)
gc()


#### Predictive Assessment ####

#Quick posterior summary for tau2
mean(results_rw$tau)
quantile(results_rw$tau,c(0.025,0.975))

mean(results_uni$tau)
quantile(results_uni$tau,c(0.025,0.975))

mean(results_multi$tau[3,])
quantile(results_multi$tau[3,],c(0.025,0.975))

#calculate y_preds
Phi = results_multi$Phi
k = dim(Phi)[2]
ypred_rw = apply(results_rw$alpha[-1,,],1:2,mean)%*%t(Phi)
ypred_uni = apply(results_uni$alpha[-1,,],1:2,mean)%*%t(Phi)
ypred_multi = apply(results_multi$alpha[-1,k*2+1:k,],1:2,mean)%*%t(Phi)

#calculate residuals
eps_rw = y_strat - ypred_rw
eps_uni = y_strat - ypred_uni
eps_multi = y_strat - ypred_multi

#train set rmse
rmse = function(x){sqrt(mean(x^2))}
rmse(eps_rw[-t_s2,-s2_uni])
rmse(eps_uni[-t_s2,-s2_uni])
rmse(eps_multi[-t_s2,-s2_uni])

#test set rmse
rmse(eps_rw[t_s2,s2_uni])
rmse(eps_uni[t_s2,s2_uni])
rmse(eps_multi[t_s2,s2_uni])

test_rmse = data.frame(
  RMSE = c(apply(eps_rw[t_s2,s2_uni],1,rmse),
           apply(eps_uni[t_s2,s2_uni],1,rmse),
           apply(eps_multi[t_s2,s2_uni],1,rmse)),
  fit = rep(c('Univariate-RW','Univariate','Multivariate'),each=length(t_s2)),
  time = rep(seq(as.Date("1991/8/1"), by = "month", length.out = length(t_s2)),3)
)

ggplot(test_rmse)+ 
  geom_line(aes(x=time,y=RMSE,color=fit))


# CRPS univariate-RW
alpha_samps = results_rw$alpha[-1,,]
y_samps = apply(alpha_samps,3,function(alpha){alpha%*%t(Phi)},simplify=F)
y_samps = abind::abind(y_samps,along=3)
y_test_samps = y_samps[t_s2,s2_uni,]
s2_samps = results_rw$sigma[t_s2,]
me_samps = aperm(a = apply(s2_samps,1:2,function(s2){rnorm(length(s2_uni),0,sqrt(s2))}), perm = c(2,1,3))
ypred_samps = y_test_samps + me_samps
d = dim(ypred_samps)[1:2]
crps_rw = matrix(NA,d[1],d[2])
ypred_test = ypred_rw[t_s2,s2_uni]
for(t in 1:d[1]){
  for(i in 1:d[2]){
    crps_rw[t,i] = scoringRules::crps_sample(ypred_test[t,i],ypred_samps[t,i,])
  }
}

# CRPS univariate
alpha_samps = results_uni$alpha[-1,,]
y_samps = apply(alpha_samps,3,function(alpha){alpha%*%t(Phi)},simplify=F)
y_samps = abind::abind(y_samps,along=3)
y_test_samps = y_samps[t_s2,s2_uni,]
s2_samps = results_uni$sigma[t_s2,]
me_samps = aperm(a = apply(s2_samps,1:2,function(s2){rnorm(length(s2_uni),0,sqrt(s2))}), perm = c(2,1,3))
ypred_samps = y_test_samps + me_samps
d = dim(ypred_samps)[1:2]
crps_uni = matrix(NA,d[1],d[2])
ypred_test = ypred_uni[t_s2,s2_uni]
for(t in 1:d[1]){
  for(i in 1:d[2]){
    crps_uni[t,i] = scoringRules::crps_sample(ypred_test[t,i],ypred_samps[t,i,])
  }
}

# CRPS multivariate
alpha_samps = results_multi$alpha[-1,2*k+1:k,]
y_samps = apply(alpha_samps,3,function(alpha){alpha%*%t(Phi)},simplify=F)
y_samps = abind::abind(y_samps,along=3)
y_test_samps = y_samps[t_s2,s2_uni,]
s2_samps = results_multi$sigma[t_s2,3,]
me_samps = aperm(a = apply(s2_samps,1:2,function(s2){rnorm(length(s2_uni),0,sqrt(s2))}), perm = c(2,1,3))
ypred_samps = y_test_samps + me_samps
d = dim(ypred_samps)[1:2]
crps_multi = matrix(NA,d[1],d[2])
ypred_test = ypred_multi[t_s2,s2_uni]
for(t in 1:d[1]){
  for(i in 1:d[2]){
    crps_multi[t,i] = scoringRules::crps_sample(ypred_test[t,i],ypred_samps[t,i,])
  }
}

#combine and plot
test_crps = data.frame(
  CRPS = c(apply(crps_rw,1,mean),apply(crps_uni,1,mean),apply(crps_multi,1,mean)),
  fit = rep(c('Univariate-RW','Univariate','Multivariate'),each=length(t_s2)),
  time = rep(seq(as.Date("1991/8/1"), by = "month", length.out = length(t_s2)),3)
)

ggplot(test_crps)+ 
  geom_line(aes(x=time,y=CRPS,color=fit))


#Joint plot for both metrics
test_metrics = data.frame(
  metric = c(rep('CRPS',nrow(test_crps)), rep('RMSPE',nrow(test_rmse))),
  fit = c(test_crps$fit,test_rmse$fit),
  time = c(test_crps$time,test_rmse$time),
  value = c(test_crps$CRPS,test_rmse$RMSE)
)

ggplot(test_metrics)+
  geom_line(aes(time,value,color=fit),linewidth=0.9)+
  facet_grid(metric~.,scales='free_y')+
  theme_bw()+
  labs(x='Date',y=NULL,color=NULL,title='Holdout Set Metrics for 50mb Stratospheric Temperature Anomalies')+
  scale_color_manual(values=c('Multivariate'='blue','Univariate'='red','Univariate-RW'='#50C878'))+
  theme(legend.position = 'top')

mean(crps_rw)
mean(crps_uni)
mean(crps_multi)


#### Extra: plots of the predicted values
plot_spatial_uni = function(yt,v = c('Strat Temp')){
  plot_data = data.frame(coords)
  names(plot_data) = c("lon","lat")
  plot_data$var = c(rep(v,n))
  plot_data$y = c(yt)
  yp = ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y)) + 
    scale_fill_gradient2(low='blue',mid='white',high='red',limits=c(-3.5,3.5))+
    theme_minimal() + 
    ggtitle(paste(v))+ 
    theme(aspect.ratio=1/1.6)
  yp
}

t=112
(plot_spatial_uni(y_strat[t,],'True Strat Temp')/
  plot_spatial_uni(ypred_uni[t,],'Univariate Predictions')/
  plot_spatial_uni(ypred_multi[t,],'Multivariate Predictions'))|
(plot_spatial_uni(y_strat[t+1,],'True Strat Temp')/
  plot_spatial_uni(ypred_uni[t+1,],'Univariate Predictions')/
  plot_spatial_uni(ypred_multi[t+1,],'Multivariate Predictions'))|
(plot_spatial_uni(y_strat[t+2,],'True Strat Temp')/
  plot_spatial_uni(ypred_uni[t+2,],'Univariate Predictions')/
  plot_spatial_uni(ypred_multi[t+2,],'Multivariate Predictions'))|
(plot_spatial_uni(y_strat[t+3,],'True Strat Temp')/
  plot_spatial_uni(ypred_uni[t+3,],'Univariate Predictions')/
  plot_spatial_uni(ypred_multi[t+3,],'Multivariate Predictions'))

