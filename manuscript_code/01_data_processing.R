#### Part 1: full data processing ####
#Prepare full dataset for all 5 P1 variables.
#Standardized anomalies for 1984-1995

#packages
library(ggplot2)
library(dplyr)
library(vars)
library(gridExtra)
library(Rcpp)
library(tictoc)
library(LatticeKrig)
library(RcppClock)
library(patchwork)
library(scales)
library(latex2exp)

#user needs to set directories
#setwd('') #set to location of repo
data_dir = "data/"

#### Helper functions (plotting) ####
plot_spatial = function(yt,v = c('AOD','Long Wave','Strat Temp')){
  plot_data = data.frame(coords)
  names(plot_data) = c("lon","lat")
  plot_data = rbind(plot_data,plot_data,plot_data)
  plot_data$var = c(rep(v[1],n),rep(v[2],n),rep(v[3],n))
  plot_data$y = c(yt)
  yp = ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y,color=y),linewidth=0) + 
    facet_grid(.~var)+
    scale_fill_gradient2(low='blue',mid='white',high='red',limits=c(-3.5,3.5),na.value='gray60')+
    scale_color_gradient2(low='blue',mid='white',high='red',limits=c(-3.5,3.5),guide='none',na.value='gray60')+
    theme_minimal() + 
    ggtitle(paste("MERRA-2 data"))+ 
    theme(aspect.ratio=1/1.6)
  print(yp)
}

plot_spatiotemporal = function(yt,v = c('AOD','Long Wave','Strat Temp'),t_labs=c()){
  if(length(t_labs)==0){
    t_labs=1:nrow(yt)
  }
  c = data.frame(coords)
  names(c) = c("lon","lat")
  dfs = lapply(1:nrow(yt),function(t){
    plot_data = rbind(c,c,c)
    plot_data$var = c(rep(v[1],n),rep(v[2],n),rep(v[3],n))
    plot_data$y = c(yt[t,])
    plot_data$t = t_labs[t]
    plot_data
  })
  plot_data = do.call("rbind", dfs)
  yp = ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y)) + 
    #geom_path(data=map_data('world'),mapping=aes(long,lat,group=group),linewidth = 0.1)+
    facet_grid(t~var)+
    scale_fill_gradient2(low='blue',mid='white',high='red',limits=c(-3.5,3.5))+
    theme_bw() + 
    theme(legend.key.height = unit(1.5, "cm"))+
    labs(title=paste("MERRA-2 Standardized Anomalies"),x=NULL,y=NULL,fill=NULL)
  print(yp)
}

plot_spatiotemporal_quarterly = function(yt,v = c('AOD','LWR','T50'),t_labs=c()){
  #make quarterly summary
  nq = nrow(yt)/3
  yq = matrix(0,nq,ncol(yt))
  for(q in 1:nq){
    yq[q,] = apply(yt[1:3+(q-1)*3,],2,mean)
  }
  yt=yq
  if(length(t_labs)==0){
    t_labs=1:nrow(yt)
  }
  c = data.frame(coords)
  names(c) = c("lon","lat")
  dfs = lapply(1:nrow(yt),function(t){
    plot_data = rbind(c,c,c)
    plot_data$var = c(rep(v[1],n),rep(v[2],n),rep(v[3],n))
    plot_data$y = c(yt[t,])
    plot_data$t = t_labs[t]
    plot_data
  })
  plot_data = do.call("rbind", dfs)
  yp = ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y,color=y),linewidth=0) + 
    geom_path(rworldmap::coastsCoarse,mapping=aes(long,lat,group=group),linewidth = 0.05)+
    facet_grid(var~t,switch='y')+
    scale_fill_gradient2( low='blue',mid='white',high='red',limits=c(-3.15,3.15))+
    scale_color_gradient2(low='blue',mid='white',high='red',limits=c(-3.15,3.15))+
    theme_bw() +
    #coord_cartesian(xlim=c(-180,180))+
    #coord_quickmap(ylim=c(-87,87))+
    xlim(-181,181)+
    #ylim(-90,90)+
    theme(legend.key.height = unit(1, "cm"))+#,aspect.ratio=1/1.6)+
    scale_y_continuous(position = "right")+#,limits = c(-91,91))+
    labs(title=paste("MERRA-2 Standardized Anomalies"),x=NULL,y=NULL,fill=NULL,color=NULL)
  print(yp)
}

plot_spatial_uni = function(yt,v = c('Strat Temp')){
  plot_data = data.frame(coords)
  names(plot_data) = c("lon","lat")
  plot_data$var = c(rep(v,n))
  plot_data$y = c(yt)
  yp = ggplot(data=plot_data)+ 
    geom_tile(aes(x=lon, y=lat, fill=y)) + 
    scale_fill_gradient2(low='blue',mid='white',high='red',limits=c(-3.5,3.5))+
    theme_minimal() + 
    ggtitle(paste(v,"MERRA-2 data"))+ 
    theme(aspect.ratio=1/1.6)
  print(yp)
}


#### Data Processing ####

#data processing function, remove climatologies to compute anomalies
removemonthly <- function(datamat, v){ #input is tibble and v: variable name
  colnames(datamat)[which(colnames(datamat)==v)] <- "v" #renaming variable of interest to neutral variable name
  mavg <- datamat %>% 
    #filter(year<=1990) %>% #optional filtering before calcing climatologies
    group_by(lat, lon, month) %>% 
    summarise(tam = mean(v, na.rm=TRUE),sdm = sd(v,na.rm=TRUE)) %>%
    dplyr::select(lat,lon,month,tam,sdm)
  datamat <- merge(datamat, mavg, all=TRUE)
  datamat <- datamat %>% mutate(tadt = (v-tam)/sdm) %>% dplyr::select(-sdm, -tam)
  
  colnames(datamat)[which(colnames(datamat)=="v")] <- v #renaming column
  return(datamat)
}

start = 1984
end = 1995

#surface temp
temp_surface <- read.csv(paste0(data_dir,"MERRA2_tavg1_2d_slv_monthly_Nx_48x24_FULL.csv")) %>%
  filter(between(year, start, end)) %>%
  rename(ta=T)
temp_surface = removemonthly(temp_surface,"ta") %>%
  arrange(lat,lon,year,month)

#strat temp
temp_strat <- read.csv(paste0(data_dir,"MERRA2_3dasm_temperature_50mb_48x24.csv")) %>%
  filter(between(year, start, end)) %>%
  rename(ta=T) 
temp_strat <- removemonthly(temp_strat, "ta") #removing monthly means at each spatial location - y_tilde
temp_strat <- temp_strat %>% 
  arrange(lat,lon,year,month)

#aod
aermat <- read.csv(paste0(data_dir,"TOTEXTTAUall_48x24.csv")) %>%
  filter(between(year, start, end)) %>%
  rename(aod=TOTEXTTAU)
aermat <- removemonthly(aermat, "aod") %>% 
  rename(aod_dt = tadt)
aermat <- aermat %>% 
  arrange(lat,lon,year,month)

#rad
radmat <- read.csv(paste0(data_dir,"MERRA2_2d_radflux_48x24.csv"))  %>%
  filter(between(year, start, end))
radmat <- removemonthly(radmat, "LWTUP") %>% rename(LWTUP_dt = tadt)
radmat <- removemonthly(radmat, "SWTNT") %>% rename(SWTNT_dt = tadt)
radmat <- radmat %>% 
  arrange(lat,lon,year,month)

# coordinate grid
coords = aermat %>% 
  filter(date==unique(aermat$date)[1]) %>%
  dplyr::select(lon,lat) %>%
  arrange(lat,lon)
n = dim(coords)[1]
nt = length(unique(temp_strat$date))

# observation matrices with no missing data
y_strat = matrix(temp_strat$tadt,nrow=nt)
y_lw = matrix(radmat$LWTUP_dt,nrow=nt)  
y_aod = matrix(aermat$aod_dt,nrow=nt)
y_sw = matrix(radmat$SWTNT_dt,nrow=nt)
y_surface = matrix(temp_surface$tadt,nrow=nt)

y_full_strat = cbind(y_aod,y_lw,y_strat)
y_full_surface = cbind(y_aod,y_sw,y_surface)

rm(aermat,radmat,temp_strat,temp_surface,aerdir,raddir,surfdir,stratdir)
gc()

#Figure 1 in paper
#plot_spatiotemporal_quarterly(y_full_strat[88:99,],t_labs=factor(c('April-June 1991','July-September 1991','October-December 1991','January-March 1992'),levels=c('April-June 1991','July-September 1991','October-December 1991','January-March 1992')))


#### Test set scenario 1 (block) ####

## Scenario 1: missing block for temperature over north america
(miss_lons = unique(coords$lon)[4:20])
(miss_lats = unique(coords$lat)[12:23])

#scenario 1 indices for univariate
s1_uni = which((coords$lon %in% miss_lons) & (coords$lat %in% miss_lats))  

#scenario 1 indices for multivariate
s1_multi = c(s1_uni+2*n)

#scenario 1 indices for bivariate
s1_bivariate = c(s1_uni+n)

#take out three years as a test set
t_s1 = 92:127

#plot to demonstrate
y_s1 = y_full_strat[92,]
y_s1[s1_multi]=NA

#Figure 6 in paper
plot_spatial(y_s1,c('AOD','LWR','T50'))+
  labs(x=NULL,y=NULL,fill=NULL,title='Holdout set for T50 anomalies (August 1991-July 1994)')+
  theme_bw()+
  coord_quickmap()+
  theme(aspect.ratio = 1/1.8)+
  geom_path(rworldmap::coastsCoarse,mapping=aes(long,lat,group=group),linewidth = 0.1)


#### Test set scenario 2 (random) ####

#scenario 1 indices for univariate
set.seed(2324)
s2_uni = sample(1:1152,115)
set.seed(NULL) #"unset" seed, we only want a fixed seed for the test set generation. 

#scenario 1 indices for multivariate
s2_multi = c(s2_uni+2*n)

#scenario 1 indices for bivariate
s2_bivariate = c(s2_uni+n)

#take out three years as a test set
t_s2 = 1:nt

#plot to demonstrate
y_s2 = y_full_strat[92,]
y_s2[s2_multi]=NA

plot_spatial(y_s2,c('AOD','LWR','T50'))+
  labs(x=NULL,y=NULL,fill=NULL,title='Holdout set for T50 anomalies (August 1991-July 1994)')+
  theme_bw()+
  coord_quickmap()+
  theme(aspect.ratio = 1/1.8)+
  geom_path(rworldmap::coastsCoarse,mapping=aes(long,lat,group=group),linewidth = 0.1)

