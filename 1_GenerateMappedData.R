#include packages
library(renv)
library(tidyverse)
library(readxl)
library(writexl)
library(parallel)
library(PowerTOST) #for OwensT function
library(doParallel)
library(foreach)
library(parallel)  


data_path <- "Data/"
chart_path <- "Charts/"

#source functions from folder 'functions'
list.files("functions", full.names = TRUE) %>% map(source)

#read in outcome data 
df_1218 <- read_excel(paste(data_path, "12-18_Data.xlsx", sep=""))
df_511 <- read_excel(paste(data_path, "5-11_Data.xlsx", sep="")) 

# Set fixed variables
icc <- 0.02 #ICC 
rho <- 0.9 #correlation coefficient for calculating CS

#combine age groups
df_1218['age']<-"12-18"
df_511['age']<-"5-11"
df <- rbind(df_1218, df_511)

df$studyID <- paste(df$study, df$age) #create studyID including age

#change measure to "BMI-z from proportion" where appropriate
for(i in 1:nrow(df)){
  if(!is.na(df$Map_outcome[i]) && df$Map_outcome[i]=="Proportion"){
    df$measure[i] <- "BMI-z from proportion"
  }
}

##############
# Step 1 - create a dataset to be mapped
##############

df_tomap <- create_ds_maptest(df, icc, rho)

#create unique list by study, arm, time point
df_map <- df_tomap
df_map$uniq <- paste(df_map$studyID, df_map$Aarm, df_map$time)

#subset by outcome
df_z <- subset(df_map, measure=="BMI-z")
df_b <- subset(df_map, measure=="BMI")
df_p <- subset(df_map, measure=="Percentile")


#############
# Step 2 - map arm level values
#############

##PERCENTILE====================================================================
# Perform mapping in parallel
cl <- parallel::makeCluster(detectCores(), outfile="Log_mapping_percentile.txt")
doParallel::registerDoParallel(cl)

## Analytic method--------------------------------------------------------------
p.cdf <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind, .packages = c("PowerTOST", "nleqslv")) %dopar% {
      map.p.cdf(df_p[i,])}


## Sampling method-------------------------------------------------------------- 
p.samp.Beta <- foreach::foreach(i = 1:nrow(df_p), .combine = rbind, .packages = c()) %dopar% {
  map.p.sample(df_p[i,], distribution = "beta", num_samps = 10000)}


## Optimization method----------------------------------------------------------
# Delta-tol = 0.0005 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0001
list.p.opt.min1 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.0005, z_step = 0.0001)}


# Delta-step = 2/5*Delta-tol = 0.0002
list.p.opt.min2 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.0005, z_step = 0.0002)}


# Delta-step = 3/5*Delta-tol = 0.0003
list.p.opt.min3 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.0005, z_step = 0.0003)}


# Delta-tol = 0.001 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0002
list.p.opt.mid1 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.001, z_step = 0.0002)}


# Delta-step = 2/5*Delta-tol = 0.0004
list.p.opt.mid2 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.001, z_step = 0.0004)}


# Delta-step = 1/5*Delta-tol = 0.0003
list.p.opt.mid3 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.001, z_step = 0.0006)}


# Delta-tol = 0.005 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.001
list.p.opt.max1 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.005, z_step = 0.001)}


# Delta-step = 2/5*Delta-tol = 0.002
list.p.opt.max2 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.005, z_step = 0.002)}


# Delta-step = 1/5*Delta-tol = 0.003
list.p.opt.max3 <- foreach::foreach(i = 1:nrow(df_p), .packages = c()) %dopar% {
  map.perc.optimise(df_p[i,], num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.005, z_step = 0.003)}


parallel::stopCluster(cl)#======================================================


#BMI============================================================================
# Perform mapping in parallel
cl <- parallel::makeCluster(detectCores(), outfile="Log_mapping_BMI.txt")
doParallel::registerDoParallel(cl)

## Sampling method--------------------------------------------------------------
# BMI - lognormal, age - normal
b.samp.bLNorm.aNorm <- foreach::foreach(i = 1:nrow(df_b), .combine = rbind, .packages = c("readxl")) %dopar% {
  map.bmi.sample(df_b[i,], bmi_dist = "lognormal", age_dist = "normal", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)}

# BMI - lognormal, age - uniform
b.samp.bLNorm.aUnif <- foreach::foreach(i = 1:nrow(df_b), .combine = rbind, .packages = c("readxl")) %dopar% {
  map.bmi.sample(df_b[i,], bmi_dist = "lognormal", age_dist = "uniform", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)}

# BMI - normal, age - normal
b.samp.bNorm.aNorm <- foreach::foreach(i = 1:nrow(df_b), .combine = rbind, .packages = c("readxl")) %dopar% {
  map.bmi.sample(df_b[i,], bmi_dist = "normal", age_dist = "normal", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)}

# BMI - normal, age - uniform
b.samp.bNorm.aUnif <- foreach::foreach(i = 1:nrow(df_b), .combine = rbind, .packages = c("readxl")) %dopar% {
  map.bmi.sample(df_b[i,], bmi_dist = "normal", age_dist = "uniform", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)}


## Optimization method----------------------------------------------------------
# Delta-tol = 0.01 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.005
list.b.opt.aNorm.minHalf <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.01, z_step = 0.005)}

# Delta-step = 1/5*Delta-tol = 0.002 
list.b.opt.aNorm.min5th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.01, z_step = 0.002)}

# Delta-step = 1/10*Delta-tol = 0.001 
list.b.opt.aNorm.min10th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.01, z_step = 0.001)}


# Delta-tol = 0.05 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.025 
list.b.opt.aNorm.midHalf <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.05, z_step = 0.025)}

# Delta-step = 1/5*Delta-tol = 0.01 
list.b.opt.aNorm.mid5th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.05, z_step = 0.01)}

# Delta-step = 1/10*Delta-tol = 0.005 
list.b.opt.aNorm.mid10th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.05, z_step = 0.005)}


# Delta-tol = 0.1 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.05 
list.b.opt.aNorm.maxHalf <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.1, z_step = 0.05)}

# Delta-step = 1/5*Delta-tol = 0.02
list.b.opt.aNorm.max5th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.1, z_step = 0.02)}

# Delta-step = 1/10*Delta-tol = 0.01 
list.b.opt.aNorm.max10th <- foreach::foreach(i = 1:nrow(df_b), .packages = c("readxl")) %dopar% {
  map.bmi.optimise(df_b[i,], age_dist = "normal", chart_path = chart_path, num_samps = 1000, 
                   muz_init=0, sigz_init=1, zlim=3, max_it=5000, bmi_thresh = 0.1, z_step = 0.01)}

parallel::stopCluster(cl)#======================================================

#remove unnecessary variables
rm(cl, df_1218, df_511, chart_path, data_path, i, icc, rho)

## output results 
save.image(file = "MappedData.RData")
