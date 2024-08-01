#Example of code used for simulation with SD = 1 for zBMI and SD=3.5 for BMI

#include packages
library(tidyverse)
library(readxl)
library(writexl)
library(parallel)
library(PowerTOST) #for OwensT function
library("nleqslv") #for solving non-linear eqs
#Parallisation libraries: https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745
library(doParallel)
library(foreach)
library(parallel)  
library(ggplot2)

#set working directory
chart_path <- "~\\Charts\\"

#source functions from folder
list.files("simulation_funcs", full.names = TRUE) %>% map(source)

###############
# Simulate from a normal distribution for zBMI (perc to zBMI & scenario 1)
###############

## Generate data----------------

cl <- parallel::makeCluster(detectCores(), outfile="DataGen-ZNorm.txt")
doParallel::registerDoParallel(cl)

df1 <- foreach::foreach(i = 0:100, .combine = rbind, .packages = c("readxl")) %dopar% {
  data_gen(muz = -2 + (i*0.04), sigz = 1, N = 200, chart_path, age_dist = "uniform", 
           minAge = 5, maxAge = 18, chart = "IOTF", prop.boys = 0.5, adj=0.99)}

parallel::stopCluster(cl)

## Run simulation-----------------

cl <- parallel::makeCluster(detectCores(), outfile="RunSim5_15-05-24.txt")
doParallel::registerDoParallel(cl)

results1 <- foreach::foreach(i = 1:nrow(df1), .combine = rbind, .packages = c("PowerTOST", "nleqslv","readxl")) %dopar% {
  run_estimation(df1[i,], chart = "IOTF", chart_path = chart_path)}

parallel::stopCluster(cl)


# Plot results------------------------------
results1$EstSDZ_p_cdf <- abs(results1$EstSDZ_p_cdf)

## percentile to zBMI==========================

## Mean
g<-ggplot(results1, aes(x = MeanZ)) +
  geom_point(aes(y = EstMeanZ_p_cdf, color = "Analytical")) +
  geom_point(aes(y = EstMeanZ_p_samp, color = "Sampling")) +
  geom_point(aes(y = EstMeanZ_p_opt_samp, color = "Optimization")) +
  geom_point(aes(y = NaiveP, color = "Naive")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Analytical" = "forestgreen", "Sampling" = "royalblue", 
                                "Optimization" = "red", "Naive" = "purple")) +
  labs(x = "Simulated mean zBMI", y = "Mean zBMI estimated from percentile",
       title = "Transformation: Percentile to zBMI \nSimulation: Normal distribution for zBMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/PtoZ_FromZ.png", g, width = 8, height = 7, dpi = 300)

## SD
g<-ggplot(results1, aes(x = SDZ)) +
  geom_point(aes(y = EstSDZ_p_cdf, color = "Analytical")) +
  geom_point(aes(y = EstSDZ_p_samp, color = "Sampling")) +
  geom_point(aes(y = EstSDZ_p_opt_samp, color = "Optimization")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Analytical" = "forestgreen", "Sampling" = "royalblue", 
                                "Optimization" = "red")) +
  labs(x = "Simulated SD zBMI", y = "SD zBMI estimated from percentile",
       title = "Transformation: Percentile to zBMI \nSimulation: Normal distribution for zBMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/PtoZ_FromZ_SD.png", g, width = 8, height = 7, dpi = 300)


# BMI to zBMI (scenario 1)=====================

## Mean
g<-ggplot(results1, aes(x = MeanZ)) +
  geom_point(aes(y = EstMeanZ_b_samp, color = "Sampling")) +
  geom_point(aes(y = EstMeanZ_b_opt_samp, color = "Optimization")) +
  geom_point(aes(y = NaiveB, color = "Naive")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red", "Naive" = "purple")) +
  labs(x = "Simulated mean zBMI", y = "Mean zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Normal distribution for zBMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromZ.png", g, width = 8, height = 7, dpi = 300)


## SD
g<-ggplot(results1, aes(x = SDZ)) +
  geom_point(aes(y = EstSDZ_b_samp, color = "Sampling")) +
  geom_point(aes(y = EstSDZ_b_opt_samp, color = "Optimization")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red")) +
  labs(x = "Simulated SD zBMI", y = "SD zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Normal distribution for zBMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromZ_SD.png", g, width = 8, height = 7, dpi = 300)


###############
# Simulate from a log normal distribution for BMI (scenario 2)
###############

## Generate data----------------

cl <- parallel::makeCluster(detectCores(), outfile="DataGen4_15-05-24.txt")
doParallel::registerDoParallel(cl)

df2 <- foreach::foreach(i = 0:100, .combine = rbind, .packages = c("readxl")) %dopar% {
  data_gen_BMI(meanB = 15 + (i*0.13), sdB = 3.5, N = 200, chart_path, age_dist = "uniform", 
               minAge = 5, maxAge=18, chart = "IOTF", prop.boys = 0.5, adj=0.99)}

parallel::stopCluster(cl)

## run simulation-----------------

cl <- parallel::makeCluster(detectCores(), outfile="RunSim4_15-05-24.txt")
doParallel::registerDoParallel(cl)

results2 <- foreach::foreach(i = 1:nrow(df2), .combine = rbind, .packages = c("PowerTOST", "nleqslv","readxl")) %dopar% {
  run_estimation(df2[i,], chart = "IOTF", chart_path = chart_path)}

parallel::stopCluster(cl)

# Plot results------------------------------

# BMI to zBMI (scenario 2)=====================

## Mean
g<-ggplot(results2, aes(x = MeanZ)) +
  geom_point(aes(y = EstMeanZ_b_samp, color = "Sampling")) +
  geom_point(aes(y = EstMeanZ_b_opt_samp, color = "Optimization")) +
  geom_point(aes(y = NaiveB, color = "Naive")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red", "Naive" = "purple")) +
  labs(x = "Simulated mean zBMI", y = "Mean zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Lognormal distribution for BMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromB.png", g, width = 8, height = 7, dpi = 300)

## SD
g<-ggplot(results2, aes(x = SDZ)) +
  geom_point(aes(y = EstSDZ_b_samp, color = "Sampling")) +
  geom_point(aes(y = EstSDZ_b_opt_samp, color = "Optimization")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red")) +
  labs(x = "Simulated SD zBMI", y = "SD zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Lognormal distribution for BMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromB_SD.png", g, width = 8, height = 7, dpi = 300)

## CONVERGED results only

## Mean
g<-ggplot(results2, aes(x = MeanZ)) +
  geom_point(aes(y = EstMeanZ_b_samp, color = "Sampling")) +
  geom_point(data = subset(results2, ConvergeB == "Y"), aes(y = EstMeanZ_b_opt_samp, color = "Optimization")) +
  geom_point(aes(y = NaiveB, color = "Naive")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red", "Naive" = "purple")) +
  labs(x = "Simulated mean zBMI", y = "Mean zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Lognormal distribution for BMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromB_converged.png", g, width = 8, height = 7, dpi = 300)

## SD
g<-ggplot(results2, aes(x = SDZ)) +
  geom_point(aes(y = EstSDZ_b_samp, color = "Sampling")) +
  geom_point(data = subset(results2, ConvergeB == "Y"), aes(y = EstSDZ_b_opt_samp, color = "Optimization")) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") +
  scale_color_manual(values = c("Sampling" = "royalblue", 
                                "Optimization" = "red")) +
  labs(x = "Simulated SD zBMI", y = "SD zBMI estimated from BMI",
       title = "Transformation: BMI to zBMI \nSimulation: Lognormal distribution for BMI \nChart: IOTF") +
  theme(legend.title = element_blank(), legend.position = "right", 
        text = element_text(size = 14), plot.title = element_text(size = 14))
ggsave("/SD_MID/BtoZ_FromB_converged_SD.png", g, width = 8, height = 7, dpi = 300)
