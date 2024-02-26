library(tidyverse)

#set working directory
setwd("/Mapping")

#load data (generated in 1_GenerateMappedData.R)
load("MappedData.RData")

#source functions from folder 'functions'
list.files("functions", full.names = TRUE) %>% map(source)


#############
# STEP 3: extract data from optimisation method
#############

#convert the lists to dataframes for each result (distribution vs sample estimates)

##PERCENTILE====================================================================
# Delta-tol = 0.0005 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0001
p.opt.min1.dist <- list.to.df(list.p.opt.min1, "zdist")
p.opt.min1.samp <- list.to.df(list.p.opt.min1, "zsamp")
p.opt.min1.p <- list.to.df(list.p.opt.min1, "psamp")

# Delta-step = 2/5*Delta-tol = 0.0002
p.opt.min2.dist <- list.to.df(list.p.opt.min2, "zdist")
p.opt.min2.samp <- list.to.df(list.p.opt.min2, "zsamp")
p.opt.min2.p <- list.to.df(list.p.opt.min2, "psamp")

# Delta-step = 3/5*Delta-tol = 0.0003
p.opt.min3.dist <- list.to.df(list.p.opt.min3, "zdist")
p.opt.min3.samp <- list.to.df(list.p.opt.min3, "zsamp")
p.opt.min3.p <- list.to.df(list.p.opt.min3, "psamp")

# Delta-tol = 0.001 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0002
p.opt.mid1.dist <- list.to.df(list.p.opt.mid1, "zdist")
p.opt.mid1.samp <- list.to.df(list.p.opt.mid1, "zsamp")
p.opt.mid1.p <- list.to.df(list.p.opt.mid1, "psamp")

# Delta-step = 2/5*Delta-tol = 0.0004
p.opt.mid2.dist <- list.to.df(list.p.opt.mid2, "zdist")
p.opt.mid2.samp <- list.to.df(list.p.opt.mid2, "zsamp")
p.opt.mid2.p <- list.to.df(list.p.opt.mid2, "psamp")

# Delta-step = 1/5*Delta-tol = 0.0003
p.opt.mid3.dist <- list.to.df(list.p.opt.mid3, "zdist")
p.opt.mid3.samp <- list.to.df(list.p.opt.mid3, "zsamp")
p.opt.mid3.p <- list.to.df(list.p.opt.mid3, "psamp")

# Delta-tol = 0.005 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.001
p.opt.max1.dist <- list.to.df(list.p.opt.max1, "zdist")
p.opt.max1.samp <- list.to.df(list.p.opt.max1, "zsamp")
p.opt.max1.p <- list.to.df(list.p.opt.max1, "psamp")

# Delta-step = 2/5*Delta-tol = 0.002
p.opt.max2.dist <- list.to.df(list.p.opt.max2, "zdist")
p.opt.max2.samp <- list.to.df(list.p.opt.max2, "zsamp")
p.opt.max2.p <- list.to.df(list.p.opt.max2, "psamp")

# Delta-step = 1/5*Delta-tol = 0.003
p.opt.max3.dist <- list.to.df(list.p.opt.max3, "zdist")
p.opt.max3.samp <- list.to.df(list.p.opt.max3, "zsamp")
p.opt.max3.p <- list.to.df(list.p.opt.max3, "psamp")


#BMI============================================================================
# Delta-tol = 0.01 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.005
b.opt.aNorm.minHalf.dist <- list.to.df(list.b.opt.aNorm.minHalf, "zdist")
b.opt.aNorm.minHalf.samp <- list.to.df(list.b.opt.aNorm.minHalf, "zsamp")
b.opt.aNorm.minHalf.bmi <- list.to.df(list.b.opt.aNorm.minHalf, "bsamp")

# Delta-step = 1/5*Delta-tol = 0.002
b.opt.aNorm.min5th.dist <- list.to.df(list.b.opt.aNorm.min5th, "zdist")
b.opt.aNorm.min5th.samp <- list.to.df(list.b.opt.aNorm.min5th, "zsamp")
b.opt.aNorm.min5th.bmi <- list.to.df(list.b.opt.aNorm.min5th, "bsamp")

# Delta-step = 1/10*Delta-tol = 0.001 
b.opt.aNorm.min10th.dist <- list.to.df(list.b.opt.aNorm.min10th, "zdist")
b.opt.aNorm.min10th.samp <- list.to.df(list.b.opt.aNorm.min10th, "zsamp")
b.opt.aNorm.min10th.bmi <- list.to.df(list.b.opt.aNorm.min10th, "bsamp")

# Delta-tol = 0.05 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.025 
b.opt.aNorm.midHalf.dist <- list.to.df(list.b.opt.aNorm.midHalf, "zdist")
b.opt.aNorm.midHalf.samp <- list.to.df(list.b.opt.aNorm.midHalf, "zsamp")
b.opt.aNorm.midHalf.bmi <- list.to.df(list.b.opt.aNorm.midHalf, "bsamp")

# Delta-step = 1/5*Delta-tol = 0.01 
b.opt.aNorm.mid5th.dist <- list.to.df(list.b.opt.aNorm.mid5th, "zdist")
b.opt.aNorm.mid5th.samp <- list.to.df(list.b.opt.aNorm.mid5th, "zsamp")
b.opt.aNorm.mid5th.bmi <- list.to.df(list.b.opt.aNorm.mid5th, "bsamp")

# Delta-step = 1/10*Delta-tol = 0.005 
b.opt.aNorm.mid10th.dist <- list.to.df(list.b.opt.aNorm.mid10th, "zdist")
b.opt.aNorm.mid10th.samp <- list.to.df(list.b.opt.aNorm.mid10th, "zsamp")
b.opt.aNorm.mid10th.bmi <- list.to.df(list.b.opt.aNorm.mid10th, "bsamp")

# Delta-tol = 0.1 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.05 
b.opt.aNorm.maxHalf.dist <- list.to.df(list.b.opt.aNorm.maxHalf, "zdist")
b.opt.aNorm.maxHalf.samp <- list.to.df(list.b.opt.aNorm.maxHalf, "zsamp")
b.opt.aNorm.maxHalf.bmi <- list.to.df(list.b.opt.aNorm.maxHalf, "bsamp")

# Delta-step = 1/5*Delta-tol = 0.02
b.opt.aNorm.max5th.dist <- list.to.df(list.b.opt.aNorm.max5th, "zdist")
b.opt.aNorm.max5th.samp <- list.to.df(list.b.opt.aNorm.max5th, "zsamp")
b.opt.aNorm.max5th.bmi <- list.to.df(list.b.opt.aNorm.max5th, "bsamp")

# Delta-step = 1/10*Delta-tol = 0.01 
b.opt.aNorm.max10th.dist <- list.to.df(list.b.opt.aNorm.max10th, "zdist")
b.opt.aNorm.max10th.samp <- list.to.df(list.b.opt.aNorm.max10th, "zsamp")
b.opt.aNorm.max10th.bmi <- list.to.df(list.b.opt.aNorm.max10th, "bsamp")

#-------------------------------------------------------------------------------

###############
# Step 4 - create vectors for plotting
###############

#ORIG---------------------------------------------------------------------------
arm.df_z <- arm_vec(df_z)
arm.df_p <- arm_vec(df_p)
arm.df_b <- arm_vec(df_b)

#Mapped percentile--------------------------------------------------------------
arm.p.cdf <- arm_vec(p.cdf)
arm.p.samp.Beta <- arm_vec(p.samp.Beta)

# Delta-tol = 0.0005 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0001
arm.p.opt.min1.dist <- arm_vec(p.opt.min1.dist)
arm.p.opt.min1.samp <- arm_vec(p.opt.min1.samp)
arm.p.opt.min1.p <- arm_vec(p.opt.min1.p)

# Delta-step = 2/5*Delta-tol = 0.0002
arm.p.opt.min2.dist <- arm_vec(p.opt.min2.dist)
arm.p.opt.min2.samp <- arm_vec(p.opt.min2.samp)
arm.p.opt.min2.p <- arm_vec(p.opt.min2.p)

# Delta-step = 3/5*Delta-tol = 0.0003
arm.p.opt.min3.dist <- arm_vec(p.opt.min3.dist)
arm.p.opt.min3.samp <- arm_vec(p.opt.min3.samp)
arm.p.opt.min3.p <- arm_vec(p.opt.min3.p)

# Delta-tol = 0.001 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.0002
arm.p.opt.mid1.dist <- arm_vec(p.opt.mid1.dist)
arm.p.opt.mid1.samp <- arm_vec(p.opt.mid1.samp)
arm.p.opt.mid1.p <- arm_vec(p.opt.mid1.p)

# Delta-step = 2/5*Delta-tol = 0.0004
arm.p.opt.mid2.dist <- arm_vec(p.opt.mid2.dist)
arm.p.opt.mid2.samp <- arm_vec(p.opt.mid2.samp)
arm.p.opt.mid2.p <- arm_vec(p.opt.mid2.p)

# Delta-step = 1/5*Delta-tol = 0.0003
arm.p.opt.mid3.dist <- arm_vec(p.opt.mid3.dist)
arm.p.opt.mid3.samp <- arm_vec(p.opt.mid3.samp)
arm.p.opt.mid3.p <- arm_vec(p.opt.mid3.p)

# Delta-tol = 0.005 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/5*Delta-tol = 0.001
arm.p.opt.max1.dist <- arm_vec(p.opt.max1.dist)
arm.p.opt.max1.samp <- arm_vec(p.opt.max1.samp)
arm.p.opt.max1.p <- arm_vec(p.opt.max1.p)

# Delta-step = 2/5*Delta-tol = 0.002
arm.p.opt.max2.dist <- arm_vec(p.opt.max2.dist)
arm.p.opt.max2.samp <- arm_vec(p.opt.max2.samp)
arm.p.opt.max2.p <- arm_vec(p.opt.max2.p)

# Delta-step = 1/5*Delta-tol = 0.003
arm.p.opt.max3.dist <- arm_vec(p.opt.max3.dist)
arm.p.opt.max3.samp <- arm_vec(p.opt.max3.samp)
arm.p.opt.max3.p <- arm_vec(p.opt.max3.p)

#Mapped BMI---------------------------------------------------------------------
arm.b.samp.bLNorm.aNorm <- arm_vec(b.samp.bLNorm.aNorm)
arm.b.samp.bLNorm.aUnif <- arm_vec(b.samp.bLNorm.aUnif)
arm.b.samp.bNorm.aNorm <- arm_vec(b.samp.bNorm.aNorm)
arm.b.samp.bNorm.aUnif <- arm_vec(b.samp.bNorm.aUnif)

# Delta-tol = 0.01 (min)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.005
arm.b.opt.aNorm.minHalf.dist <- arm_vec(b.opt.aNorm.minHalf.dist)
arm.b.opt.aNorm.minHalf.samp <- arm_vec(b.opt.aNorm.minHalf.samp)
arm.b.opt.aNorm.minHalf.bmi <- arm_vec(b.opt.aNorm.minHalf.bmi)

# Delta-step = 1/5*Delta-tol = 0.002 
arm.b.opt.aNorm.min5th.dist <- arm_vec(b.opt.aNorm.min5th.dist)
arm.b.opt.aNorm.min5th.samp <- arm_vec(b.opt.aNorm.min5th.samp)
arm.b.opt.aNorm.min5th.bmi <- arm_vec(b.opt.aNorm.min5th.bmi)

# Delta-step = 1/10*Delta-tol = 0.001 
arm.b.opt.aNorm.min10th.dist <- arm_vec(b.opt.aNorm.min10th.dist)
arm.b.opt.aNorm.min10th.samp <- arm_vec(b.opt.aNorm.min10th.samp)
arm.b.opt.aNorm.min10th.bmi <- arm_vec(b.opt.aNorm.min10th.bmi)

# Delta-tol = 0.05 (mid)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.025
arm.b.opt.aNorm.midHalf.dist <- arm_vec(b.opt.aNorm.midHalf.dist)
arm.b.opt.aNorm.midHalf.samp <- arm_vec(b.opt.aNorm.midHalf.samp)
arm.b.opt.aNorm.midHalf.bmi <- arm_vec(b.opt.aNorm.midHalf.bmi)

# Delta-step = 1/5*Delta-tol = 0.01 
arm.b.opt.aNorm.mid5th.dist <- arm_vec(b.opt.aNorm.mid5th.dist)
arm.b.opt.aNorm.mid5th.samp <- arm_vec(b.opt.aNorm.mid5th.samp)
arm.b.opt.aNorm.mid5th.bmi <- arm_vec(b.opt.aNorm.mid5th.bmi)

# Delta-step = 1/10*Delta-tol = 0.005 
arm.b.opt.aNorm.mid10th.dist <- arm_vec(b.opt.aNorm.mid10th.dist)
arm.b.opt.aNorm.mid10th.samp <- arm_vec(b.opt.aNorm.mid10th.samp)
arm.b.opt.aNorm.mid10th.bmi <- arm_vec(b.opt.aNorm.mid10th.bmi)

# Delta-tol = 0.1 (max)~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Delta-step = 1/2*Delta-tol = 0.05
arm.b.opt.aNorm.maxHalf.dist <- arm_vec(b.opt.aNorm.maxHalf.dist)
arm.b.opt.aNorm.maxHalf.samp <- arm_vec(b.opt.aNorm.maxHalf.samp)
arm.b.opt.aNorm.maxHalf.bmi <- arm_vec(b.opt.aNorm.maxHalf.bmi)

# Delta-step = 1/5*Delta-tol = 0.02
arm.b.opt.aNorm.max5th.dist <- arm_vec(b.opt.aNorm.max5th.dist)
arm.b.opt.aNorm.max5th.samp <- arm_vec(b.opt.aNorm.max5th.samp)
arm.b.opt.aNorm.max5th.bmi <- arm_vec(b.opt.aNorm.max5th.bmi)

# Delta-step = 1/10*Delta-tol = 0.01 
arm.b.opt.aNorm.max10th.dist <- arm_vec(b.opt.aNorm.max10th.dist)
arm.b.opt.aNorm.max10th.samp <- arm_vec(b.opt.aNorm.max10th.samp)
arm.b.opt.aNorm.max10th.bmi <- arm_vec(b.opt.aNorm.max10th.bmi)


##########################
# Step 5 - check convergence of optimization
##########################

#Percentile (optimization)------------------------------------------------------

#out of 54
len <- nrow(arm.p.opt.min1.samp)

min1 <- sum(arm.p.opt.min1.samp$converge=="Y")/len #0.17
min2 <- sum(arm.p.opt.min2.dist$converge=="Y")/len#0.67
min3 <- sum(arm.p.opt.min3.samp$converge=="Y")/len#0.87

mid1 <- sum(arm.p.opt.mid1.samp$converge=="Y")/len#0.62
mid2 <- sum(arm.p.opt.mid2.samp$converge=="Y")/len#0.98 
mid3 <- sum(arm.p.opt.mid3.samp$converge=="Y")/len#1

max1 <- sum(arm.p.opt.max1.samp$converge=="Y")/len#1
max2 <- sum(arm.p.opt.max2.samp$converge=="Y")/len#1 ## choose this one to present
max3 <- sum(arm.p.opt.max3.samp$converge=="Y")/len#1 

#BMI (optimization)-------------------------------------------------------------

#out of 205
len <- nrow(arm.b.opt.aNorm.minHalf.samp)

minHalf2 <- sum(arm.b.opt.aNorm.minHalf.samp$converge=="Y")/len #0.56
min5th2 <- sum(arm.b.opt.aNorm.min5th.samp$converge=="Y")/len#0.83
min10th2 <- sum(arm.b.opt.aNorm.min10th.samp$converge=="Y")/len#0.89

midHalf2 <- sum(arm.b.opt.aNorm.midHalf.samp$converge=="Y")/len#0.60
mid5th2 <- sum(arm.b.opt.aNorm.mid5th.samp$converge=="Y")/len#0.83 
mid10th2 <- sum(arm.b.opt.aNorm.mid10th.samp$converge=="Y")/len#0.88

maxHalf2 <- sum(arm.b.opt.aNorm.maxHalf.samp$converge=="Y")/len#0.58
max5th2 <- sum(arm.b.opt.aNorm.max5th.samp$converge=="Y")/len#0.85
max10th2 <- sum(arm.b.opt.aNorm.max10th.samp$converge=="Y")/len#0.90 ## choose this one to present

##########################
# Step 6 - output data to plot
##########################

## Percentile ------------------------------------------------------------------
#keep only z measurements in percentile vectors
arm.df_z.p <- comp_dfs(arm.df_z, arm.p.cdf)

write.csv(arm.df_z.p, "PlottingData/z_for_perc.csv")
write.csv(arm.p.samp.Beta, "PlottingData/perc_sampling.csv")
write.csv(arm.p.opt.max2.dist, "PlottingData/perc_opt_dist.csv")
write.csv(arm.p.opt.max2.samp, "PlottingData/perc_opt_samp.csv")
write.csv(arm.p.cdf, "PlottingData/perc_analytic.csv")



## BMI -------------------------------------------------------------------------
#keep only z measurements in bmi vectors
arm.df_z.b <- comp_dfs(arm.df_z, arm.b.samp.bNorm.aNorm)

write.csv(arm.df_z.b, "PlottingData/z_for_bmi.csv")
write.csv(arm.b.samp.bLNorm.aUnif, "PlottingData/bmi_sampling_aUnif.csv")
write.csv(arm.b.samp.bLNorm.aNorm, "PlottingData/bmi_sampling_aNorm.csv")
write.csv(arm.b.opt.aNorm.max10th.dist, "PlottingData/bmi_opt_dist.csv")
write.csv(arm.b.opt.aNorm.max10th.samp, "PlottingData/bmi_opt_samp.csv")


# save data
save.image(file = "PlottingVecs.RData")

