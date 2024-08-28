#Function to map a row of BMI data onto zBMI scale via optimization method
map.bmi.optimise <- function(df_row, age_dist, chart_path, num_samps, muz_init=0, sigz_init=1, adj=0.99,
                             max_it=5000, bmi_thresh, z_step){
  #chart_path = directory where LMS reference charts are saved (must end in /)
  #num_samps is the number of times we sample from the z distribution
  #max_it is the maximum number of iterations allowed in the optimisation
  #muz_init and sigz_init are initial guesses for mean/sd of z distribution
  #adj is the factor by which we reduce z if it exceeds its maximum
  #bmi_thresh is the tolerance threshold at which we stop the optimisation (Delta-tol)
  #z_step is the step size for adjusting the mean/sd of z distribution in optimisation (Delta-step)
  
  
  #convert age in years to months
  mAge0 <- df_row$Mean_Age*12 #mean age at baseline
  sdAge <- df_row$SD_Age*12 #SD age (same at baseline and FU)
  mAge1 <- mAge0 + df_row$fu_months #mean age at FU = baseline + FU time in months
  prop.boys <- df_row$Prop_Male/100 #proportion of boys (convert from percent to proportion)
  
  #chart with LMS values
  chart <- df_row$chart #name of chart to use
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  
  #3 rows of results will be output: 
  #z_dist (final mu_z, sigma_z) = distribution estimates
  #z_sample (final mean(z), sd(z) from sample) = sample estimates
  #bmi_sample (final mean(b), sd(b) from sample) = corresponding values of BMI sampled during optimization
  #initialise as df_row
  df_row_zdist <- df_row
  df_row_zsamp <- df_row
  df_row_bsamp <- df_row
  
  ## For each input row we map baseline & FU values for the intervention and reference arm (where available)
  
  #BASELINE=====================================================================
  if(df_row$test_type == "Baseline" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    
    #Intervention arm-----------------------------------------------------------
    #generate initial sample
    initsA0 <- init.sample(mAge0, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, muz_init, sigz_init, adj, num_samps)
    agesA0 <- initsA0[[1]]
    sexsA0 <- initsA0[[2]]
    zsA0 <- initsA0[[3]]
    bsA0 <- initsA0[[4]]
    max.zsA0 <- initsA0[[5]]
    
    #initialised mean/SD BMI
    meanb_initA0 <- mean(bsA0)
    sdb_initA0 <- sd(bsA0)
    
    #Reported BMI values (target)
    meanb_targetA0 <- df_row$BFmean0A #mean bmi
    sdb_targetA0 <- df_row$BFsd0A #SD bmi
    
    #optimise
    resultsA0 <- bmi.optim.routine(meanb_initA0, sdb_initA0, meanb_targetA0, sdb_targetA0, 
                                   agesA0, sexsA0, zsA0, max.zsA0, muz_init, sigz_init, adj,
                                   chart.boys, chart.girls, max_it, bmi_thresh, z_step)
    
    #results are output as: list(muz, sigz, meanz, sdz, meanb, sdb)
    
    #use results to re-define means/SDs in zdist, zsamp, and bsamp
    df_row_zdist$BFmean0A <- resultsA0[[1]] #muz
    df_row_zdist$BFsd0A <- resultsA0[[2]] #sigz
    
    df_row_zsamp$BFmean0A <- resultsA0[[3]] #meanz
    df_row_zsamp$BFsd0A <- resultsA0[[4]] #sdz
    
    df_row_bsamp$BFmean0A <- resultsA0[[5]] #meanb
    df_row_bsamp$BFsd0A <- resultsA0[[6]] #sdb
    
    #check if optimization converged:
    if(abs(resultsA0[[5]]-meanb_targetA0)<= bmi_thresh && abs(resultsA0[[6]]-sdb_targetA0)<= bmi_thresh){
      df_row_zdist$convergeA0 <- "Y"
      df_row_zsamp$convergeA0 <- "Y"
      df_row_bsamp$convergeA0 <- "Y"
    }else{
      df_row_zdist$convergeA0 <- "N"
      df_row_zsamp$convergeA0 <- "N"
      df_row_bsamp$convergeA0 <- "N"
    }
    
    #record number of z values that needed to be adjusted
    df_row_zdist$zadjA0 <- resultsA0[[7]]#n.adj
    df_row_zsamp$zadjA0 <- resultsA0[[7]]#n.adj
    df_row_bsamp$zadjA0 <- resultsA0[[7]]#n.adj
    
    #Reference arm--------------------------------------------------------------
    #generate initial sample
    initsB0 <- init.sample(mAge0, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, muz_init, sigz_init, adj, num_samps)
    agesB0 <- initsB0[[1]]
    sexsB0 <- initsB0[[2]]
    zsB0 <- initsB0[[3]]
    bsB0 <- initsB0[[4]]
    max.zsB0 <- initsB0[[5]]
    
    #initialised mean/SD BMI
    meanb_initB0 <- mean(bsB0)
    sdb_initB0 <- sd(bsB0)
    
    #Reported BMI values (target)
    meanb_targetB0 <- df_row$BFmean0B
    sdb_targetB0 <- df_row$BFsd0B
    
    #optimise
    resultsB0 <- bmi.optim.routine(meanb_initB0, sdb_initB0, meanb_targetB0, sdb_targetB0, 
                                   agesB0, sexsB0, zsB0, max.zsB0, muz_init, sigz_init, adj,
                                   chart.boys, chart.girls, max_it, bmi_thresh, z_step)
    
    #results are output as: list(muz, sigz, meanz, sdz, meanb, sdb)
    
    #use results to re-define means/SDs in zdist, zsamp, and bsamp
    df_row_zdist$BFmean0B <- resultsB0[[1]] #muz
    df_row_zdist$BFsd0B <- resultsB0[[2]] #sigz
    
    df_row_zsamp$BFmean0B <- resultsB0[[3]] #meanz
    df_row_zsamp$BFsd0B <- resultsB0[[4]] #sdz
    
    df_row_bsamp$BFmean0B <- resultsB0[[5]] #meanb
    df_row_bsamp$BFsd0B <- resultsB0[[6]] #sdb
    
    #check if optimization converged:
    if(abs(resultsB0[[5]]-meanb_targetB0)<= bmi_thresh && abs(resultsB0[[6]]-sdb_targetB0)<= bmi_thresh){
      df_row_zdist$convergeB0 <- "Y"
      df_row_zsamp$convergeB0 <- "Y"
      df_row_bsamp$convergeB0 <- "Y"
    }else{
      df_row_zdist$convergeB0 <- "N"
      df_row_zsamp$convergeB0 <- "N"
      df_row_bsamp$convergeB0 <- "N"
    }
    
    #record number of z values that needed to be adjusted
    df_row_zdist$zadjB0 <- resultsB0[[7]]#n.adj
    df_row_zsamp$zadjB0 <- resultsB0[[7]]#n.adj
    df_row_bsamp$zadjB0 <- resultsB0[[7]]#n.adj
    
  } #end mapping baseline values
  
  #FOLLOW-UP====================================================================
  if(df_row$test_type == "FU" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    #Intervention arm-----------------------------------------------------------
    #generate initial sample
    initsA1 <- init.sample(mAge1, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, muz_init, sigz_init, adj, num_samps)
    agesA1 <- initsA1[[1]]
    sexsA1 <- initsA1[[2]]
    zsA1 <- initsA1[[3]]
    bsA1 <- initsA1[[4]]
    max.zsA1 <- initsA1[[5]]
    
    #initialised mean/SD BMI
    meanb_initA1 <- mean(bsA1)
    sdb_initA1 <- sd(bsA1)
    
    #Reported BMI values (target)
    meanb_targetA1 <- df_row$BFmean1A #mean bmi
    sdb_targetA1 <- df_row$BFsd1A #SD bmi
    
    #optimise
    resultsA1 <- bmi.optim.routine(meanb_initA1, sdb_initA1, meanb_targetA1, sdb_targetA1, 
                                   agesA1, sexsA1, zsA1, max.zsA1, muz_init, sigz_init, adj,
                                   chart.boys, chart.girls, max_it, bmi_thresh, z_step)
    
    #results are output as: list(muz, sigz, meanz, sdz, meanb, sdb)
    
    #use results to re-define means/SDs in zdist, zsamp, and bsamp
    df_row_zdist$BFmean1A <- resultsA1[[1]] #muz
    df_row_zdist$BFsd1A <- resultsA1[[2]] #sigz
    
    df_row_zsamp$BFmean1A <- resultsA1[[3]] #meanz
    df_row_zsamp$BFsd1A <- resultsA1[[4]] #sdz
    
    df_row_bsamp$BFmean1A <- resultsA1[[5]] #meanb
    df_row_bsamp$BFsd1A <- resultsA1[[6]] #sdb
    
    #check if optimization converged:
    if(abs(resultsA1[[5]]-meanb_targetA1)<= bmi_thresh && abs(resultsA1[[6]]-sdb_targetA1)<= bmi_thresh){
      df_row_zdist$convergeA1 <- "Y"
      df_row_zsamp$convergeA1 <- "Y"
      df_row_bsamp$convergeA1 <- "Y"
    }else{
      df_row_zdist$convergeA1 <- "N"
      df_row_zsamp$convergeA1 <- "N"
      df_row_bsamp$convergeA1 <- "N"
    }
    
    #record number of z values that needed to be adjusted
    df_row_zdist$zadjA1 <- resultsA1[[7]]#n.adj
    df_row_zsamp$zadjA1 <- resultsA1[[7]]#n.adj
    df_row_bsamp$zadjA1 <- resultsA1[[7]]#n.adj
    
    #Reference arm---------------------------------------------------------------
    
    #generate initial sample
    initsB1 <- init.sample(mAge1, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, muz_init, sigz_init, adj, num_samps)
    agesB1 <- initsB1[[1]]
    sexsB1 <- initsB1[[2]]
    zsB1 <- initsB1[[3]]
    bsB1 <- initsB1[[4]]
    max.zsB1 <- initsB1[[5]]
    
    #initialised mean/SD BMI
    meanb_initB1 <- mean(bsB1)
    sdb_initB1 <- sd(bsB1)
    
    #Reported BMI values (target)
    meanb_targetB1 <- df_row$BFmean1B #mean bmi
    sdb_targetB1 <- df_row$BFsd1B #SD bmi
    
    #optimise
    resultsB1 <- bmi.optim.routine(meanb_initB1, sdb_initB1, meanb_targetB1, sdb_targetB1, 
                                   agesB1, sexsB1, zsB1, max.zsB1, muz_init, sigz_init, adj,
                                   chart.boys, chart.girls, max_it, bmi_thresh, z_step)
    
    #results are output as: list(muz, sigz, meanz, sdz, meanb, sdb)
    
    #use results to re-define means/SDs in zdist, zsamp, and bsamp
    df_row_zdist$BFmean1B <- resultsB1[[1]] #muz
    df_row_zdist$BFsd1B <- resultsB1[[2]] #sigz
    
    df_row_zsamp$BFmean1B <- resultsB1[[3]] #meanz
    df_row_zsamp$BFsd1B <- resultsB1[[4]] #sdz
    
    df_row_bsamp$BFmean1B <- resultsB1[[5]] #meanb
    df_row_bsamp$BFsd1B <- resultsB1[[6]] #sdb
    
    #check if optimization converged:
    if(abs(resultsB1[[5]]-meanb_targetB1)<= bmi_thresh && abs(resultsB1[[6]]-sdb_targetB1)<= bmi_thresh){
      df_row_zdist$convergeB1 <- "Y"
      df_row_zsamp$convergeB1 <- "Y"
      df_row_bsamp$convergeB1 <- "Y"
    }else{
      df_row_zdist$convergeB1 <- "N"
      df_row_zsamp$convergeB1 <- "N"
      df_row_bsamp$convergeB1 <- "N"
    }  
    
    #record number of z values that needed to be adjusted
    df_row_zdist$zadjB1 <- resultsB1[[7]]#n.adj
    df_row_zsamp$zadjB1 <- resultsB1[[7]]#n.adj
    df_row_bsamp$zadjB1 <- resultsB1[[7]]#n.adj
    
  } #end mapping FU values
  
  df_row_zdist$measure <- "BMI-z from BMI (dist)"
  df_row_zsamp$measure <- "BMI-z from BMI (samp)"
  df_row_bsamp$measure <- "BMI from BMI (check)"
  
  df_rows <- list(df_row_zdist, df_row_zsamp, df_row_bsamp)
  df_rows
}

# function to sample the initial sample of zBMI values
init.sample<-function(mAge, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, mu0_z, sig0_z, adj, num_samps){
  #sample age (using functions in file map_bmi_sample.R)
  if(age_dist=="normal"){
    ages <- age.sample.norm(mAge, sdAge, chart, num_samps)
  }else if(age_dist=="uniform"){
    ages <- age.sample.unif(mAge, sdAge, chart, num_samps)
  }else{
    print("Not a valid distribution for age")
  }
  
  #sample sex
  sexs <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
  
  #find max z for each combination of sampled age/sex
  max.zs <- find.max.zs(ages, sexs, chart.boys, chart.girls)
  
  #sample zBMI from normal distribution within limit maxzs
  zs <- z.sample.norm(mu0_z, sig0_z, max.zs, adj, num_samps)
  
  #convert to BMI
  bs <- bmi.from.z(zs, ages, sexs, chart.boys, chart.girls)
  
  inits <- list(ages, sexs, zs, bs, max.zs)
  inits
}


#function to find the maximum value of z that will allow us to calculate BMI at each sampled age and sex
find.max.zs <- function(ages, sexs, chart.boys, chart.girls){
  
  max.zs <- c()
  for(i in 1:length(sexs)){
    index <- which.min(abs(chart.boys$AgeMonths - ages[i])) #same for boys and girls
    
    if(sexs[i]==1){#boys
      L <- chart.boys$L[index]
      M <- chart.boys$M[index]
      S <- chart.boys$S[index]
    }else{#girls
      L <- chart.girls$L[index]
      M <- chart.girls$M[index]
      S <- chart.girls$S[index]
    }
    
    maxz <- -1/(L*S)
    max.zs <- append(max.zs, maxz)
    
  }
  max.zs
}

#sample zBMI from a normal distribution (adjusting values that are too large)
z.sample.norm <- function(muz, sigz, max.zs, adj, num_samps){
  #adjust is the fraction by which we adjust the Z sample if it is out of range
  zs <- c()
  for(i in 1:num_samps){
    z <- rnorm(1, muz, sigz)
    if(z > max.zs[i]){
      z <- adj*max.zs[i]
    }
    zs <- append(zs, z)
  }
  zs
}

#function to calculate BMI from zBMI (& age & sex)
bmi.from.z <- function(zs, ages, sexs, chart.boys, chart.girls){
  #NB: all charts must list age in months as 'AgeMonths', and LMS values as L, M and S
  
  bs <- c()
  for(i in 1:length(zs)){
    #find the index for LMS values by finding the age in df closest to sampled age
    index <- which.min(abs(chart.boys$AgeMonths - ages[i])) #same for boys and girls
    
    if(sexs[i]==1){#boys
      L <- chart.boys$L[index]
      M <- chart.boys$M[index]
      S <- chart.boys$S[index]
    }else{#girls
      L <- chart.girls$L[index]
      M <- chart.girls$M[index]
      S <- chart.girls$S[index]
    }
    b <- M*(1 + L*S*zs[i])^(1/L)
    bs <- append(bs, b)
  }
  bs
}


#function to perform the optimization routine for BMI
bmi.optim.routine <- function(meanb_init, sdb_init, meanb_target, sdb_target, 
                              ages, sexs, zs, max.zs, muz_init, sigz_init, adj,
                              chart.boys, chart.girls, max_it, bmi_thresh, z_step){
  
  #initialise muz, sigz (parameters of z distribution)
  muz <- muz_init
  sigz <- sigz_init
  
  #initialise meanz, sdz (of the samples)
  meanz <- mean(zs, na.rm = TRUE)
  sdz <- sd(zs, na.rm = TRUE)
  
  #initialise meanb, sdb
  meanb <- meanb_init
  sdb <- sdb_init
  
  num.adj <- 0 #count the number of z samples we need to adjust
  
  num_it <- 0 #track number of iterations
  #while not converged
  while( (abs(meanb-meanb_target) > bmi_thresh || abs(sdb - sdb_target) > bmi_thresh) && num_it < max_it){
    
    #CHOOSE DIRECTION OF CHANGE--------------------------
    #choose delta mu based on mean BMI
    if(abs(meanb-meanb_target) <= bmi_thresh){ #if mean BMI is within threshold - fix muz
      del_mu <- 0.0
    }else if(meanb > meanb_target){ #if mean BMI is too big - reduce mu z
      del_mu <- -z_step #negative
    }else if(meanb < meanb_target){ #if mean BMI is too small - increase mu z
      del_mu <- z_step #positive
    }
    
    #choose delta sigma based on SD BMI
    if(abs(sdb-sdb_target) <= bmi_thresh){ #if sd BMI is within threshold - fix sigz
      del_sig <- 0.0
    }else if(sdb > sdb_target){ #if sd BMI is too big - reduce sig z
      del_sig <- -z_step #negative
    }else if(sdb < sdb_target){ #if mean BMI is too small - increase mu z
      del_sig <- z_step #positive
    }
    
    #UPDATE PARAMETERS-----------------------------------
    #count number of times z is adjusted
    num.adj <- 0
    #update zs
    #zi(t+1) = (1 + delsig/sig(t))*(zi(t)-mu(t)) + (mu(t)+delmu)
    for(i in 1:length(zs)){
      zs[i] <- (1 + del_sig/sigz)*(zs[i] - muz) + (muz + del_mu)
      
      #restrict range of z based on maxz
      if(zs[i] > max.zs[i]){
        zs[i] <- adj*max.zs[i]
        num.adj <- num.adj + 1
      }
    }
    
    #update mu and sigma z
    muz <- muz + del_mu
    sigz <- sigz + del_sig
    
    #calculate mean/sd of sampled zs
    meanz <- mean(zs, na.rm = TRUE)
    sdz <- sd(zs, na.rm = TRUE)
    
    #update bs (using new zs, same age/sex/LMS)
    bs <- bmi.from.z(zs, ages, sexs, chart.boys, chart.girls)
    
    #update mean/sd bmi
    meanb <- mean(bs, na.rm = TRUE)
    sdb <- sd(bs, na.rm = TRUE)
    
    num_it <- num_it + 1
  }
  
  result <- list(muz, sigz, meanz, sdz, meanb, sdb, num.adj)
  result
} 

#function to convert output from optimisation into dataframe for each outcome
list.to.df <- function(ls, outcome){
  #outcome = zdist (distribution estimate), zsamp (sample estimate), bsamp (BMI samples)
  if(outcome == "zdist"){
    ind <- 1
  }else if(outcome == "zsamp"){
    ind <- 2
  }else if(outcome == "bsamp" || outcome == "psamp"){
    ind <- 3
  }
  else{
    print("Error: outcome must be one from: zdist, zsamp, bsamp, psamp")
  }
  df <- data.frame(ls[[1]][ind])
  for(i in 2:length(ls)){
    print(i)
    newrow <- data.frame(ls[[i]][ind])
    df <- rbind(df, newrow)
  }
  
  df
}