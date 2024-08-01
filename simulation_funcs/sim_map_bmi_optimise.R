## functions to map bmi to zBMI using the optimization method

map.bmi.optimise <- function(df_row, age_dist, chart_path, num_samps, muz_init=0, sigz_init=1, adj=0.99,
                             max_it=5000, bmi_thresh, z_step){
  #chart_path = directory where charts are saved (must end in /)
  #num_samps is the number of times we sample from the z distribution
  #max_it is the maximum number of iterations allowed in the optimisation
  #muz_init and sigz_init are initial guesses for mean/sd of z distribution
  #adj is the factor by which we reduce z if it exceeds its maximum
  #bmi_thresh is the threshold at which we stop the optimisation
  #z_step is the step size for adjusting the mean/sd of z distribution in optimisation
  
  #mean and sd age (in  months)
  mAge <- df_row$Mean_Age
  sdAge <- df_row$SD_Age
  #proportion of boys
  prop.boys <- df_row$Prop_Male  
  
  #chart (must be saved in wd and must be named as follows)
  chart <- df_row$chart
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  #generate initial sample
  inits <- init.sample(mAge, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, muz_init, sigz_init, adj, num_samps)
  ages <- inits[[1]]
  sexs <- inits[[2]]
  zs <- inits[[3]]
  bs <- inits[[4]]
  max.zs <- inits[[5]]
  
  
  #initialised mean/SD BMI
  meanb_init <- mean(bs)
  sdb_init <- sd(bs)
  
  #Reported BMI values (target)
  meanb_target <- df_row$MeanB
  sdb_target <- df_row$SDB
  
  #optimise
  results <- bmi.optim.routine(meanb_init, sdb_init, meanb_target, sdb_target, 
                               ages, sexs, zs, max.zs, muz_init, sigz_init, adj,
                               chart.boys, chart.girls, max_it, bmi_thresh, z_step)
  
  #order of results: list(muz, sigz, meanz, sdz, meanb, sdb)
  
  #re-define means/SDs in df_row -> zdist, zsamp, bsamp
  est.mzdist <- results[[1]] #muz
  est.szdist <- results[[2]] #sigz
  
  est.mzsamp <- results[[3]] #muz
  est.szsamp <- results[[4]] #sigz
  
  est.mbsamp <- results[[5]] #muz
  est.sbsamp <- results[[6]] #sigz
  
  if(abs(results[[5]]-meanb_target)<= bmi_thresh && abs(results[[6]]-sdb_target)<= bmi_thresh){
    converge <- "Y"
  }else{
    converge <- "N"
  }
  
  zadj <- results[[7]]#n.adj
  
  ests <- list(c(est.mzdist,est.szdist), c(est.mzsamp, est.szsamp), c(est.mbsamp, est.sbsamp), converge, zadj)  
  ests
}



bmi.optim.routine <- function(meanb_init, sdb_init, meanb_target, sdb_target, 
                              ages, sexs, zs, max.zs, muz_init, sigz_init, adj,
                              chart.boys, chart.girls, max_it, bmi_thresh, z_step){
  
  #initialise muz, sigz (parameters of z distribution)
  muz <- muz_init
  sigz <- sigz_init
  
  #initialise meanz, sdz (parameters of z distribution)
  meanz <- mean(zs, na.rm = TRUE)
  sdz <- sd(zs, na.rm = TRUE)
  
  #initialise meanb, sdb
  meanb <- meanb_init
  sdb <- sdb_init
  
  num.adj <- 0
  
  num_it <- 0
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

bmi.from.z <- function(zs, ages, sexs, chart.boys, chart.girls){
  bs <- c()
  for(i in 1:length(zs)){
    #index for LMS values - find value in df closest to age
    #NB: all charts must list age in months as 'AgeMonths', LMS as L, M and S
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

init.sample<-function(mAge, sdAge, age_dist, chart, chart.boys, chart.girls, prop.boys, mu0_z, sig0_z, adj, num_samps){
  #sample age (using functions in map_bmi_sample)
  if(age_dist=="normal"){
    ages <- age.sample.norm(mAge, sdAge, chart, num_samps)
  }else if(age_dist=="uniform"){
    ages <- age.sample.unif(mAge, sdAge, chart, num_samps)
  }else{
    print("Not a valid distribution for age")
  }
  
  #sample sex
  sexs <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
  
  #find max z for each age/sex
  max.zs <- find.max.zs(ages, sexs, chart.boys, chart.girls)
  
  #sample zBMI from normal distribution within limit maxzs
  zs <- z.sample.norm(mu0_z, sig0_z, max.zs, adj, num_samps)
  
  #convert to BMI
  bs <- bmi.from.z(zs, ages, sexs, chart.boys, chart.girls)
  
  inits <- list(ages, sexs, zs, bs, max.zs)
  inits
}


list.to.df <- function(ls, outcome){
  #converts output from optimisation into dataframe for eah outcome
  #outcome = zdist, zsamp, bsamp
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
    newrow <- data.frame(ls[[i]][ind])
    df <- rbind(df, newrow)
  }
  
  df
}