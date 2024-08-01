#Function to generate synthetic data by sampling zBMI from a normal distribution 
#with mean muz and SD sigz
data_gen <- function(muz, sigz, N, chart_path, age_dist = "uniform", 
                     minAge = 5, maxAge=18, chart = "CDC", prop.boys = 0.5, adj=0.99){
  
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  sdAge.yrs <- (maxAge-minAge)/4
  mAge.yrs <- minAge + (2*sdAge.yrs)
  
  sdAge.months <- sdAge.yrs*12
  mAge.months <- mAge.yrs*12
  
  #generate sample of zs, age, sex & associated bmis
  inits <- init.sample(mAge.months, sdAge.months, age_dist, chart, chart.boys, chart.girls, 
                       prop.boys, muz, sigz, adj, N)
  age <- inits[[1]]
  sex <- inits[[2]]
  zBMI <- inits[[3]]
  BMI <- inits[[4]]
  max.zs <- inits[[5]]
  #calculate percentile from zs
  perc <- pnorm(zBMI)
  
  data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  dat <- c(mean(age), sd(age), sum(sex)/length(sex), mean(zBMI, na.rm = TRUE),
           sd(zBMI, na.rm = TRUE), mean(BMI, na.rm=TRUE), sd(BMI, na.rm=TRUE),
           mean(perc, na.rm=TRUE), sd(perc, na.rm=TRUE), chart)
  data <- rbind(data, dat)
  colnames(data) <- c("Mean_Age", "SD_Age", "Prop_Male", "MeanZ", "SDZ", 
                      "MeanB", "SDB", "MeanP", "SDP", "chart")
  
  #convert columns 1 to 9 to numeric
  data[, 1:9] <- lapply(data[, 1:9], as.numeric)
  
  data
  #UNITS:
  #age = months
  #prop male = proportion (decimal)
  #percentile = decimal
}

#Function to generate synthetic data by sampling BMI from a log normal distribution 
#with mean meanb and SD sdB
data_gen_BMI <- function(meanB, sdB, N, chart_path, age_dist = "uniform", 
                         minAge = 5, maxAge=18, chart = "CDC", prop.boys = 0.5, adj=0.99){
  
  BMI <- b.sample.lnorm(meanB, sdB, N)
  
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  sdAge.yrs <- (maxAge-minAge)/4
  mAge.yrs <- minAge + (2*sdAge.yrs)
  
  sdAge.months <- sdAge.yrs*12
  mAge.months <- mAge.yrs*12
  
  age <- age.sample.unif(mAge.months, sdAge.months, chart, N)
  sex <- rbinom(N, 1, prop.boys)#1=boy, 0=girl
  
  #convert to BMI-z
  zBMI <- z.from.bmi(BMI, age, sex, chart.boys, chart.girls)
  
  
  #calculate percentile from zs
  perc <- pnorm(zBMI)
  
  data <- data.frame(matrix(ncol = 10, nrow = 0))
  
  dat <- c(mean(age), sd(age), sum(sex)/length(sex), mean(zBMI, na.rm = TRUE),
           sd(zBMI, na.rm = TRUE), mean(BMI, na.rm=TRUE), sd(BMI, na.rm=TRUE),
           mean(perc, na.rm=TRUE), sd(perc, na.rm=TRUE), chart)
  data <- rbind(data, dat)
  colnames(data) <- c("Mean_Age", "SD_Age", "Prop_Male", "MeanZ", "SDZ", 
                      "MeanB", "SDB", "MeanP", "SDP", "chart")
  
  #convert columns 1 to 9 to numeric
  data[, 1:9] <- lapply(data[, 1:9], as.numeric)
  
  data
  #UNITS:
  #age = months
  #prop male = proportion (decimal)
  #percentile = decimal
}



#Function to estimate zBMI from percentile and BMI aggregate data sets
run_estimation <- function(df_row, chart, chart_path){
  
  ## Percentile to zBMI--------------------------------------
  # Analytical method
  est.p.cdf <- map.p.cdf(df_row)
  df_row$EstMeanZ_p_cdf <- est.p.cdf[1]
  df_row$EstSDZ_p_cdf <- est.p.cdf[2]
  
  # Sampling method
  est.p.samp <- map.p.sample(df_row, distribution = "beta", num_samps = 10000)
  df_row$EstMeanZ_p_samp <- est.p.samp[1]
  df_row$EstSDZ_p_samp <- est.p.samp[2]
  
  # Optimisation method
  est.p.opt <- map.perc.optimise(df_row, num_samps = 1000, muz_init=0, sigz_init=1, max_it=5000, 
                    perc_thresh = 0.005, z_step = 0.002)
  df_row$EstMeanZ_p_opt_dist <- est.p.opt[[1]][1]
  df_row$EstSDZ_p_opt_dist <- est.p.opt[[1]][2]
  
  df_row$EstMeanZ_p_opt_samp <- est.p.opt[[2]][1]
  df_row$EstSDZ_p_opt_samp <- est.p.opt[[2]][2]
  
  df_row$ConvergeP <- est.p.opt[[4]]
  
  # Naive method
  df_row$NaiveP <- qnorm(df_row$MeanP)
  df_row$NaiveP_SD <- qnorm(df_row$SDP) #don't use this
  
  ## BMI to zBMI------------------------------------------------
  # Sampling method
  est.b.samp <- map.bmi.sample(df_row, bmi_dist = "lognormal", age_dist = "uniform", chart_path = chart_path,
                 bmi_age_corr = FALSE, num_samps = 10000)
  df_row$EstMeanZ_b_samp <- est.b.samp[1]
  df_row$EstSDZ_b_samp <- est.b.samp[2]
  
  # Optimisation method
  m_init <- rnorm(1, df_row$MeanZ, sd=0.2)
  s_init <- abs(rnorm(1, df_row$SDZ, sd=0.2))
  est.b.opt <- map.bmi.optimise(df_row, age_dist = "uniform", chart_path = chart_path, num_samps = 1000, 
                   muz_init=m_init, sigz_init=s_init, adj=0.99, max_it=5000, bmi_thresh = 0.1, z_step = 0.01)
  df_row$EstMeanZ_b_opt_dist <- est.b.opt[[1]][1]
  df_row$EstSDZ_b_opt_dist <- est.b.opt[[1]][2]
  
  df_row$EstMeanZ_b_opt_samp <- est.b.opt[[2]][1]
  df_row$EstSDZ_b_opt_samp <- est.b.opt[[2]][2]
  
  df_row$ConvergeB <- est.b.opt[[4]]
  df_row$Num_ZAdjB <- est.b.opt[[5]]
  
  # Naive method
  df_row$NaiveB <- naive_BMI(df_row, df_row$MeanB)
  df_row$NaiveB_SD <- naive_BMI(df_row, df_row$SDB) #don't use this
  
  df_row
}

#function to perform the naive method to map from BMI to zBMI
naive_BMI <- function(df_row, b){
  
  #b, age, sex, chart
  age <- df_row$Mean_Age
  
  if(df_row$Prop_Male >= 0.5){
    sex <- 1
  }else{
    sex <- 0
  }
  
  chart <- df_row$chart
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  #index for LMS values - find value in df closest to age
  #NB: all charts must list age in months as 'AgeMonths', LMS as L, M and S
  index <- which.min(abs(chart.boys$AgeMonths - age)) #same for boys and girls
  
  if(sex==1){#boys
    L <- chart.boys$L[index]
    M <- chart.boys$M[index]
    S <- chart.boys$S[index]
  }else{#girls
    L <- chart.girls$L[index]
    M <- chart.girls$M[index]
    S <- chart.girls$S[index]
  }
  z <- ((b/M)^L-1)/(L*S)
    
  z
}


