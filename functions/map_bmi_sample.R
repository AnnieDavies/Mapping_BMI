#direct to correct code
map.bmi.sample <- function(df_row, bmi_dist = "lognormal", age_dist = "normal", 
                           chart_path, bmi_age_corr = FALSE, num_samps = 10000){
  if(bmi_age_corr==FALSE){
    df_row <- map.bmi.sample.nocorr(df_row, bmi_dist, age_dist, chart_path, num_samps)
  }else{
    if(age_dist=="uniform"){
      print("Can't sample with correlation when age_dist is uniform")
    }else{
      df_row <- map.bmi.sample.corr(df_row, bmi_dist, age_dist, chart_path, num_samps)
    }
  }
  df_row
}

#Sampling without correlation---------------------------------------------------

#Map from BMI to zBMI using the sampling method
map.bmi.sample.nocorr <- function(df_row, bmi_dist, age_dist, chart_path, num_samps){
  
  #convert age in years to months
  mAge0 <- df_row$Mean_Age*12 #mean age at baseline
  sdAge <- df_row$SD_Age*12 #SD age (same at baseline and FU)
  mAge1 <- mAge0 + df_row$fu_months #mean age at FU = baseline + FU time in months
  prop.boys <- df_row$Prop_Male/100 #proportion of boys (convert from percent to proportion)
  
  #chart with LMS values
  chart <- df_row$chart #name of chart to use
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  ## For each input row we map baseline & FU values for the intervention and reference arm (where available)
  
  #BASELINE=====================================================================
  if(df_row$test_type == "Baseline" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    #Baseline
    #Intervention arm
    bA0 <- df_row$BFmean0A #mean bmi
    sbA0 <- df_row$BFsd0A #SD bmi
    
    #Reference arm
    bB0 <- df_row$BFmean0B
    sbB0 <- df_row$BFsd0B
    
    #sample bmi
    if(bmi_dist=="lognormal"){
      bsA0 <- b.sample.lnorm(bA0, sbA0, num_samps)
      bsB0 <- b.sample.lnorm(bB0, sbB0, num_samps)
    }else if(bmi_dist=="normal"){
      bsA0 <- b.sample.norm(num_samps, bA0, sbA0)
      bsB0 <- b.sample.norm(num_samps, bB0, sbB0)
    }else{
      print("Not a valid distribution for BMI")
    }
    
    #sample age
    if(age_dist=="normal"){
      agesA0 <- age.sample.norm(mAge0, sdAge, chart, num_samps)
      agesB0 <- age.sample.norm(mAge0, sdAge, chart, num_samps)
    }else if(age_dist=="uniform"){
      agesA0 <- age.sample.unif(mAge0, sdAge, chart, num_samps)
      agesB0 <- age.sample.unif(mAge0, sdAge, chart, num_samps)
    }else{
      print("Not a valid distribution for age")
    }
    
    #sample sex
    sexsA0 <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
    sexsB0 <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
    
    #convert from BMI to zBMI
    zsA0 <- z.from.bmi(bsA0, agesA0, sexsA0, chart.boys, chart.girls)
    zsB0 <- z.from.bmi(bsB0, agesB0, sexsB0, chart.boys, chart.girls)
    
    #use results to re-define means/SDs in df_row
    df_row$BFmean0A <- mean(zsA0)
    df_row$BFsd0A <- sd(zsA0)
    
    df_row$BFmean0B <- mean(zsB0)
    df_row$BFsd0B <- sd(zsB0)
    
    
  }
  #FOLLOW-UP=====================================================================
  if(df_row$test_type == "FU" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    #Intervention arm
    bA1 <- df_row$BFmean1A #mean bmi
    sbA1 <- df_row$BFsd1A #SD bmi
    
    #Reference arm
    bB1 <- df_row$BFmean1B
    sbB1 <- df_row$BFsd1B
    
    #sample bmi
    if(bmi_dist=="lognormal"){
      bsA1 <- b.sample.lnorm(bA1, sbA1, num_samps)
      bsB1 <- b.sample.lnorm(bB1, sbB1, num_samps)
    }else if(bmi_dist=="normal"){
      bsA1 <- b.sample.norm(num_samps, bA1, sbA1)
      bsB1 <- b.sample.norm(num_samps, bB1, sbB1)
    }else{
      print("Not a valid distribution for BMI")
    }
    
    #sample age
    if(age_dist=="normal"){
      agesA1 <- age.sample.norm(mAge1, sdAge, chart, num_samps)
      agesB1 <- age.sample.norm(mAge1, sdAge, chart, num_samps)
    }else if(age_dist=="uniform"){
      agesA1 <- age.sample.unif(mAge1, sdAge, chart, num_samps)
      agesB1 <- age.sample.unif(mAge1, sdAge, chart, num_samps)
    }else{
      print("Not a valid distribution for age")
    }
    
    #sample sex
    sexsA1 <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
    sexsB1 <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
    
    #convert from BMI to zBMI
    zsA1 <- z.from.bmi(bsA1, agesA1, sexsA1, chart.boys, chart.girls)
    zsB1 <- z.from.bmi(bsB1, agesB1, sexsB1, chart.boys, chart.girls)
    
    #use results to re-define means/SDs in df_row
    df_row$BFmean1A <- mean(zsA1)
    df_row$BFsd1A <- sd(zsA1)
    
    df_row$BFmean1B <- mean(zsB1)
    df_row$BFsd1B <- sd(zsB1)
  }
  
  df_row$measure <- "BMI-z from BMI"
  df_row
  
}

#sample BMI from a lognormal distribution
b.sample.lnorm <- function(meanb, SDb, num_samps){
  
  #calculate mean and variance for log normal from mean/SD of normal
  mu <- log(meanb^2/sqrt(meanb^2 + SDb^2)) #mean for lognormal
  var <- log(1+(SDb^2/meanb^2)) #variance for log normal
  sig <- sqrt(var) #sig for log normal
  
  bs <- rlnorm(n = num_samps, meanlog = mu, sdlog = sig)
  bs
}

#sample BMI from a normal distribution truncated at 0
b.sample.norm <- function(meanb, SDb, num_samps){
  bs <- c()
  for(i in 1:num_samps){
    b <- -1
    while(b<0){
      b <- rnorm(1, meanb, SDb)
    }
    bs <- append(bs, b)
  }
  bs
}

#sample age from a normal distribution (within range of reference chart)
age.sample.norm <- function(mAge, sdAge, chart, num_samps){
  #sample age range depending on chart
  if(chart=="CDC"){
    low <- 24
    high <- 240.5
  }else if(chart=="IOTF"){
    low <- 24
    high <- 216
  }else if(chart=="WHO"){
    low <- 61
    high <- 228
  }else if(chart=="Flemish"){
    low <- 24
    high <- 252
  }
  
  ages <- c()
  for(i in 1:num_samps){
    age <- 0
    while(age<low || age>high){
      age <- rnorm(1, mAge, sdAge)
    }
    ages <- append(ages, age)
  }
  
  ages
}

#sample age from a uniform distribution (within range of reference chart)
age.sample.unif <- function(mAge, sdAge, chart, num_samps){
  min_age <- mAge - (2*sdAge)
  max_age <- mAge + (2*sdAge)
  
  #sample age range depending on chart
  if(chart=="CDC"){
    low <- 24
    high <- 240.5
  }else if(chart=="IOTF"){
    low <- 24
    high <- 216
  }else if(chart=="WHO"){
    low <- 61
    high <- 228
  }else if(chart=="Flemish"){
    low <- 24
    high <- 252
  }
  
  ages <- c()
  for(i in 1:num_samps){
    age <- 0
    while(age<low || age>high){
      age <- runif(1, min_age, max_age)
    }
    ages <- append(ages, age)
  }
  
  ages
}

#calculate zBMI from BMI
z.from.bmi <- function(bs, ages, sexs, chart.boys, chart.girls){
  #NB: all charts must list age in months as 'AgeMonths', and LMS values as L, M and S
  
  zs <- c()
  for(i in 1:length(bs)){
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
    z <- ((bs[i]/M)^L-1)/(L*S)
    
    zs <- append(zs, z)
  }
  zs
}



#Sampling with correlation------------------------------------------------------
#not yet written
map.bmi.sample.corr <- function(df_row, bmi_dist, age_dist, chart_path, num_samps){
  
  print("Correlation code not ready")
  
  df_row
}