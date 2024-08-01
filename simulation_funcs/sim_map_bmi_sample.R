## functions to map bmi to zBMI using the sampling method

#direct to correct code
map.bmi.sample <- function(df_row, bmi_dist = "lognormal", age_dist = "normal", 
                           chart_path, bmi_age_corr = FALSE, num_samps = 10000){
  if(bmi_age_corr==FALSE){
    est <- map.bmi.sample.nocorr(df_row, bmi_dist, age_dist, chart_path, num_samps)
  }else{
    if(age_dist=="uniform"){
      print("Can't sample with correlation when age_dist is uniform")
    }else{
      est <- map.bmi.sample.corr(df_row, bmi_dist, age_dist, chart_path, num_samps)
    }
  }
  est
}

#No correlation-----------------------------------------------------------------
map.bmi.sample.nocorr <- function(df_row, bmi_dist, age_dist, chart_path, num_samps){
  
  #mean and SD age
  mAge <- df_row$Mean_Age
  sdAge <- df_row$SD_Age
  prop.boys <- df_row$Prop_Male #proportion of boys
  
  #chart (must be saved in wd and must be named as follows)
  chart <- df_row$chart
  chart.boys <- read_excel(paste(chart_path, chart,"-LMS-Boys.xlsx", sep = ""))
  chart.girls <- read_excel(paste(chart_path, chart,"-LMS-Girls.xlsx", sep = ""))
  
  #mean and SD B
  b <- df_row$MeanB 
  sb <- df_row$SDB 
  
  #sample bmi
  if(bmi_dist=="lognormal"){
    bs <- b.sample.lnorm(b, sb, num_samps)
  }else if(bmi_dist=="normal"){
    bs <- b.sample.norm(b, sb, num_samps)
  }else{
    print("Not a valid distribution for BMI")
  }
  
  #sample age
  if(age_dist=="normal"){
    ages <- age.sample.norm(mAge, sdAge, chart, num_samps)
  }else if(age_dist=="uniform"){
    ages <- age.sample.unif(mAge, sdAge, chart, num_samps)
  }else{
    print("Not a valid distribution for age")
  }
  
  #sample sex
  sexs <- rbinom(num_samps, 1, prop.boys)#1=boy, 0=girl
  
  #convert to BMI-z
  zs <- z.from.bmi(bs, ages, sexs, chart.boys, chart.girls)
  
  #re-define means/SDs in df_row
  mz.est <- mean(zs)
  sz.est <- sd(zs)
  
  z.est <- c(mz.est, sz.est)
  z.est
}


z.from.bmi <- function(bs, ages, sexs, chart.boys, chart.girls){
  zs <- c()
  for(i in 1:length(bs)){
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
    z <- ((bs[i]/M)^L-1)/(L*S)
    
    zs <- append(zs, z)
  }
  zs
}

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

b.sample.lnorm <- function(meanb, SDb, num_samps){
  
  #calculate mean and variance for log normal from mean/SD of normal
  mu <- log(meanb^2/sqrt(meanb^2 + SDb^2)) #mean for lognormal
  var <- log(1+(SDb^2/meanb^2)) #variance for log normal
  sig <- sqrt(var) #sig for log normal
  
  bs <- rlnorm(n = num_samps, meanlog = mu, sdlog = sig)
  bs
}

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



#With correlation---------------------------------------------------------------
map.bmi.sample.corr <- function(df_row, bmi_dist, age_dist, chart_path, num_samps){
  
  #later (check nocorr is working)
  print("Correlation code not ready")
  
  df_row
}