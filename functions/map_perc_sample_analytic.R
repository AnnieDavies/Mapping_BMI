## ANALYTIC METHOD

#Function to map a row of percentile data onto zBMI scale via analytic method
map.p.cdf <- function(df_row){
  
  ## Baseline
  if(df_row$test_type == "Baseline" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    
    #Intervention arm----------------------
    pA0 <- df_row$BFmean0A/100
    spA0 <- df_row$BFsd0A/100
    
    #solve simultaneous equations via Newton method
    resA0 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pA0, spA0, method = "Newton")
    
    #save results
    df_row$BFmean0A <- resA0$x[1]
    df_row$BFsd0A <- resA0$x[2]
    
    #Reference arm------------------------
    pB0 <- df_row$BFmean0B/100
    spB0 <- df_row$BFsd0B/100
    
    #solve simultaneous equations via Newton method
    resB0 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pB0, spB0, method = "Newton")
    
    #save results
    df_row$BFmean0B <- resB0$x[1]
    df_row$BFsd0B <- resB0$x[2]
  }
  
  ## Follow-up
  if(df_row$test_type == "FU" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    
    #Intervention arm-----------------------
    pA1 <- df_row$BFmean1A/100
    spA1 <- df_row$BFsd1A/100
    
    resA1 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pA1, spA1, method = "Newton")
    
    #save results
    df_row$BFmean1A <- resA1$x[1]
    df_row$BFsd1A <- resA1$x[2]
    
    #Reference arm--------------------------
    pB1 <- df_row$BFmean1B/100
    spB1 <- df_row$BFsd1B/100
    
    resB1 <- nleqslv(c(1,1), sim_eqs, jac = NULL, pB1, spB1, method = "Newton")
    
    #save results 
    df_row$BFmean1B <- resB1$x[1]
    df_row$BFsd1B <- resB1$x[2]
  }
  
  df_row$measure <- "BMI-z from percentile"
  
  df_row
  
}

#simultaneous equations for expectation and variance of percentile
sim_eqs <- function(x, meanp, SDp){
  z = x[1]
  s = x[2]
  eq1 = meanp - pnorm(z/sqrt(1+s^2))
  
  h <- z/sqrt(1 + s^2)
  a <- 1/sqrt(1 + 2*s^2)
  
  eq2 = SDp^2 - (pnorm(z/sqrt(1+s^2))-2*OwensT(h,a)-pnorm(z/sqrt(1+s^2))^2)
  result <- c(eq1, eq2)
  result
}


## SAMPLING METHOD

#Function to map a row of percentile data onto zBMI scale via the sampling method
map.p.sample <- function(df_row, distribution = "beta", num_samps = 10000){
  #atm distribution can be beta or normal
  
  ## BASELINE
  if(df_row$test_type == "Baseline" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){

    #Intervention arm
    pA0 <- df_row$BFmean0A/100 #mean percentile
    spA0 <- df_row$BFsd0A/100 #SD percentile
    
    #Reference arm
    pB0 <- df_row$BFmean0B/100
    spB0 <- df_row$BFsd0B/100
    
    if(distribution=="beta"){
      #sample ps from a beta dist
      psA0 <- p.sample.beta(pA0, spA0, num_samps) 
      psB0 <- p.sample.beta(pB0, spB0, num_samps)
      
    }else if(distribution=="normal"){
      #sample ps from a normal dist
      psA0 <- p.sample.norm(pA0, spA0, num_samps) 
      psB0 <- p.sample.norm(pB0, spB0, num_samps)
      
    }else{
      print("Not a valid distribution")
    }
    
    #convert sampled ps to zs
    zsA0 <- qnorm(psA0)
    zsB0 <- qnorm(psB0)
    
    #save results
    df_row$BFmean0A <- mean(zsA0)
    df_row$BFsd0A <- sd(zsA0)
    
    df_row$BFmean0B <- mean(zsB0)
    df_row$BFsd0B <- sd(zsB0)
    
  }
  ## FOLLOW-UP
  if(df_row$test_type == "FU" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){

    #Intervention arm
    pA1 <- df_row$BFmean1A/100 #mean percentile
    spA1 <- df_row$BFsd1A/100 #SD percentile
    
    #Reference arm
    pB1 <- df_row$BFmean1B/100
    spB1 <- df_row$BFsd1B/100
    
    if(distribution=="beta"){
      #sample ps from a beta dist
      psA1 <- p.sample.beta(pA1, spA1, num_samps) 
      psB1 <- p.sample.beta(pB1, spB1, num_samps)
      
    }else if(distribution=="normal"){
      #sample ps from a normal dist
      psA1 <- p.sample.norm(pA1, spA1, num_samps) 
      psB1 <- p.sample.norm(pB1, spB1, num_samps)
      
    }else{
      print("Not a valid distribution")
    }
    
    #convert sampled ps to zs
    zsA1 <- qnorm(psA1)
    zsB1 <- qnorm(psB1)
    
    #save results
    df_row$BFmean1A <- mean(zsA1)
    df_row$BFsd1A <- sd(zsA1)
    
    df_row$BFmean1B <- mean(zsB1)
    df_row$BFsd1B <- sd(zsB1)
  }
  
  df_row$measure <- "BMI-z from percentile"
  df_row
  
}

#sample percentiles from a beta distribution
p.sample.beta <- function(mu, sd, num_samps){
  #mu = mean percentile, sd = SD percentile
  var <- sd^2
  
  #calculate beta and alpha from mean and variance
  alpha <- (mu^2)*(((1-mu)/var) - 1/mu)
  beta <- alpha*((1/mu) - 1)
  
  ps <- c()
  for(i in 1:num_samps){
    pit <- 10
    while(pit>=1 || pit<=0){ #make sure not sampling exactly 0 or 1
      pit <- rbeta(1, alpha, beta, ncp = 0)
    }
    ps <- append(ps, pit)
  }
  ps
}

#sample percentile from a normal distribution
p.sample.norm <- function(mu, sd, num_samps){
  #mu = mean percentile, sd = SD percentile
  ps <- c()
  for(i in 1:num_samps){
    pit <- 10
    while(pit>=1 || pit<=0){ #make sure not sampling exactly 0 or 1
      pit <- rnorm(1, mu, sd)
    }
    ps <- append(ps, pit)
  }
  ps
}



