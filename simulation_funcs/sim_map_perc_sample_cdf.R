## functions to map percentile to zBMI using the analytical method

map.p.cdf <- function(df_row){
  
  #mean and SD P
  p <- df_row$MeanP
  sp <- df_row$SDP
  
  res <- nleqslv(c(1,1), sim_eqs, jac = NULL, p, sp, method = "Newton")
  
  #estimate mean and sdz
  mz.est <- res$x[1]
  sz.est <- abs(res$x[2])
  
  z.est <- c(mz.est, sz.est)
  z.est
}

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


## functions to map percentile to zBMI using the sampling method

map.p.sample <- function(df_row, distribution = "beta", num_samps = 10000){
  
  #mean & SD perc
  p <- df_row$MeanP 
  sp <- df_row$SDP 
  
  if(distribution=="beta"){
    #sample ps from a beta dist
    ps <- p.sample.beta(p, sp, num_samps) 
  }else if(distribution=="normal"){
    #sample ps from a normal dist
    ps <- p.sample.norm(p, sp, num_samps) 
  }else{
    print("Not a valid distribution")
  }
  
  #convert sampled ps to zs
  zs <- qnorm(ps)
  
  #re-define 
  mz.est<- mean(zs)
  sz.est <- sd(zs)
  
  z.est <- c(mz.est, sz.est)
  z.est
  
}


p.sample.beta <- function(mu, sd, num_samps){
  #mu = mean percentile, sd = SD percentile
  var <- sd^2
  #calculate beta and alpha
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

