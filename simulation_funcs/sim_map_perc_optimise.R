## functions to map percentile to zBMI using the optimization method

map.perc.optimise <- function(df_row, num_samps, muz_init=0, sigz_init=1,
                              max_it=5000, perc_thresh, z_step){
  #num_samps is the number of times we sample from the z distribution
  #max_it is the maximum number of iterations allowed in the optimisation
  #muz_init and sigz_init are initial guesses for mean/sd of z distribution
  #perc_thresh is the threshold at which we stop the optimisation
  #z_step is the step size for adjusting the mean/sd of z distribution in optimisation
  
  
  #generate initial sample
  zs <- rnorm(num_samps, muz_init, sigz_init)
  ps <- pnorm(zs)
  
  #initialised mean/SD perc
  meanp_init <- mean(ps)
  sdp_init <- sd(ps)
  
  #Reported perc values (target)
  meanp_target <- df_row$MeanP
  sdp_target <- df_row$SDP
  
  results <- perc.optim.routine(meanp_init, sdp_init, meanp_target, sdp_target,
                                zs, muz_init, sigz_init, max_it, perc_thresh, z_step)
  
  #results -> zdist, zsamp, bsamp
  est.mzdist <- results[[1]] #muz
  est.szdist <- results[[2]] #sigz
  
  est.mzsamp <- results[[3]] #meanz
  est.szsamp <- results[[4]] #sdz
  
  est.mpsamp <- results[[5]] #meanp
  est.spsamp <- results[[6]] #sdp
  
  if(abs(results[[5]]-meanp_target)<= perc_thresh && abs(results[[6]]-sdp_target)<= perc_thresh){
    converge <- "Y"
  }else{
    converge <- "N"
  }
  
  ests <- list(c(est.mzdist,est.szdist), c(est.mzsamp, est.szsamp), c(est.mpsamp, est.spsamp), converge)  
  ests  
  
}

perc.optim.routine <- function(meanp_init, sdp_init, meanp_target, sdp_target,
                               zs, muz_init, sigz_init, max_it, perc_thresh, z_step){
  
  #initialise muz, sigz (parameters of z distribution)
  muz <- muz_init
  sigz <- sigz_init
  
  #initialise meanp, sdp
  meanp <- meanp_init
  sdp <- sdp_init
  
  
  num_it <- 0
  while( (abs(meanp-meanp_target) > perc_thresh || abs(sdp - sdp_target) > perc_thresh) && num_it < max_it){
    
    #CHOOSE DIRECTION OF CHANGE--------------------------
    #choose delta mu based on mean percentile
    if(abs(meanp-meanp_target) <= perc_thresh){ #if mean p is within threshold - fix muz
      del_mu <- 0.0
    }else if(meanp > meanp_target){ #if mean p is too big - reduce mu z
      del_mu <- -z_step #negative
    }else if(meanp < meanp_target){ #if mean p is too small - increase mu z
      del_mu <- z_step #positive
    }
    
    #choose delta sigma based on SD p
    if(abs(sdp-sdp_target) <= perc_thresh){ #if sd p is within threshold - fix sigz
      del_sig <- 0.0
    }else if(sdp > sdp_target){ #if sd p is too big - reduce sig z
      del_sig <- -z_step #negative
    }else if(sdp < sdp_target){ #if mean p is too small - increase mu z
      del_sig <- z_step #positive
    }
    
    #UPDATE PARAMETERS-----------------------------------
    #update zs
    #zi(t+1) = (1 + delsig/sig(t))*(zi(t)-mu(t)) + (mu(t)+delmu)
    for(i in 1:length(zs)){
      zs[i] <- (1 + del_sig/sigz)*(zs[i] - muz) + (muz + del_mu)
    }
    
    #update mu and sigma z
    muz <- muz + del_mu
    sigz <- sigz + del_sig
    
    
    #calculate mean/sd of sampled zs
    meanz <- mean(zs, na.rm = TRUE)
    sdz <- sd(zs, na.rm = TRUE)
    
    
    #update ps (using new zs)
    ps <- pnorm(zs)
    
    #update mean/sd bmi
    meanp <- mean(ps, na.rm = TRUE)
    sdp <- sd(ps, na.rm = TRUE)
    
    num_it <- num_it + 1
    
  }
  
  result <- list(muz, sigz, meanz, sdz, meanp, sdp)
  result
}
