#Function to map a row of percentile data onto zBMI scale via optimization method
map.perc.optimise <- function(df_row, num_samps, muz_init=0, sigz_init=1,
                              max_it=5000, perc_thresh, z_step){
  #num_samps is the number of times we sample from the z distribution
  #max_it is the maximum number of iterations allowed in the optimisation
  #muz_init and sigz_init are initial guesses for mean/sd of z distribution
  #perc_thresh is the threshold at which we stop the optimisation (Delta-tol)
  #z_step is the step size for adjusting the mean/sd of z distribution in optimisation (Delta-step)
  
  #3 df rows will be output: 
  #z_dist (final mu_z, sigma_z) = distribution estimates
  #z_sample (final mean(z), sd(z) from sample) = sample estimates
  #p_sample (final mean(p), sd(p) from sample) = corresponding values of percentile sampled during optimization
  #initialise each as = df_row
  df_row_zdist <- df_row
  df_row_zsamp <- df_row
  df_row_psamp <- df_row
  
  ## For each input row we map baseline & FU values for the intervention and reference arm (where available)
  
  #BASELINE=====================================================================
  if(df_row$test_type == "Baseline" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    
    #Intervention arm---------------------------------------------------------------
    #generate initial sample
    zsA0 <- rnorm(num_samps, muz_init, sigz_init)
    psA0 <- pnorm(zsA0)
    
    #initialised mean/SD perc
    meanp_initA0 <- mean(psA0)
    sdp_initA0 <- sd(psA0)
    
    #Reported perc values (target)
    meanp_targetA0 <- df_row$BFmean0A/100 #mean percentile
    sdp_targetA0 <- df_row$BFsd0A/100 #SD percentile
    
    #optimise
    resultsA0 <- perc.optim.routine(meanp_initA0, sdp_initA0, meanp_targetA0, sdp_targetA0,
                                    zsA0, muz_init, sigz_init, max_it, perc_thresh, z_step)
    
    #use results to re-define means/SDs in zdist, zsamp, and psamp
    df_row_zdist$BFmean0A <- resultsA0[[1]] #muz
    df_row_zdist$BFsd0A <- resultsA0[[2]] #sigz
    
    df_row_zsamp$BFmean0A <- resultsA0[[3]] #meanz
    df_row_zsamp$BFsd0A <- resultsA0[[4]] #sdz
    
    df_row_psamp$BFmean0A <- resultsA0[[5]] #meanp
    df_row_psamp$BFsd0A <- resultsA0[[6]] #sdp
    
    #check if optimization converged:
    if(abs(resultsA0[[5]]-meanp_targetA0)<= perc_thresh && abs(resultsA0[[6]]-sdp_targetA0)<= perc_thresh){
      df_row_zdist$convergeA0 <- "Y"
      df_row_zsamp$convergeA0 <- "Y"
      df_row_psamp$convergeA0 <- "Y"
    }else{
      df_row_zdist$convergeA0 <- "N"
      df_row_zsamp$convergeA0 <- "N"
      df_row_psamp$convergeA0 <- "N"
    }
    
    
    #Reference arm------------------------------------------------------------------
    #generate initial sample
    zsB0 <- rnorm(num_samps, muz_init, sigz_init)
    psB0 <- pnorm(zsB0)
    
    #initialised mean/SD perc
    meanp_initB0 <- mean(psB0)
    sdp_initB0 <- sd(psB0)
    
    #Reported BMI values (target)
    meanp_targetB0 <- df_row$BFmean0B/100 #mean percentile
    sdp_targetB0 <- df_row$BFsd0B/100 #SD percentile
    
    #optimise
    resultsB0 <- perc.optim.routine(meanp_initB0, sdp_initB0, meanp_targetB0, sdp_targetB0,
                                    zsB0, muz_init, sigz_init, max_it, perc_thresh, z_step)
    
    #use results to re-define means/SDs in zdist, zsamp, and psamp
    df_row_zdist$BFmean0B <- resultsB0[[1]] #muz
    df_row_zdist$BFsd0B <- resultsB0[[2]] #sigz
    
    df_row_zsamp$BFmean0B <- resultsB0[[3]] #meanz
    df_row_zsamp$BFsd0B <- resultsB0[[4]] #sdz
    
    df_row_psamp$BFmean0B <- resultsB0[[5]] #meanp
    df_row_psamp$BFsd0B <- resultsB0[[6]] #sdp
    
    #check if optimization converged:
    if(abs(resultsB0[[5]]-meanp_targetB0)<= perc_thresh && abs(resultsB0[[6]]-sdp_targetB0)<= perc_thresh){
      df_row_zdist$convergeB0 <- "Y"
      df_row_zsamp$convergeB0 <- "Y"
      df_row_psamp$convergeB0 <- "Y"
    }else{
      df_row_zdist$convergeB0 <- "N"
      df_row_zsamp$convergeB0 <- "N"
      df_row_psamp$convergeB0 <- "N"
    }
    
  }#end mapping baseline values
  
  #FOLLOW-UP====================================================================
  if(df_row$test_type == "FU" || df_row$test_type == "Baseline and FU" || df_row$test_type == "Baseline and CS"){
    #Intervention arm---------------------------------------------------------------
    #generate initial sample
    zsA1 <- rnorm(num_samps, muz_init, sigz_init)
    psA1 <- pnorm(zsA1)
    
    #initialised mean/SD perc
    meanp_initA1 <- mean(psA1)
    sdp_initA1 <- sd(psA1)
    
    #Reported BMI values (target)
    meanp_targetA1 <- df_row$BFmean1A/100 #mean percentile
    sdp_targetA1 <- df_row$BFsd1A/100 #SD percentile
    
    #optimise
    resultsA1 <- perc.optim.routine(meanp_initA1, sdp_initA1, meanp_targetA1, sdp_targetA1,
                                    zsA1, muz_init, sigz_init, max_it, perc_thresh, z_step)
    
    #use results to re-define means/SDs in zdist, zsamp, and psamp
    df_row_zdist$BFmean1A <- resultsA1[[1]] #muz
    df_row_zdist$BFsd1A <- resultsA1[[2]] #sigz
    
    df_row_zsamp$BFmean1A <- resultsA1[[3]] #meanz
    df_row_zsamp$BFsd1A <- resultsA1[[4]] #sdz
    
    df_row_psamp$BFmean1A <- resultsA1[[5]] #meanp
    df_row_psamp$BFsd1A <- resultsA1[[6]] #sdp
    
    #check if optimization converged:
    if(abs(resultsA1[[5]]-meanp_targetA1)<= perc_thresh && abs(resultsA1[[6]]-sdp_targetA1)<= perc_thresh){
      df_row_zdist$convergeA1 <- "Y"
      df_row_zsamp$convergeA1 <- "Y"
      df_row_psamp$convergeA1 <- "Y"
    }else{
      df_row_zdist$convergeA1 <- "N"
      df_row_zsamp$convergeA1 <- "N"
      df_row_psamp$convergeA1 <- "N"
    }
    
    #Reference arm------------------------------------------------------------------
    #generate initial sample
    zsB1 <- rnorm(num_samps, muz_init, sigz_init)
    psB1 <- pnorm(zsB1)
    
    #initialised mean/SD perc
    meanp_initB1 <- mean(psB1)
    sdp_initB1 <- sd(psB1)
    
    #Reported BMI values (target)
    meanp_targetB1 <- df_row$BFmean1B/100 #mean percentile
    sdp_targetB1 <- df_row$BFsd1B/100 #SD percentile
    
    #optimise
    resultsB1 <- perc.optim.routine(meanp_initB1, sdp_initB1, meanp_targetB1, sdp_targetB1,
                                    zsB1, muz_init, sigz_init, max_it, perc_thresh, z_step)
    
    #use results to re-define means/SDs in zdist, zsamp, and psamp
    df_row_zdist$BFmean1B <- resultsB1[[1]] #muz
    df_row_zdist$BFsd1B <- resultsB1[[2]] #sigz
    
    df_row_zsamp$BFmean1B <- resultsB1[[3]] #meanz
    df_row_zsamp$BFsd1B <- resultsB1[[4]] #sdz
    
    df_row_psamp$BFmean1B <- resultsB1[[5]] #meanp
    df_row_psamp$BFsd1B <- resultsB1[[6]] #sdp
    
    #check if optimization converged:
    if(abs(resultsB1[[5]]-meanp_targetB1)<= perc_thresh && abs(resultsB1[[6]]-sdp_targetB1)<= perc_thresh){
      df_row_zdist$convergeB1 <- "Y"
      df_row_zsamp$convergeB1 <- "Y"
      df_row_psamp$convergeB1 <- "Y"
    }else{
      df_row_zdist$convergeB1 <- "N"
      df_row_zsamp$convergeB1 <- "N"
      df_row_psamp$convergeB1 <- "N"
    }
    
  }#end mapping FU values
  
  
  df_row_zdist$measure <- "BMI-z from percentile (dist)"
  df_row_zsamp$measure <- "BMI-z from percentile (samp)"
  df_row_psamp$measure <- "percentile from percentile (check)" #may need to multiple by 100
  df_rows <- list(df_row_zdist, df_row_zsamp, df_row_psamp)
  df_rows
}

#function to perform the optimization routine for percentile
perc.optim.routine <- function(meanp_init, sdp_init, meanp_target, sdp_target,
                               zs, muz_init, sigz_init, max_it, perc_thresh, z_step){
  
  #initialise muz, sigz (parameters of z distribution)
  muz <- muz_init
  sigz <- sigz_init
  
  #initialise meanp, sdp
  meanp <- meanp_init
  sdp <- sdp_init
  
  num_it <- 0 #track number of iterations
  #while not converged
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
    
    #update mean/sd percentile
    meanp <- mean(ps, na.rm = TRUE)
    sdp <- sd(ps, na.rm = TRUE)
    
    num_it <- num_it + 1
    
  }
  
  result <- list(muz, sigz, meanz, sdz, meanp, sdp)
  result
}

