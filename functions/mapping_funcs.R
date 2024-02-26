#Function to create dataset on which to perform the mapping tests
create_ds_maptest <- function(df, icc, rho){
  #subset by MAP_test (include in mapping test)
  df <- subset(df, MAP_test == "Y")
  
  #exclude Adjusted FU
  df <- subset(df, test_type != "Adjusted FU")
  
  #exclude MD only
  df <- subset(df, test_type != "MD only")
  
  for(i in 1:nrow(df)){
    
    #unadjust for clustering-----------------------------------------------------
    if(df$UnAdj_cluster[i]!="N" & !is.na(df$UnAdj_cluster[i])){
      
      #Define ICC
      if(df$calcICC[i]=="Y"){ # If Calculated ICC = Y use ICC (reported)
        del <- df$ICC[i]
      }else if(df$calcICC[i]=="N"){ # Else if Calculated ICC= N use icc (imputed)
        del <- icc
      }else{print(paste("CalcICC error:", df$study[i]))}
      
      #Check nAc and nBc (no. of clusters per group) are defined
      if(is.na(df$nAc[i]) || is.na(df$nBc[i])){print(paste("n cluster not defined:", df$study[i]))}
      
      #Calculate mean cluster size-------------------------------------------
      nA1 <- df$nA1[i] #no. of participants in group A at follow-up
      nB1 <- df$nB1[i] #no. of participants in group B at follow-up
      nA0 <- df$nA0[i] #no. of participants in group A at baseline
      nB0 <- df$nB0[i] #no. of participants in group B at baseline
      nAc <- as.numeric(df$nAc[i]) #no. of clusters in group A
      nBc <- as.numeric(df$nBc[i]) #no. of clusters in group B
      mB <- nB1/nBc
      mA <- nA1/nAc
      mB0<- nB0/nBc
      mA0 <- nA0/nAc
      
      #Un-adjust SDs---------------------------------------------------
      #sqrt(design effect) = SQRT(1+icc(m-1))
      desA1 <- sqrt(1+del*(mA-1)) 
      desB1 <- sqrt(1+del*(mB-1))
      desA0 <- sqrt(1+del*(mA0-1)) 
      desB0 <- sqrt(1+del*(mB0-1))
      
      #adjusted SD: SD' = SD*sqrt(design effect)
      #un-adjusted SD: SD = SD'/sqrt(design effect)
      
      #unadjust CS
      if(df$UnAdj_cluster[i]=="CS"||df$UnAdj_cluster[i]=="Baseline and CS"){
        df$CSsdA[i] <- df$CSsdA[i]/desA1
        df$CSsdB[i] <- df$CSsdB[i]/desB1
      }
      #unadjust baseline
      if(df$UnAdj_cluster[i]=="Baseline" || df$UnAdj_cluster[i]=="Baseline and FU" || df$UnAdj_cluster[i]=="Baseline and CS"){
        df$BFsd0A[i] <- df$BFsd0A[i]/desA0
        df$BFsd0B[i] <- df$BFsd0B[i]/desB0
      }
      #unadjust FU
      if(df$UnAdj_cluster[i]=="FU"||df$UnAdj_cluster[i]=="Baseline and FU"){
        df$BFsd1A[i] <- df$BFsd1A[i]/desA1
        df$BFsd1B[i] <- df$BFsd1B[i]/desB1
      }
      
    }#end unadjust for clustering
    
    #calculate baseline & FU from CS--------------------------------------------
    if(df$test_type[i]=="Baseline and CS"){
      #FU = baseline + CS
      df$BFmean1A[i] <- df$BFmean0A[i] + df$CSmeanA[i]
      df$BFmean1B[i] <- df$BFmean0B[i] + df$CSmeanB[i]
      
      #set FU SD = baseline SD
      df$BFsd1A[i] <- df$BFsd0A[i] 
      df$BFsd1B[i] <- df$BFsd0B[i] 
      
    }#end calc FU
    
  }
  

  df <- df[, c("study", "studyID", "measure", "time", "Aarm", "Barm","nA0", "nA1","nB0", "nB1", "test_type", 
           "BFmean0A", "BFsd0A","BFmean0B", "BFsd0B","BFmean1A", "BFsd1A","BFmean1B", "BFsd1B",
           "Mean_Age", "SD_Age", "Prop_Male", "UnAdj_cluster", "fu_months", "chart", "chart_ref")]
  df
}


#function to format results
arm_vec <- function(data){
  cols <- c("studyID", "arm", "time", "mean", "SD", "n", "test_type", "UnAdj_cluster", "chart_ref", "chart", "converge")
  new_df <- data.frame(matrix(nrow = 0, ncol = length(cols)))
  
  for(i in 1:nrow(data)){
    if(data$test_type[i]=="Baseline and FU" || data$test_type[i]=="Baseline and CS"){
      
      row0A <- c(data$studyID[i], data$Aarm[i], "baseline", data$BFmean0A[i], data$BFsd0A[i], data$nA0[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeA0[i])
      row1A <- c(data$studyID[i], data$Aarm[i], data$time[i], data$BFmean1A[i], data$BFsd1A[i], data$nA1[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeA1[i])
      row0B <- c(data$studyID[i], data$Barm[i], "baseline", data$BFmean0B[i], data$BFsd0B[i], data$nB0[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeB0[i])
      row1B <- c(data$studyID[i], data$Barm[i], data$time[i], data$BFmean1B[i], data$BFsd1B[i], data$nB1[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeB1[i])
      
      new_df <- rbind(new_df, row0A)
      new_df <- rbind(new_df, row1A)
      new_df <- rbind(new_df, row0B)
      new_df <- rbind(new_df, row1B)
      
    }else if(data$test_type[i]=="Baseline"){
      
      row0A <- c(data$studyID[i], data$Aarm[i], "baseline", data$BFmean0A[i], data$BFsd0A[i], data$nA0[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeA0[i])
      row0B <- c(data$studyID[i], data$Barm[i], "baseline", data$BFmean0B[i], data$BFsd0B[i], data$nB0[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeB0[i])
      
      new_df <- rbind(new_df, row0A)
      new_df <- rbind(new_df, row0B)
      
    }else if(data$test_type[i]=="FU"){
      
      row1A <- c(data$studyID[i], data$Aarm[i], data$time[i], data$BFmean1A[i], data$BFsd1A[i], data$nA1[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeA1[i])
      row1B <- c(data$studyID[i], data$Barm[i], data$time[i], data$BFmean1B[i], data$BFsd1B[i], data$nB1[i], data$test_type[i], data$UnAdj_cluster[i], data$chart_ref[i], data$chart[i], data$convergeB1[i])
      
      new_df <- rbind(new_df, row1A)
      new_df <- rbind(new_df, row1B)
      
    }else{
      print(paste("Error test type", data$studyID[i]))
    }
    
  }
  
  colnames(new_df) <- cols
  
  #remove duplicate rows
  new_df <- new_df %>% distinct(studyID, arm, time, n, .keep_all = TRUE) 
  
  new_df
  
}

#function to only keep rows of df1 that appear in df2
comp_dfs <- function(df1, df2){
  df1$uniq <- paste(df1$studyID, df1$arm, df1$time)
  df2$uniq <- paste(df2$studyID, df2$arm, df2$time)
  
  df1 <- df1[df1$uniq %in% df2$uniq, ]
  
  df1 <- subset(df1, select = -uniq)
  df1
}

