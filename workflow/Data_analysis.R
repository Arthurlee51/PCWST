#Data_analysis with tuning. 50% training, 20%validation and 30%testing
# Main script to run simulations
# --------------------------------------------------------------------------------

# Load required libraries
# Load required libraries
library(expm)
library(Rcpp)
library(RcppEigen)
library(truncnorm)
library(BradleyTerryScalable)
library(PCWST)

# Source functions from external R scripts
source("support.R")
sourceCpp("rcpp_support.cpp")

#load data(ATP for HotS)
example = "ATP"
if(example =="ATP"){
  load("ATPdata.Rdata")
  data <- rm_extreme_data.func(full_data)
  #Remove full_data
  rm(full_data)
}

if(example =="HotS"){#WoL data from SC2
  load("HotS.Rdata")
  data <- rm_extreme_data.func(HotSdata)
}

#======================================================================================
#Analysis part.
numScales <- 20
scaleVals <- 10^seq(-1, 1, length.out = numScales)#Potential values of C\sqrt{2K} for comparison
set.seed(123)
nmethod <- 2
#split_ratio <- 0.5 #ratio of training data. The rest is left as test data. #Tricky in the sense that missing rate is too high.
ndata<- nrow(data)
n_train <- round(ndata*0.5)
n_valid <- round(ndata*0.2)#Number of items for validation

#Get training, validation and test data.
  train_ind <-  sample(1:ndata,n_train,replace = FALSE)
  data_train <-data[train_ind,]
  
  data_rest <- data[-train_ind, ]
  n_rest <- nrow(data_rest)
  #Get validation set 
  valid_ind <- sample(1:n_rest,n_valid,replace = FALSE)
  data_valid <- data_rest[valid_ind,]
  data_test <-  data_rest[-valid_ind,]
  
  #Drop extreme data in training set. 
  data_train <- rm_extreme_data.func(data_train)
  
  #Get Player and Player id from training data and form W. Note that it is different in each split.
  players_train <- unique(data_train$Winner)
  n_players_train  <- length(players_train)
  #Get W_train
  W_train <- GetW.func(data_train, n_players_train,   players_train)
  W_train <- as.matrix(W_train)
  
  
  #Get W_valid
  data_valid_final <- data_valid[which(data_valid$Winner %in%  players_train & data_valid$Loser %in%  players_train),]
  #compute number of players removed 
  nremoved_valid <- nrow(data_valid) - nrow(data_valid_final)
  #Now create W_test for easier analysis. Only consider players that appeared in te training set. Since otherwise the prediction is just 0.5
  W_valid <- GetW.func(data_valid_final, n_players_train,   players_train )
  
  
  #Estimate based on proposed method
  start_time <- Sys.time()
  out_full <- pcwst_est_tune(W_train, W_valid,scaleVals)
  end_time <- Sys.time()
  time_proposed <-  end_time -  start_time 
  #out <- out_full$out
  scale <- out_full$scale
  valid_like <- out_full$valid_like
  
  #Combine W_valid and W_train and fit the model again.
  start_time <- Sys.time()
  W_comb <- W_train + W_valid
  out_comb <- pcwst_est(W_comb, scale)
  Mhat <- g.func(out_comb$M_mat)
  end_time <- Sys.time()
  time_proposed_comb <-  end_time -  start_time 
  #BT version
  start_time <- Sys.time()
  BT_object <- BT_est.func(W_comb)
  end_time <- Sys.time()
  time_BT <-  end_time -  start_time 
  Mhat_BT <- BT_object$Pihat_bt
  
  #Evaluation part.
  data_test_final <- data_test[which(data_test$Winner %in%  players_train & data_test$Loser %in%  players_train),]
  #compute number of players removed 
  nremoved <- nrow(data_test) - nrow(data_test_final)
  #Now create W_test for easier analysis. Only consider players that appeared in te training set. Since otherwise the prediction is just 0.5
  W_test <- GetW.func(data_test_final, n_players_train,   players_train )
  

  #like_loss: loss based on likelihood calculations. Needs scaling (n_players_train^2) to get back the result in manuscript.
  like_loss <- rcpp_like(W_test, Mhat, n_players_train)
  like_loss_BT <- rcpp_like(W_test, Mhat_BT, n_players_train)
  like_losses <- c( like_loss,like_loss_BT )
  
  #acc: accuracy 
  acc <- rcpp_acc(W_test, Mhat, n_players_train)
  acc_BT <- rcpp_acc(W_test, Mhat_BT, n_players_train)
  acc_s <- c(acc,  acc_BT)
  save(time_proposed,time_BT, time_proposed_comb,example, n_players_train, like_losses, acc_s,  scale, out_comb,valid_like=valid_like,   file = sprintf("dataanalysis%sn%dnumscale%d.Rdata", example, n_players_train,numScales))
 
