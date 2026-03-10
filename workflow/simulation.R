#Code to reproduce simulation results 

# Load required libraries
library(expm)
library(Rcpp)
library(truncnorm)
library(BradleyTerryScalable)
library(PCWST)
# Source supporting functions for simulations
 source("support.R")
# source("functions_supp.R")
# source("nuclear_est.R")
# sourceCpp("Rcpp_functions.cpp")
# sourceCpp("Rcpp_nuclear.cpp")
# Change data generation settings
nsim <- 50 #nsim: number of simulations
#Get results for different K, corresponding to ranks from 2 to 20. Use 1:10 to reproduce the results from simulations.
for (K in 2:2){
  for (sparse_lv in 3:3){#Sparse level: 1(sparse), 2(less sparse), 3(Dense,constant)
    #Choose n, the number of subjects (500,1000,1500 and 2000)
    for ( n in c(500)){
      #set.seed for reproducible outcomes 
      set.seed(123)
      
      # Generate parameters based on the settings
      parameters <- pcwst_genpar.func( n, K)
      
      #Store true values
      Pi_star <- parameters$Pi
      T_max <- parameters$T_max
      # Initialize storage arrays and variables
      BTcoefstore = matrix(0, nrow = n, ncol = nsim)
      Pihatstore = array(0, dim = c(n,n,nsim))
      niterstore = rep(0, nsim)
      Improvementsobjstore = rep(0, nsim)
      objstore = rep(0, nsim)
      exitstore = rep(0, nsim)
      
      
      # Storing loss information
      Totalloss = rep(0, nsim)
      Totalloss_BT= rep(0,nsim)
      #Initiate time based on BT and proposed.
      BTtime <- 0
      Prposedtime <- 0
      #Start running simulations
      for (z in 1:nsim) {
        print(sprintf("%dth simulation", z))
        
        # Data generation. 
        Y= pcwst_gendata.func(T_max, Pi_star, n,  sparse_lv )
        
        #Set scale based on rank. In this particular setting we use 2*K since we know the true nuclear norm from our construction.
        scale <-2*K 
        
        # Estimation process
        #Proposed
        start.time <- Sys.time()
        out <- pcwst_est(Y, scale)
        end.time <- Sys.time()
        Prposedtime <- Prposedtime + end.time - start.time
        exitstore[z] = out$exit
        
        niterstore[z] = out$iteration
        objstore[z] = out$objective
        
        Pihatstore[,,z] <- g.func(out$M_mat)
        Totalloss[z] <- loss.func(Pihatstore[,,z],Pi_star,n)
        
        #BT
        start.time <- Sys.time()
        BT_object <- BT_est.func(Y) 
        end.time <- Sys.time()
        BTtime <-  BTtime + end.time - start.time
        Pihat_BT <- BT_object$Pihat_bt
        BTcoefstore[,z] <- BT_object$coeff_bt
        Totalloss_BT[z] <- loss.func(Pihat_BT,Pi_star,n)
      }
      
      #Now generate test data and compute predictive likelihood
      #Save data 
      out <- list(Prposedtime = Prposedtime, BTtime = BTtime, Totalloss_BT =Totalloss_BT, Totalloss= Totalloss, objstore = objstore, niterstore = niterstore, exitstore = exitstore,  BTcoefstore=BTcoefstore, Pi_star = Pi_star ,Pihatstore = Pihatstore  )
      save(out, file = sprintf("slv%dnsim%dn%dK%d.Rdata",sparse_lv, nsim,n,K))
    }
  }
}