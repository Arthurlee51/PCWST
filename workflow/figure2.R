#Script to generate figure 2, where we compute the predictive log-likelihood on newly generated test set before plotting. 
library(Rcpp)
library(PCWST)
source("support.R")

#Get predictive likelihood
nlong <- c(500, 1000, 1500, 2000)
sparse_long <- c( 1,2, 3)
K_ranges <- 1:10
t=0
nsim <- 50
T_max <- 5
len_report <- length(nlong) * length(sparse_long) * length(K_ranges)
lik_results <- list()
for (sparse_lv in sparse_long) {
  for (n in nlong) {
    for (k in K_ranges) {
      print(t)
      t <- t + 1
      load(sprintf("slv%dnsim%dn%dK%d.Rdata", sparse_lv, nsim, n, k))
      #load BT and proposed parameters 
      BTcoefstore <- out$BTcoefstore
      Pihatstore <- out$Pihatstore
      #After loading, create new data to evaluate the predictive log-likelihood
      #Get the correct Pistar
      Pi_star <- out$Pi_star
      #set.seed for reproducible outcomes 
      set.seed(1234)#use a different seed to make sure that test data is independent to training data.
      lik_BT <- rep(0, nsim)
      lik_proposed <- rep(0, nsim)
      for(z in 1:nsim){
        Y_test = pcwst_gendata.func(T_max, Pi_star, n,  sparse_lv, test =TRUE)
        #Now compute the predictive likelihood on Y to evaluate the performance of the method 
        #BT
        BTcoef <- BTcoefstore[,z]
        #Construct Pihat based on the BTcoef
        Pihat_BT <- g.func(matrix( BTcoef  ,nrow=n, ncol=n ) - matrix(  BTcoef  ,nrow=n, ncol=n, byrow=TRUE))
        #Get BTcoef
        #Proposed
        Pihat_proposed <- Pihatstore[,,z]
        
        count <- sum(Y_test)
        #Compute predictive likelihood
        lik_BT[z] <-  PCWST:::rcpp_like(Y_test,  Pihat_BT, n)/  count 
        lik_proposed[z] <- PCWST:::rcpp_like(Y_test,  Pihat_proposed, n)  /count 
      }
      
      # store everything needed for plotting
      lik_results[[t]] <- list(
        sparse_lv = sparse_lv,
        n = n,
        K = k,
        lik_BT = lik_BT,
        lik_proposed = lik_proposed, 
        totaltime =  test_time
      )
    }
  }
}

#======================================================================================
#Plotting part 
library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(ggtext)

data <- c()
nsim <- 50
nlong <- c(500, 1000, 1500, 2000)
sparse_long <- c(1, 2, 3)
K_ranges <- 1:10

# Aim: Plot lik for various levels of K.
# Final figure: 2 x 4 panels (rows = Method, cols = n)

len_report <- length(nlong) * length(sparse_long) * length(K_ranges)

# ---- report stores mean and sd for both methods (7 columns) ----
report <- data.frame(matrix(0, nrow = len_report, ncol = 7))
colnames(report) <- c("K", "n", "sparse_lv",
                      "lik_mean", "lik_sd",
                      "lik_BT_mean", "lik_BT_sd")

z <- 0

for (sparse_lv in sparse_long) {
  for (n in nlong) {
    for (k in K_ranges) {
      z <- z + 1
      result <- lik_results[[z]]
      report[z, ] <- c(
        result$K, result$n, result$sparse_lv,
        mean(result$lik_proposed ),     sd(result$lik_proposed),
        mean(result$lik_BT),  sd(result$lik_BT)
      )
    }
  }
}

# Step 1: Reshape to long format, keeping mean+sd together
report_long <- report %>%
  pivot_longer(
    cols = c("lik_mean", "lik_BT_mean", "lik_sd", "lik_BT_sd"),
    names_to = c("Method", ".value"),
    names_pattern = "(lik|lik_BT)_(mean|sd)"
  )

# Step 2: Update factor levels and labels + compute mean ± 2SD bounds
report_long <- report_long %>%
  mutate(
    K = as.integer(K),
    n = factor(n, levels = c(500, 1000, 1500, 2000)),
    sparse_lv = factor(sparse_lv),
    Method = factor(Method, levels = c("lik", "lik_BT"),
                    labels = c("Proposed", "BT")),
    lik = mean,
    lik_low  = mean - 2 * sd,
    lik_high = mean + 2 * sd
  )

#poltting
ggplot(report_long,
       aes(x = K, y = lik,
           color = sparse_lv,
           linetype = sparse_lv,
           group = sparse_lv)) +
  
  geom_line(linewidth = 0.7) +
  
  geom_errorbar(
    aes(ymin = lik_low, ymax = lik_high),
    width = 0.30,
    linewidth = 1.0,
    linetype = "solid",     
    show.legend = FALSE
  ) +
  
  facet_grid(Method ~ n,
             labeller = labeller(n = function(x) paste0("n = ", x))) +
  
  labs(
    title = "Predictive Likelihood Comparison by Sparsity Level and Sample Size",
    x = "K",
    y = "Mean Likelihood",
    color = "Sparsity Level",
    linetype = "Sparsity Level"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(size = 12),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.spacing = unit(1.1, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  
  scale_color_manual(
    values = c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3"),
    labels = c("Sparse", "Less Sparse", "Dense")
  ) +
  
  scale_linetype_manual(
    values = c("1" = "dotted", "2" = "dashed", "3" = "solid"),
    labels = c("Sparse", "Less Sparse", "Dense")
  ) +
  
  scale_x_continuous(
    breaks = K_ranges,
    minor_breaks = NULL
  )
