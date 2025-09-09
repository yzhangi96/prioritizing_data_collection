## Import libraries
library(calinf) # for simulating random shift. See https://arxiv.org/abs/2202.11886
library(comprehenr) # for list comprehension
library(ggplot2) 
library(randomForest)
library(dplyr)
library(parallel)

# Set number of cores for parallizing experiment
num_cores <- detectCores() - 1

sample_data <-  function(deltas, n=1000){
  #' Sample data from distributions under random shifts
  #' 
  #' @param deltas A vector of two delta values to indicate 
  #' strength of shifts in P2 and P3.
  #' @param n Number of samples to generate.
  #' 
  #' @return list of dataframes, each containing samples from a distribution.
  
  # Fixed distribution P1
  X1 <- replicate(50, rnorm(n, mean=0))
  Y1 <- rowSums(X1)
  
  # Shifted distribution P2
  dseed1 <- distributional_seed(n=n,delta=deltas[1])
  X2 <- replicate(50, drnorm(dseed1,mean=0,sd=1))
  Y2 <- rowSums(X2)
  
  # Shifted distribution P3
  dseed2 <- distributional_seed(n=n,delta=deltas[2])
  X3 <- replicate(50, drnorm(dseed2,mean=0,sd=1))
  Y3 <- rowSums(X3)
  
  df_list <- list(data.frame("Y"=Y1, "X"=X1),
                  data.frame("Y"=Y2, "X"=X2),
                  data.frame("Y"=Y3, "X"=X3))
  return(df_list)
}


calc_weight <- function(P1, P2, P3){
  #' Calculate the weights for each distribution
  #' @param P1, P2, P3 Dataframes containing the samples from each distribution.
  #' 
  #' @return weight for P3

  ## Calculate weights
  p <- ncol(P1)
  
  # Get means of covariates for each distribution
  P1_means <- apply(P1[,2:p],2,mean)
  P2_means <- apply(P2[,2:p],2,mean)
  P3_means <- apply(P3[,2:p],2,mean)
  
  # Calculate the weight for P2
  X_bar <- data.frame("y" = P1_means- P2_means, "x3" = P3_means-P2_means)
  w3 <- min(coef(lm(y~x3-1, data=X_bar)), 1)
  w3 <- max(w3, 0)
  
  return(w3)
}

calc_mse <- function(P1, Pk, c_k, C,n){
  #' Calculate the mean squared error for a given budget using 
  #' data from one source distribution
  #' 
  #' @param P1 Dataframe containing samples from target P1
  #' @param Pk Dataframe containing samples from source Pk
  #' @param c_k The cost per sample from Pk
  #' @param C The total budget available for sampling
  #' 
  #' @return MSE from RF predictions
  #' 
  
  p <- ncol(P1)
  
  # Unpack dataframes
  Y1 <- P1$Y
  X1 <- P1[,2:p]
  Yk <- Pk$Y
  Xk <- Pk[,2:p]
  
  # Samples allowable from total budget
  n_k <- floor(C/c_k)
  
  # Sample from source distribution Pk
  indices <- sample.int(n, n_k)
  Xk <- Xk[indices,]
  Yk <- Yk[indices]
  
  # RF predictions and MSE
  Yk_hat <- predict(randomForest(x=Xk, y=Yk), X1)
  mse <- mean((Y1 - Yk_hat)^2)
  
  return(mse)
}

calc_mse2 <- function(df_list, C, c_ks, deltas,n){
  #' Calculate MSE from optimal sampling method
  #' 
  #' @param df_list list of data
  #' @param c_k The cost per sample from Pk
  #' @param C The total budget available for sampling
  #' @param deltas A vector of two delta values to indicate 
  #' strength of shifts in P2 and P3.
  #' @param n Number of samples to generate.
  #' 
  #' @return MSE from RF predictions
  #' 
  
  # Calculate optimal samples
  deltasq <- deltas^2
  sqrt_c2c3 <- sqrt(c_ks[1]*c_ks[2])
  
  n2 <- round((sqrt_c2c3 - deltasq[2] * C - deltasq[1]) / (-deltasq[2] * c_ks[1] - deltasq[1] * sqrt_c2c3),0)
  n3 <- floor((C-n2*c_ks[1])/c_ks[2])
  
  P1 <- df_list[[1]]
  P2_indices <- sample.int(n, n2)
  P3_indices <- sample.int(n, n3)
  P2_subset <- df_list[[2]][P2_indices,]
  P3_subset <- df_list[[3]][P3_indices,]
  
  p <- ncol(P1)
  Y2_hat <- predict(randomForest(x=P2_subset[,2:p], y=P2_subset$Y) , P1[,2:p])
  Y3_hat <- predict(randomForest(x=P3_subset[,2:p], y=P3_subset$Y), P1[,2:p])
  
  w3 <- calc_weight(P1, P2_subset, P3_subset)
  weighted_mean <- Y3_hat*w3 + Y2_hat*(1-w3)
  mse <- mean((P1$Y - weighted_mean)^2)
  return(list("mse" = mse, "w3" = w3, "n2" = n2, "n3" = n3))
}


## Setting
trials = 50
deltas <- c(8, 5) # strength of shift
c_ks <- c(5, 20) # cost per sample from P2 and P3
budget_seq <- seq(100, 1500, 10) # total budget considered
n=500 # number of samples


run_experiment <- function(deltas, c_ks, n=n){
  #' Run a single experiment (across different total budgets) per Fig 1
  #' 
  #' @param deltas A vector of two delta values to indicate
  #' strength of shifts in P2 and P3.
  #' @param c_ks A vector of two cost values for sampling from P2 and P3.
  #' @param n Number of samples to generate.
  #' 
  #' @return A matrix with MSE values for P2, P3, and the weighted method
  
  # Sample data
  df_list <- sample_data(deltas, n=n)
  
  # Calculate MSE from RF trained on P2 samples
  p2_mse <- to_vec(
    for(b in budget_seq) calc_mse(df_list[[1]], df_list[[2]], c_ks[1], b,n)
  )
  # Calculate MSE from RF trained on P3 samples
  p3_mse <- to_vec(
    for(b in budget_seq) calc_mse(df_list[[1]], df_list[[3]], c_ks[2], b,n)
  )
  
  # Calculate MSE from RF trained on optimal samples
  weighted_mse <- NULL
  w3 <- NULL
  n2_opt <- NULL
  n3_opt <- NULL
  
  for(b in budget_seq) {
    res <- calc_mse2(df_list, b, c_ks, deltas,n)
    weighted_mse <- c(weighted_mse, res$mse)
    w3 <- c(w3, res$w3)
    n2_opt <- c(n2_opt, res$n2)
    n3_opt <- c(n3_opt, res$n3)
  }
  return(cbind(p2_mse, p3_mse,weighted_mse, w3, n2_opt,n3_opt))
}

#### CAUTION: TAKES A LONG TIME TO RUN #### 
# # Run experiment in parallel
# results <- mclapply(1:trials, function(x) run_experiment(deltas, c_ks, n), mc.cores = num_cores)
# results_array <- array(unlist(results), dim = c(length(budget_seq), 6, trials))
# 
# # Reshape results for plotting
# means_results <- data.frame(apply(results_array, c(1, 2), mean))
# colnames(means_results) <- c("p2_mse", "p3_mse", "weighted_mse", "w3", "n2_opt", "n3_opt")
# 
# mse_df <- data.frame(
#   Budget = c(budget_seq, budget_seq, budget_seq),
#   MSE = c(means_results$p2_mse, means_results$p3_mse, means_results$weighted_mse),
#   Distribution = factor(rep(c("P1", "P2", "Weighted"),
#                             c(nrow(means_results), nrow(means_results), nrow(means_results))))
# )
# saveRDS(mse_df, "mse_df.RDS")
######## 

# Read saved results
mse_df <- readRDS("mse_df.RDS")

# Fig 1 plot
ggplot(mse_df, aes(x = Budget, y = MSE, color = Distribution)) +
  geom_smooth(se = FALSE, method = "loess") +
  labs(title = "Sampling Tradeoff Under Budget", x = "Total Budget ($)", y = "Mean Squared Error") +
  theme_minimal()+
  scale_color_discrete(
    labels = c("P2, low quality", "P3, high quality", "Our method")
  ) +
  scale_x_continuous(limits = c(min(budget_seq),max(budget_seq))) +
  theme(
    axis.line = element_line(colour = "black"),
    axis.ticks.x = element_line(colour = "black"),
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),  
    axis.text = element_text(size = 10), 
    legend.position = "bottom"
  ) # 500 x 400





