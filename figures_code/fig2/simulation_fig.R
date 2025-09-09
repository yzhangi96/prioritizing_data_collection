## Import libraries
library(calinf) # for simulating random shift. See https://arxiv.org/abs/2202.11886
library(comprehenr) # for list comprehension
library(ggplot2) 
library(dplyr)
library(parallel)

get_P1 <- function(n, d=15){
  #' Get samples from target distribution P1
  #' @param n Number of samples to generate
  #' @param d 1/2 Number of covariates 
  #' @return Dataframe containing samples from target distribution P1

  # Generate X1–Xd ~ Normal(0,1)
  X_norm <- replicate(d, rnorm(n, mean=0, sd=1))
  # Generate X_d+1–X_2d ~ Bern(0.5)
  X_bin <- replicate(d, rbinom(n, size=1, prob=0.5))
  noise <- rnorm(n, mean=0, sd=1)

  X <- cbind(X_norm, X_bin)
  # Apply nonlinear transforms
  X[,1] <- X[,1]^2   # X1^2
  X[,2] <- X[,2]^2   # X2^2
  Y <- rowSums(X) + noise
  df <- data.frame(Y = Y, X)
  return(df)
}

get_random_shifts <- function(delta,n, d=15){
  #' Get randomly shifted X,y for a given delta
  #' @param delta Strength of shift
  #' @param n Number of samples to generate
  #' @param d 1/2 Number of covariates 
  #' @return Dataframe containing samples from a randomly shifted distribution
  #' 

  dseed <- distributional_seed(n=n,delta=delta)
  # Generate X1–Xd ~ Normal(0,1)
  X_norm <- replicate(d, drnorm(dseed, mean=0, sd=1))
  # Generate X_d+1–X_2d ~ Bern(0.5)
  X_bin <- replicate(d, drbinom(dseed, size=1, prob=0.5))
  noise <- drnorm(dseed, mean=0, sd=1)
  
  X <- cbind(X_norm, X_bin)
  # Apply nonlinear transforms
  X[,1] <- X[,1]^2   # X1^2
  X[,2] <- X[,2]^2   # X2^2
  Y <- rowSums(X) + noise
  df <- data.frame(Y = Y, X)
  return(df)
}

sample_data <- function(deltas, nk_vec) {
  #' Sample data from distributions under random shifts
  #' 
  #' @param deltas A vector of distr_var values to indicate 
  #' strength of shifts 
  #' 
  #' @return list of dataframes, each containing samples from a distribution.
  
  ## Sample Pk
  source_list <- list()
  for(i in 1:length(deltas)){
    source_list[[i]] <- get_random_shifts(deltas[i], nk_vec[i])
  }
  names(source_list) <- paste0("P", 1:length(deltas))
  
  return(source_list)
}


calc_weights_and_duc <- function(P1_train, P2_train, Pk_train, case=3){
  #' Calculate the optimal weights and duc for each distribution
  #' @param P1_train Subset of samples from target distribution P1
  #' @param P2_train Subset of samples from P2
  #' @param Pk_train Subset of samples from Pk
  #' @param case Case 3 refers to 3 distributions (P1, P2, Pk), 
  #' case 2 refers to 2 distributions (P1, P2)
  #' 
  #' @return list of weight vector and duc
  
  p <- ncol(P1_train)
  P1_pop_means <- c(c(1,1), rep(0, (p-1)/2-2), rep(0.5, (p-1)/2)) #population means

  # Get means of covariates for each distribution
  P1_samp_means  <- apply(P1_train[,2:p],2,mean)
  P2_samp_means <- apply(P2_train[,2:p],2,mean)
  Pk_samp_means <- apply(Pk_train[,2:p],2,mean)
  
  X_bar <- data.frame("y" = P1_pop_means- P1_samp_means, 
                      "x2" = P2_samp_means-P1_samp_means,
                      "xk" = Pk_samp_means-P1_samp_means)
  if(case == 3){
    # Calculate the weights
    weights <- coef(lm(y~x2+xk-1, data=X_bar))
    weights[weights > 1] <- 1
    weights[weights < 0] <- 0
    
    # Calculate duc
    df <- data.frame("y" = P1_pop_means- P1_samp_means,
                     "x2" = P2_samp_means- P1_samp_means,
                     "xk" = Pk_samp_means- P1_samp_means)
    
    r1 <- residuals(lm(y~x2-1, data=df))
    r2 <- residuals(lm(xk~x2-1, data=df))
    
    duc <- cor(r1,r2)^2
    
  }else if(case==2){
    # Calculate the weights
    weights <- coef(lm(y~x2-1, data=X_bar))
    weights[weights > 1] <- 1
    weights[weights < 0] <- 0
    duc <- NA
    
  }else{
    weights <- 1
    duc <- NA
    }
  
  res <- list("weights" = weights, "duc" = duc)
  return(res)
}

calc_mse <- function(P1_test, P1_train, P2_train, Pk_train, weights, case=3){
  #' Calculate the mean squared error based on weighted predictions
  #' 
  #' @param P1_test Test samples from target distribution P1
  #' @param P1_train Train samples from target distribution P1
  #' @param P2_train Train samples from P2
  #' @param Pk_train Train samples from Pk
  #' @param weights Weights for P1, P2 and Pk predictions
  #' @param case Case 3 refers to 3 distributions (P1, P2, Pk), 
  #' case 2 refers to 2 distributions (P1, P2), case 1 refers to only P1 samples
  #' 
  #' @return Excess MSE from weighted predictions
  
  p <- ncol(P1_test)
  
  # Test points for evaluation
  X_test <- P1_test[,2:p]
  
  fit1 <- lm(P1_train$Y ~ .-1, data = P1_train[, 2:p])
  fit2 <- lm(P2_train$Y ~ .-1, data = P2_train[, 2:p])
  fitk <- lm(Pk_train$Y ~ .-1, data = Pk_train[, 2:p])
  
  w1 <- max(min(1 - sum(weights), 1), 0)
  # Linear reg predictions and MSE
  if (case == 3) {
    weighted_mean <- w1 * predict(fit1, newdata = X_test) +
      weights[1] * predict(fit2, newdata = X_test) +
      weights[2] * predict(fitk, newdata = X_test)
  }
  if (case == 2) {
    weighted_mean <- w1 * predict(fit1, newdata = X_test) +
      weights[1] * predict(fit2, newdata = X_test)
  }
  if (case == 1) {
    weighted_mean <- predict(fit1, newdata = X_test)
  }
  
  # Predict on test set
  mse <- mean((P1_test$Y - weighted_mean)^2)
  pop_mse <- mean((P1_test$Y - rowSums(X_test))^2) #pop coefs is 1
  
  # Return excess risk
  return(mse-pop_mse)
}

run_experiment <- function(deltas, nk_vec, n_test, trials, num_cores,seed = 123){
  #' Run an experiment across multiple trials in parallel
  #'
  #' @param deltas Numeric vector: shift strength for sources, .. 
  #' @param nk_vec Number of samples for training from P1, P2, ...
  #' @param n_test Number of test samples from P1
  #' @param trials Number of trials to run
  #' @param num_cores Number of cores for parallel processing
  #' @param seed Random seed for reproducibility
  #'
  #' @return A list containing mse_df for different source distributions
  
  set.seed(seed)
  
  # Single trial function
  single_run <- function(trial_id) {
    # Keep trial-specific seed for reproducibility in parallel
    set.seed(seed + trial_id)
    
    # P1
    P1_train <- get_P1(n1_train)
    P1_test <- get_P1(n_test)
    
    # Sample data
    source_list <- sample_data(deltas, nk_vec)
    
    # Samples only from P1 and P2
    weights_and_ducs_P1P2 <- calc_weights_and_duc(
      P1_train, source_list[[1]], source_list[[2]], case=2
      )
    
    # Weights and ducs
    weights_and_ducs <- to_list(for (s in 2:length(deltas)) 
        calc_weights_and_duc(
          P1_train, source_list[[1]], source_list[[s]][1:nk_vec[s],], case=3)
        )
    
    # MSE calculations
    mse_P1 <- calc_mse(P1_test, P1_train, source_list[[1]], source_list[[2]], 
                        1, case=1) 
    mse_P1P2 <- calc_mse(P1_test, P1_train, source_list[[1]], source_list[[2]], 
                         weights_and_ducs_P1P2$weights, case=2) 
    all_mse <- to_vec(for (s in 2:length(deltas)) 
        calc_mse(P1_test, P1_train, source_list[[1]], source_list[[s]][1:nk_vec[s],], 
                 weights_and_ducs[[s-1]]$weights, case=3)
        )
    
    list("mse" = c(mse_P1, mse_P1P2, all_mse), 
         "weights_P2" = weights_and_ducs_P1P2$weights, 
         "weights_Pk" = to_list(for (s in weights_and_ducs) s$weights),
         "ducs" = to_vec(for (s in weights_and_ducs) s$duc))
  }
  
  # Run trials in parallel (reproducible with set.seed)
  results <- parallel::mclapply(
    1:trials, single_run, mc.cores = num_cores
  )
  return(results)
}

trials = 1000
n1_train = 300

# shift strengths taken in by calinf package. See https://arxiv.org/abs/2202.11886
deltas <- c(1.8,   2, 2.1,  2.15, 2.2,  2.5,  2.7,  2.8,  3,    3.2,  3.5,
            4,     4.05, 4.15, 4.17, 4.2,   4.3,   5,     5.1) 
nk_vec <- c(400, 800, 900, 1000, 1300, 1700, 1900, 2000, 2500, 3000, 3200, 5000,
            8000, 9000, 11000, 15000, 19000, 20000, 40000, 50000, 70000) # number of samples

n_test = 50000 # number of test samples
num_cores <- parallel::detectCores() - 1 # Set number of cores for paralleling experiment

# Run experiment in parallel
results <- run_experiment(
  deltas = deltas,
  nk_vec = nk_vec,
  n_test = n_test,
  trials = trials,
  num_cores = num_cores,
)
# saveRDS(results, "results_simulation.RDS")
########

# Avg mse and ducs across trials
all_mse <- lapply(results, `[[`, "mse")
avg_mse <- Reduce("+", all_mse) / length(all_mse)
all_ducs <- lapply(results, `[[`, "ducs")
avg_ducs <- Reduce("+", all_ducs) / length(all_ducs)
avg_rel_imp <- 1-avg_mse[3:length(avg_mse)]/avg_mse[2]

# Avg weights
all_w2 <- lapply(results, `[[`, "weights_P2")
avg_w2 <- Reduce("+", all_w2) / length(all_w2)
all_wk <- lapply(results, `[[`, "weights_Pk")
avg_wk <- lapply(seq_along(all_wk[[1]]), function(k) {
  mats <- lapply(all_wk, `[[`, k)
  Reduce("+", mats) / length(mats)
})

# Print results
avg_ducs
avg_rel_imp
# avg_mse
# avg_w2
# avg_wk

########
# Plotting

df <- data.frame(
  AvgDUC = avg_ducs,
  RelImp = avg_rel_imp
)
max_idx <- which.max(df$AvgDUC)
min_idx <- which.min(df$AvgDUC)

p <- ggplot(df, aes(x = AvgDUC, y = RelImp)) +
  geom_point(size = 3, color = "black") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Predicting model performance without outcome data",
    x = "Avg Estimated DUC",
    y = "1 - Avg Relative Excess MSE "
  ) +
  annotate(
    "text", x = 0.82, y = 0.94, label = "y = x",
    color = "gray40", size = 4, hjust = 0
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) 

print(p) ## 700x400

