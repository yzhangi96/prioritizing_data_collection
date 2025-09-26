## Import libraries
library(calinf) # for simulating random shift. See https://arxiv.org/abs/2202.11886
library(comprehenr) # for list comprehension
library(ggplot2) 
library(dplyr)
library(parallel)
library(caret)       
library(stats)      
library(ks)        
library(arrow)
library(patchwork)
library(MASS)       
EPS <- 1e-6

get_P1 <- function(n, seed, d=15){
  #' Get samples from target distribution P1
  #' @param n Number of samples to generate
  #' @param seed Random seed for reproducibility
  #' @param d 1/2 Number of covariates 
  #' 
  #' @return Dataframe containing samples from target distribution P1
  
  set.seed(seed)
  
  # Generate Normal(0,1)
  X_norm <- replicate(d, rnorm(n, mean=0, sd=1))
  
  # Apply nonlinear transforms
  X_norm[,1] <- X_norm[,1]^2   # X1^2
  X_norm[,2] <- X_norm[,2]^2   # X2^2
  
  # Generate Bern(0.5)
  X_bin <- replicate(d, rbinom(n, size=1, prob=0.5))
  noise <- runif(n, -1, 1)

  X <- cbind(X_bin, X_norm)
  Y <- rowSums(X) + noise
  df <- data.frame(Y = Y, X)
  return(df)
}

get_random_shifts <- function(delta,n, seed, d=15){
  #' Get randomly shifted X,y for a given delta
  #' @param delta Strength of shift
  #' @param n Number of samples to generate
  #' @param seed Random seed for reproducibility
  #' @param d 1/2 Number of covariates 
  #' 
  #' @return Dataframe containing samples from a randomly shifted distribution
  
  set.seed(seed)
  dseed <- distributional_seed(n=n,delta=delta)
  
  # Generate X1–Xd ~ Normal(0,1)
  X_norm <- replicate(d, drnorm(dseed, mean=0, sd=1))
  
  # Apply nonlinear transforms
  X_norm[,1] <- X_norm[,1]^2   # X1^2
  X_norm[,2] <- X_norm[,2]^2   # X2^2
  
  # Generate X_d+1–X_2d ~ Bern(0.5)
  X_bin <- replicate(d, drbinom(dseed, size=1, prob=0.5))
  noise <- drunif(dseed, -1, 1)
  
  X <- cbind(X_bin, X_norm)
  Y <- rowSums(X) + noise
  df <- data.frame(Y = Y, X)
  return(df)
}

sample_data <- function(deltas, nk_vec, seed) {
  #' Sample data from distributions under random shifts
  #' 
  #' @param deltas A vector of distr_var values to indicate 
  #' strength of shifts 
  #' @param nk_vec A vector of number of samples to draw from each distribution
  #' @param seed Random seed for reproducibility
  #' 
  #' @return list of dataframes, each containing samples from a distribution.
  
  ## Sample Pk
  
  source_list <- list()
  for(i in 1:length(deltas)){
    source_list[[i]] <- get_random_shifts(deltas[i], nk_vec[i], seed)
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
  P1_pop_means <- c(rep(0.5, (p-1)/2),c(1,1), rep(0, (p-1)/2-2)) #population means

  # Get means of covariates for each distribution
  P1_samp_means  <- apply(P1_train[,2:p],2,mean)
  P2_samp_means <- apply(P2_train[,2:p],2,mean)
  Pk_samp_means <- apply(Pk_train[,2:p],2,mean)
  
  X_bar <- cbind(
    x2 = P2_samp_means - P1_samp_means,
    xk = Pk_samp_means - P1_samp_means
  )
  y_bar <- P1_pop_means - P1_samp_means
  
  if(case == 3){
    # Calculate the weights
    
    ### Using lm is much slower
    XtX <- crossprod(X_bar)
    Xty <- crossprod(X_bar, y_bar)
    weights <- as.numeric(solve(XtX, Xty))
    weights <- pmin(pmax(weights, 0), 1)
    
    # Calculate duc
    r1 <- y_bar - X_bar[,1] * (sum(X_bar[,1] * y_bar) / sum(X_bar[,1]^2))
    r2 <- X_bar[,2] - X_bar[,1] * (sum(X_bar[,1] * X_bar[,2]) / sum(X_bar[,1]^2))
    
    duc <- cor(r1,r2)^2
    
  }else if(case==2){
    # Calculate the weights
    weights <- sum((X_bar[,1]*y_bar)) / sum(X_bar[,1]^2)
    weights <- pmin(pmax(weights, 0), 1)
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
  
  ## Test points for evaluation
  X_test <- as.matrix(P1_test[, 2:p])
  
  ## Using lm function much slower
  get_betas <- function(df) {
    X <- as.matrix(df[, 2:p])
    y <- df$Y
    as.numeric(solve(crossprod(X), crossprod(X, y)))
  }
  
  b1 <- get_betas(P1_train)
  b2 <- get_betas(P2_train)
  bk <- get_betas(Pk_train)
  
  w1 <- max(min(1 - sum(weights), 1), 0)
  
  # Linear reg predictions and MSE
  if (case == 3) {
    weighted_mean <- w1*(X_test %*% b1) + weights[1]*(X_test %*% b2) + weights[2]*(X_test %*% bk)
  } else if (case == 2) {
    weighted_mean <- w1*(X_test %*% b1) + weights[1]*(X_test %*% b2)
  } else {
    weighted_mean <- X_test %*% b1
  }
  
  # Predict on test set
  mse <- mean((P1_test$Y - weighted_mean)^2)
  pop_mse <- mean((P1_test$Y - rowSums(X_test))^2) #pop coefs is 1
  
  # Return excess risk
  return(mse-pop_mse)
}

compute_pca <- function(P1_train, P2_train, Pk_train, seed) {
  #' Compute PCA projection for input of KL calculations as in 
  #' Shen et. al (2024)
  #' 
  #' @param P1_train Subset of samples from target distribution P1
  #' @param P2_train Subset of samples from P2
  #' @param Pk_train Subset of samples from Pk
  #' @param seed Random seed for reproducibility
  #' 
  #' @return A list containing PCA projections of X1, X2 and Xk

  set.seed(seed)
  p <- ncol(P1_train)
  X1 <- P1_train[,2:p]
  X2 <- P2_train[,2:p]
  Xk <- Pk_train[,2:p]
  
  combined <- rbind(X1, X2, Xk)
  
  # Scale and apply PCA to 3 components
  combined_scaled <- scale(combined)
  pca <- prcomp(combined_scaled, center = FALSE, scale. = FALSE)
  
  # Project X1X2 and Xk
  X1X2_proj <- predict(pca, scale(rbind(X1, X2), center = attr(combined_scaled, "scaled:center"),
                                scale  = attr(combined_scaled, "scaled:scale")))[, 1:3]
  Xk_proj <- predict(pca, scale(Xk, center = attr(combined_scaled, "scaled:center"),
                                scale  = attr(combined_scaled, "scaled:scale")))[, 1:3]
  
  return(list(X1X2_proj = X1X2_proj, Xk_proj = Xk_proj))
}

fit_density <- function(X, seed, bandwidths = 10^seq(-1, 1, length.out = 10)) {
  #' Fit a kernel density estimator with CV for 
  #' for bandwidth selection as in Shen et. al (2024)
  #' 
  #' @param X Data matrix for which to fit the density
  #' @param seed Random seed for reproducibility
  #' @param bandwidths A vector of bandwidths to consider for CV. 
  #' Default to that used in Shen et. al (2024)
  #' 
  #' @return Fitted kde object from ks package
  
  set.seed(seed)
  
  k <- 5
  n <- nrow(X)
  d <- ncol(X)
  folds <- sample(rep(1:k, length.out = n))
  
  best_bw <- NULL
  best_ll <- -Inf
  
  for (bw in bandwidths) {
    fold_lls <- numeric(k)
    
    # Pre-build H matrix once
    H <- diag(bw, d)
    
    for (fold in 1:k) {
      idx_train <- folds != fold
      idx_valid <- !idx_train
      
      X_train <- X[idx_train, , drop = FALSE]
      X_valid <- X[idx_valid, , drop = FALSE]
      
      # Directly evaluate density only at validation points
      kde_fit <- ks::kde(X_train, H = H, eval.points = X_valid, binned = TRUE)
      dens_valid <- kde_fit$estimate + EPS
      
      fold_lls[fold] <- mean(log(dens_valid))
    }
    
    mean_ll <- mean(fold_lls)
    if (mean_ll > best_ll) {
      best_ll <- mean_ll
      best_bw <- bw
    }
  }
  
  # Fit final KDE on all data using best bandwidth
  return(ks::kde(X, H = diag(best_bw, d), binned = TRUE))
}
  

compute_kl_and_score <- function(P1_train, P2_train, Pk_train, seed) {
  #' Compute KL divergence and score fn between densities estimated from
  #' (X1, X2) and Xk after PCA projection
  #' @param P1_train Subset of samples from target distribution P1
  #' @param P2_train Subset of samples from P2
  #' @param Pk_train Subset of samples from Pk
  #' @param seed Random seed for reproducibility
  #' 
  #' @return list of KL divergence value and score_x
  
  pcs <- compute_pca(P1_train, P2_train, Pk_train, seed)
  X1X2_pca <- pcs$X1X2_proj
  Xk_pca <- pcs$Xk_proj
  
  kde_t <- fit_density(X1X2_pca, seed)
  kde_s <- fit_density(Xk_pca, seed)
  
  # Evaluate densities 
  t <- predict(kde_t, x = X1X2_pca) + EPS
  s <- predict(kde_s, x = X1X2_pca) + EPS
  
  t <- t / sum(t)
  s <- s / sum(s)
  
  # KL(t||s)
  kl <- sum(t * log(t / s))
  
  log_t <- log(predict(kde_t, x = X1X2_pca) + EPS)
  log_s <- log(predict(kde_s, x = X1X2_pca) + EPS)
  
  log_ratio <- log_t - log_s
  score_x <- mean(log_ratio) / sd(log_ratio)
  
  return(list("kl" = kl, "score_x" = score_x))
}


run_experiment <- function(deltas, nk_vec, n_test, trials, num_cores, trial,
                           compute_kl=FALSE){
  #' Run an experiment across multiple trials in parallel
  #'
  #' @param deltas Numeric vector: shift strength for sources, .. 
  #' @param nk_vec Number of samples for training from P1, P2, ...
  #' @param n_test Number of test samples from P1
  #' @param trials Number of trials to run
  #' @param num_cores Number of cores for parallel processing
  #' @param trial Use as random seed for reproducibility
  #' @param compute_kl Whether to compute KL and score_x in R (default FALSE).
  #' If False, stores samples as parquets 
  #'
  #' @return A list containing mse_df for different source distributions
  
  
  # Single trial function
  single_run <- function(trial_id) {

    # P1
    P1_train <- get_P1(n1_train, seed = trial_id)
    P1_test <- get_P1(n_test, seed = trial_id)
    
    # Sample data
    source_list <- sample_data(deltas, nk_vec, seed = trial_id)
    
    # Samples only from P1 and P2
    weights_and_ducs_P1P2 <- calc_weights_and_duc(
      P1_train, source_list[[1]], source_list[[2]], case=2
      )

    # Weights and ducs (~0.05 secs)
    weights_and_ducs <- to_list(for (s in 2:length(deltas)) 
        calc_weights_and_duc(
          P1_train, source_list[[1]], source_list[[s]][1:nk_vec[s],], case=3)
        )

    if(compute_kl){
      # Calculate KL and score
      kl_score_x <- to_list(for (s in 2:length(deltas))
        compute_kl_and_score(
          P1_train, source_list[[1]], source_list[[s]][1:nk_vec[s],], seed=trial_id)
      )
      kl_x <- to_vec(for (s in kl_score_x) s$kl)
      score_x <- to_vec(for (s in kl_score_x) s$score_x)
    } else{
      kl_x <- NA
      score_x <- NA
      # Save datasets inside "data" folder
      write_parquet(P1_train, file.path("data", paste0("P1_train_trial", trial_id, ".parquet")))
      write_parquet(P1_test,  file.path("data", paste0("P1_test_trial",  trial_id, ".parquet")))
      
      # Save each source in source_list
      for (i in seq_along(source_list)) {
        write_parquet(
          source_list[[i]],
          file.path("data", paste0("source_", i, "_trial", trial_id, ".parquet"))
        )
      }
    }
    
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
         "ducs" = to_vec(for (s in weights_and_ducs) s$duc),
         "kl_x" = kl_x, "score_x" = score_x
         )
  }
  
  # Run trials in parallel (reproducible with set.seed)
  results <- parallel::mclapply(
    1:trials, single_run, mc.cores = num_cores
  )
  return(results)
}

## CAUTION: Takes a long time to run KL/score calculations 
## because of kde function in R. When compute_kl=FALSE in run_experiment
## (the defalse), it will store a in a file named "data" the samples.
## The graph below reads in KL/score calculations done in python script called
## kl_score_x.ipynb using samples stored from experiments from R
## Note that it still takes ~1 min per trial to run the KL/Score calculations
## on 24GB Apple M3 

trials = 1000
n1_train = 300

# shift strengths taken in by calinf package. See https://arxiv.org/abs/2202.11886
deltas <- c(1.8,    1.5,  2, 2.1,  2.15, 2.2,  2.5,  2.6,  2.8,  3,    3.2,  3.5,
            4,     4.10,   4.4, 4.42) 
nk_vec <- c(400, 500, 900, 1000, 1200, 1400, 1800, 1950, 2200, 2500, 3200, 4000, 7000,
            9000,  18000, 21000) # number of samples

n_test = 50000 # number of test samples
num_cores <- 4


# Run experiment in parallel
results <- run_experiment(
  deltas = deltas,
  nk_vec = nk_vec,
  n_test = n_test,
  trials = trials,
  num_cores = num_cores,
)

saveRDS(results, "results_simulation_comparison.RDS")

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

########
# Plotting

avg_results_kl_score <- read.csv("data/kl_score_x_avg_v3.csv")
df<- data.frame(
  AvgDUC = avg_ducs,
  NegRelImp = -1*avg_rel_imp,
  AvgNegKL = -1*avg_results_kl_score$kl,
  AvgNegScoreX = -1*avg_results_kl_score$score_x
)

## Calculate corr and linear fit for plots
coefs_kl <- coef(lm(avg_rel_imp ~ df$AvgNegKL))
corr_kl <- cor(df$NegRelImp, df$AvgNegKL)
corr_duc <- cor(df$NegRelImp, df$AvgDUC)
coefs_score <- coef(lm(avg_rel_imp ~ df$AvgNegScoreX))
corr_score <- cor(df$NegRelImp, df$AvgNegScoreX)


p <- ggplot(df, aes(x = AvgDUC, y = NegRelImp)) +
  geom_point(size = 3, color = "black") +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Avg Estimated DUC",
    y = "Avg Neg Relative Excess MSE "
  ) +
  annotate(
    "text", x = 0.72, y = 0, label = "y = -x",
    color = "black", size = 3.5, hjust = 0
  ) +
  annotate(
    "text", x = 0.72, y = -0.1, label = "Corr = -0.99",
    color = "black", size = 3.5, hjust = 0
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, by = 0.1)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) 

p2 <- ggplot(df, aes(x = AvgNegKL, y = NegRelImp)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
  labs(
    x = "Avg Neg KL",
  ) +
  annotate(
    "text", x = -0.20, y = 0, label = "y=0.61+0.92x",
    color = "black", size = 3.5, hjust = 0
  ) +
  annotate(
    "text", x = -0.20, y = -0.1, label = "Corr = -0.76",
    color = "black", size = 3.5, hjust = 0
  ) +
  scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, by = 0.1)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )


p3 <- ggplot(df, aes(x = AvgNegScoreX, y = NegRelImp)) +
  geom_point(size = 3, color = "black") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
  labs(
    x = "Avg Neg Score",
  ) +
  annotate(
    "text", x = -0.65, y = 0, label = "y=1+0.74x",
    color = "black", size = 3.5, hjust = 0
  ) +
  annotate(
    "text", x = -0.65, y = -0.1, label = "Corr = -0.90",
    color = "black", size = 3.5, hjust = 0
  ) +
  scale_y_continuous(limits = c(-1, 0), breaks = seq(-1, 0, by = 0.1)) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p2 <- p2 +
  theme(
    axis.title.y = element_blank(),
  )

p3 <- p3 +
  theme(
    axis.title.y = element_blank(),
  )

combined <- (p + p2 + p3) +
  plot_layout(ncol = 3) +
  plot_annotation(
    title = "Predicting model performance without outcome data"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5) 
  )

print(combined)

