################################################
# Recreating simulations from Menard et al. 2014

# Generating X1 and X2 using either MVN or
#   MVN with exp(X2)

# Cor X1 and X2 seq(0,0.9,0.1)

# Outcome y either linear regression or logistic
# with varying betas for logistic and varying 
# sigma_y's for linear

# Continuous, 2, 3, 4, 5 categories for X2

# equal sized groups or quantiles

# For categorization: hard, 
#   triangle soft (w/wout extra categories), 
#   normal soft (w/wout extra categories)

# Sample size 100, 1000

# 500 iterations for each to get coverage
# of 95% confidence intervals
# for beta_1 of 0 (Type I error)

################################################
# libraries

library(dplyr)
library(MASS)
library(EnvStats)


################################################
# function to generate X

generate_x <- function(n, sigma_x, method = "normal"){
  mu_x <- c(0,0)
  Sigma_x <- matrix(c(1, sigma_x, sigma_x, 1), byrow = T, nrow = 2)
  
  X <- mvrnorm(n, mu_x, Sigma_x)
  
  if(method == "normal"){
    return(X)
  } else{
    X[,2] <- exp(X[,2])
    return(X)
  }
}

#X <- generate_x(n, sigma_x, method = "normal")


################################################
# function to generate outcome y


generate_y <- function(X, y_param, method = "linear"){
  if(method == "linear"){
    y <- X[,2] + rnorm(nrow(X), sd = y_param)
  } else{
    p_y <- 1/(1+exp(-(y_param*X[,2])))
    y <- rbinom(nrow(X), 1, p_y)
  }
  return(y)
}

#y <- generate_y(X, beta_lin, sigma_y, method = "linear")
#y <- generate_y(X, beta_log, sigma_y, method = "logistic")


################################################
# function to categorize X2


categorize_X2 <- function(X_2, categories=1, alg = "equal", method = "hard"){
  
  if(categories == 1){
    return(X_2)
  }
  
  
  # equal sized categories or categories by quantile
  if(alg == "equal"){
    cuts_X <- seq(min(X_2), max(X_2), length = categories+1)
    groups_X <- as.numeric(cut(X_2, cuts_X, include.lowest = T))
  } else{
    cuts_X <- quantile(X_2, probs = seq(0, 1, length = categories + 1))
    groups_X <- as.numeric(cut(X_2, cuts_X, include.lowest = T))
  }
  
  if(method == "hard"){
    X2_cols <- model.matrix( ~ as.factor(groups_X) - 1)
  } else if(method == "triangle_add"){
    
    # if unequal groups, still take average group size as triangle width
    g_size = mean(diff(cuts_X))
    
    weights <- matrix(nrow = length(X_2), ncol = categories + 2)
    
    weights[,1] <- (ptri(cuts_X[1], min = X_2 - g_size, max = X_2 + g_size, mode = X_2) -
                      ptri(cuts_X[1] - g_size, min = X_2 - g_size, max = X_2 + g_size, mode = X_2))
    
    for(i in 2:length(cuts_X)){
      weights[,i] <- (ptri(cuts_X[i], min = X_2 - g_size, max = X_2 + g_size, mode = X_2)-
                          ptri(cuts_X[i-1], min = X_2 - g_size, max = X_2 + g_size, mode = X_2))
    }
    
    weights[,length(cuts_X) + 1] <- (ptri(cuts_X[length(cuts_X)] + g_size, 
                                          min = X_2 - g_size, max = X_2 + g_size, mode = X_2)-
                                       ptri(cuts_X[length(cuts_X)], 
                                            min = X_2 - g_size, max = X_2 + g_size, mode = X_2))
    
    # normalize weights to 1
    weights <- weights/rowSums(weights)
    
    
    # group means, add additional groups above and below
    group_means <- c(mean(c(cuts_X[1], cuts_X[1] - g_size)),
                            sapply(1:(length(cuts_X) - 1), function(i) mean(c(cuts_X[i], cuts_X[i + 1]))),
                     mean(c(cuts_X[length(cuts_X)], cuts_X[length(cuts_X)] + g_size)))
    
    # multiply weights by group means to keep coefficients the same scale
    X_2_cat <- t(t(weights)*group_means)
    
    return(X_2_cat)
    
  } else if(method == "triangle_no_add"){
    # if unequal groups, still take average group size as triangle width
    g_size = mean(diff(cuts_X))
    
    weights <- matrix(nrow = length(X_2), ncol = categories)
    
    for(i in 1:(length(cuts_X) - 1)){
      weights[,i] <- (ptri(cuts_X[i+1], min = X_2 - g_size, max = X_2 + g_size, mode = X_2)-
                        ptri(cuts_X[i], min = X_2 - g_size, max = X_2 + g_size, mode = X_2))
    }
    
    # normalize weights to 1
    weights <- weights/rowSums(weights)
    
    # group means
    group_means <- c(sapply(1:(length(cuts_X) - 1), function(i) mean(c(cuts_X[i], cuts_X[i + 1]))))
    
    # multiply weights by group means to keep coefficients the same scale
    X_2_cat <- t(t(weights)*group_means)
    
    return(X_2_cat)
  } else if(method == "normal_add"){
    
    # if unequal groups, still take average group size as triangle width
    g_size = mean(diff(cuts_X))
    
    weights <- matrix(nrow = length(X_2), ncol = categories + 2)
    
    # take sd of normal to be group size/2
    weights[,1] <- (pnorm(cuts_X[1], mean = X_2, sd = g_size/2) -
                      pnorm(cuts_X[1] - g_size, mean = X_2, sd = g_size/2))
    
    for(i in 2:length(cuts_X)){
      weights[,i] <- (pnorm(cuts_X[i], mean = X_2, sd = g_size/2)-
                        pnorm(cuts_X[i-1], mean = X_2, sd = g_size/2))
    }
    
    weights[,length(cuts_X) + 1] <- (pnorm(cuts_X[length(cuts_X)] + g_size, mean = X_2, sd = g_size/2)-
                                       pnorm(cuts_X[length(cuts_X)], mean = X_2, sd = g_size/2))
    
    # normalize weights to 1
    weights <- weights/rowSums(weights)
    
    
    # group means, add additional groups above and below
    group_means <- c(mean(c(cuts_X[1], cuts_X[1] - g_size)),
                     sapply(1:(length(cuts_X) - 1), function(i) mean(c(cuts_X[i], cuts_X[i + 1]))),
                     mean(c(cuts_X[length(cuts_X)], cuts_X[length(cuts_X)] + g_size)))
    
    # multiply weights by group means to keep coefficients the same scale
    X_2_cat <- t(t(weights)*group_means)
    
    return(X_2_cat)
    
  } else{
    # if unequal groups, still take average group size as triangle width
    g_size = mean(diff(cuts_X))
    
    weights <- matrix(nrow = length(X_2), ncol = categories)
    
    for(i in 1:(length(cuts_X) - 1)){
      weights[,i] <- (pnorm(cuts_X[i+1], mean = X_2, sd = g_size/2)-
                        pnorm(cuts_X[i], mean = X_2, sd = g_size/2))
    }
    
    # normalize weights to 1
    weights <- weights/rowSums(weights)
    
    # group means
    group_means <- c(sapply(1:(length(cuts_X) - 1), function(i) mean(c(cuts_X[i], cuts_X[i + 1]))))
    
    # multiply weights by group means to keep coefficients the same scale
    X_2_cat <- t(t(weights)*group_means)
    
    return(X_2_cat)
  } 
  
}


#X_2_cat <- categorize_X2(X[,2], categories = 5, method = "triangle_add")
#X_2_cat <- categorize_X2(X[,2], categories = 1)

################################################
# sim parameters

# categorical sim parameters
X_gen <- c("normal", "lognormal")
sigma_x <- seq(0,0.9,0.1)
y_gen <- c("linear", "logistic")
# 0.5 and 2 are for the logistic regression, 3.17 and 1.02 the linear regression
y_param <- c(0.5, 2, 3.17, 1.02)
categories <- c(2,3,4)
cat_alg <- c("equal", "quantile")
cat_method <- c("hard", "triangle_add", "triangle_no_add", "normal_add", "normal_no_add")
n <- c(100, 1000)

# sim params
sim_params <- expand.grid(X_gen = X_gen, sigma_x = sigma_x, y_gen = y_gen,
                          y_param = y_param, categories = categories,
                          cat_alg = cat_alg, cat_method = cat_method, n = n)

# remove parameters corresponding to other y generation model
sim_params <- sim_params %>% filter(!((y_gen == "linear" & y_param %in% c(0.5, 2)) |
                                      (y_gen == "logistic" & y_param %in% c(3.17, 1.02))))

# add parameters for keeping X_2 continuous
sim_params <- sim_params %>% bind_rows(expand.grid(X_gen = X_gen, sigma_x = sigma_x, y_gen = y_gen,
                                     y_param = y_param, categories = 1, n = n))

# add column for results
sim_params <- sim_params %>% mutate(type_1_error = NA)


################################################
# run simulation
# repeat simulation 500 (for now) times for each setting
# record 95% credible interval and see if contains 0 for beta_1
# average over 500 sims for type_1_error in each case

# num mc iter
num_iter <- 500

for(i in 1:nrow(sim_params)){
  
  print(i)
  
  sim_mc <- data.frame(beta_1_est = rep(NA, num_iter),
                       beta_1_025 = rep(NA, num_iter), beta_1_975 = rep(NA,num_iter))
  
  X_gen_sim <- sim_params[i,]$X_gen
  sigma_x_sim <- sim_params[i,]$sigma_x
  y_gen_sim <- sim_params[i,]$y_gen
  y_param_sim <- sim_params[i,]$y_param
  categories_sim <- sim_params[i,]$categories
  cat_alg_sim <- sim_params[i,]$cat_alg
  cat_method_sim <- sim_params[i,]$cat_method
  n_sim <- sim_params[i,]$n
  
  
  for(j in 1:num_iter){
    
    X <- generate_x(n_sim, sigma_x_sim, method = X_gen_sim)
    y <- generate_y(X, y_param = y_param_sim, method = y_gen_sim)
    
    # if continuous
    if(categories_sim == 1){
      if(y_gen_sim == "linear"){
        sum_model <- summary(lm(y ~ X))$coef
      } else{
        sum_model <- summary(glm(y ~ X, family = binomial(link = "logit")))$coef
      }
      
      sim_mc[j,] <- c(sum_model[2,1], sum_model[2,1] - qnorm(0.975)*sum_model[2,2],
                      sum_model[2,1] + qnorm(0.975)*sum_model[2,2])

    } else{
      
      X_2_cat <- categorize_X2(X[,2], categories = categories_sim,
                               alg = cat_alg_sim, method = cat_method_sim)
      
      if(y_gen_sim == "linear"){
        sum_model <- summary(lm(y ~ X[,1] + X_2_cat - 1))$coef
      } else{
        sum_model <- summary(glm(y ~ X[,1] + X_2_cat - 1, family = binomial(link = "logit")))$coef
      }
      
    }
    
    sim_mc[j,] <- c(sum_model[1,1], sum_model[1,1] - qnorm(0.975)*sum_model[1,2],
                    sum_model[1,1] + qnorm(0.975)*sum_model[1,2])
    
  }
  
  sim_params[i,]$type_1_error <- 1 - mean(sim_mc$beta_1_025 < 0 & sim_mc$beta_1_975 > 0)
  
}


# save sim results
write.csv(sim_params, file = "full-sim-results.csv")




