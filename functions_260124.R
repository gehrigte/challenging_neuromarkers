# Challenging the Search for Neuromarkers of Mental Disorders
# Technical Script/ Part 1

## In this script, we pre-define the simulation and plot functions, that are 
## used (i.e., called) in the main script. 
## Any mathematical or technical questions (i.e., how we operationalized a certain
## idea or construct) can be (attempted to) answer here.

#---------------------
# LIBRARIES
#---------------------
# We load libraries that are necessary for following the next steps of this script.
# If a library is not installed, you can easily install it by calling 
# install.packages('foo'), and replace 'foo' with the respective package name.

library(faux)
library(dplyr)
library(corrplot)
library(corpcor)
library(Matrix)
library(highfrequency)
library(rockchalk)
library(caret)
library(caTools)
library(ggplot2)
library(viridis)
library(broom)
library(tidyverse)

#---------------------
# SIMULATION FUNCTIONS
#---------------------

# To define the following functions, we need to set 
# N = Number of observations per individual
# cor_bs = correlation between brain region and symptom
# cor_bb = mean correlation between brain regions
# when calling the function.

# HIGH MULTIPLE REALIZABILITY (MR), i.e., activity in random brain regions can 
# elicit the same symptom.
# Our simulation is based on 8 brain regions, meaning we pre-define a vector
# of brain-symptom correlation with length 8. This can be easily adapted.

make_high_mr <- function(N, cor_bs, cor_bb, var_names){ 
  p <- length(cor_bs)
  # sample from the brain-symptom correlation vector to create different correlation patterns 
  s <- sample(cor_bs, p, replace = FALSE)
  # make a positive definite correlation matrix from the vector of 
  # correlations with the symptom variable s:
  # first create a PD matrix:
  m <- matrix(cor_bb, ncol = p, nrow = p)
  diag(m) <- 1
  s1 <- c(1, s)
  m_new <- cbind(s, m)
  m <- rbind(s1, m_new)
  colnames(m) <- rownames(m) <-  var_names
  
  # we create the data frame (T:REPRESENTING ACTIVATION PATTERNS?)
  # we draw an n-dimensional vector of observations using the PD correlation matrix 
  # to specify the relationships between time points
  data <- rnorm_multi(n = N, vars = p+1, r = m)
  
  # usually the numbers are not exactly as intended, but we can decide what 
  # deviations are acceptable
  # here we make sure the first row is not too different from the intended numbers 
  # (var 0.01, mean 0.2):
  while(((abs(mean(cor(data)[2:9,1]) - mean(cor_bs))) > 0.005) | 
        (abs(var(cor(data)[2:9,1]) - var(cor_bs)) > 0.004)) {
    data <- rnorm_multi(n = N, vars = p + 1, r = m)
  }
  
  return(data)
}


# LOW MR, i.e., similar activity in brain regions across individuals elicits
# the same symptom.
make_low_mr <- function(N, cor_bs, cor_bb, var_names) {
  # function for low MR
  p <- length(cor_bs)
  #create a positive definite correlation matrix
  # this is the vector of the correlations with the symptom variable
  # note: we do not sample here, so the order of brain-symptom corr is fixed
  s <- cor_bs
  #we first create a PD matrix
  m <- matrix(cor_bb, ncol = p+1, nrow = p+1)
  diag(m) <- 1
  s1 <- c(1, s)
  m[ ,1] <- m[1, ] <- s1
  colnames(m) <- rownames(m) <-  var_names
  
  #we create the data frame with the PD matrix
  data <- rnorm_multi(n = N, vars = p + 1, r = m)
  
  #usually the numbers are not exactly as intended, but we can decide what deviations are acceptable
  #here i make sure the first row is not too different from the intended numbers (var 0.01, mean 0.2)
  
  while(((abs(mean(cor(data)[2:9,1])-mean(cor_bs))) > 0.005) | 
        (abs((var(cor(data)[2:9,1])-var(cor_bs))) > 0.004) | 
        is.unsorted(cor(data)[2:9,1])) {
    data <- rnorm_multi(n = N, vars = p + 1, r = m)
  }
  return(data)
}


make_medium_mr <- function(N, cor_bs, cor_bb, var_names) {
  p <- length(cor_bs)
  s_neg <- sample(cor_bs[cor_bs < 0], replace = F)
  s_pos <- sample(cor_bs[cor_bs > 0], replace = F)
  s <- c(s_neg, s_pos)
  
  m <- matrix(cor_bb, ncol = p+1, nrow = p+1)
  diag(m) <- 1
  s1 <- c(1, s)
  m[1, ] <- m[, 1] <- s1
  # mat <- makePsd(m, method = "correlation")
  # mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
  
  #we create the data frame with the PD matrix 
  data <- rnorm_multi(n = N, vars = p+1, r = m, varnames = var_names)
  
  #usually the numbers are not exactly as intended, but we can decide what deviations are acceptable
  #here i make sure the first row is not too different from the intended numbers (var 0.01, mean 0.2)
  while(((abs(mean(cor(data)[2:9,1]) - mean(cor_bs))) > 0.005) | 
        (abs((var(cor(data)[2:9,1]) - var(cor_bs))) > 0.004) |
        any((cor(data)[1,2:5])>= 0) | 
        any((cor(data)[1,6:9])<= 0)) {
    data <-  rnorm_multi(n = N, vars = p +1, r = m, varnames = var_names)
  }
}




#-------------------------
#MODEL DIAGNOSTICS
#-------------------------

# The function 'make_mod_diagnostics' (below) calculates a linear regression of the brain 
# regions on the symptom variable S. From the output, R squared is calculated. 
# Next, an idicator variable D, signifying a diagnosis if the symptom variable is 
# above its median, is added to the data. This is used to train a predictive model, 
# in which activation patterns of brain regions (B1+...+B8) are used to predict 
# a diagnosis (D). The classification accuracy (CA) of the model is then calculated.
# The function returns a list of information, containing: 
# the linear models per individual (foo$lm_individ) [more detailed model info 
#                             can be accessed via summary(foo$lm_individ[[i]])],
# R squared per individual (foo$rsq_individ), 
# the mean R squared per individual (foo$mean_rsq_individ),
# predictive models per individual (foo$pm_individ),
# the classification accuracy (foo$ca_individ), and
# and mean classification accuracy per individual (foo$mean_ca_individ).



make_diagnostics <- function(data, individual = TRUE) {
  if (individual == TRUE) {
    models <- data %>% 
      group_by(ID) %>% 
      do(lin_model = lm(S ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8, data = .))
    # https://stackoverflow.com/questions/1169539/linear-regression-and-group-by-in-r
    
    # R-squared for predictive model of brain region activity on symptom 
    # on individual level:
    r_sq_id <- rep(NA, N_ind)
    for (i in 1: N_ind) {
      r_sq_id[i] <- summary(models$lin_model[[i]])$r.squared
    }
    
    # To calculate classification accuracy (CA) we train a predictive model. 
    # The model aims to predict the diagnosis (D) using the brain activity.
    # To build the model, we use 80% of the full data set as training data 
    # and test it on a test subset (20%) of the data.
    
    # Create training and test data subsets:
    data_train <- data_full %>% 
      group_by(ID) %>% 
      slice_sample(prop = .8)
    
    data_test <- data_full %>% 
      group_by(ID) %>% 
      do(subset(., !(.$rownum%in%data_train$rownum))) # select remaining observations that are NOT in training data
    
    # Train model:
    trained_models <- data_train %>% 
      do(prediction_mods = train(form = D ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8,
                                data = .,
                                trControl = trainControl(method = "cv", number = 5),
                                method = "glm",
                                family = "binomial")
      )
    
    # # Check whether function above yields same results as step-wise procedure:
    # 
    # data_train_ind <- data_train[data_train$ID == 1,]
    # data_test_ind <- data_test[data_test$ID == 1,]
    # 
    # test_model <- train(form = D ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8,
    #       data = data_train_ind,
    #       trControl = trainControl(method = "cv", number = 5),
    #       method = "glm",
    #       family = "binomial")
    # out_test <- predict(test_model, data_test_ind)
    # mean(out_test == data_test_ind$D)
    
    predicted_data <- list(rep(NA, N_ind))
    ca_id <- list(rep(NA, N_ind))
    
    # Use trained model to predict a diagnosis
    for (i in 1: N_ind) {
      predicted_data[[i]] <- predict(trained_models$prediction_mods[[i]], 
                                     newdata = data_test[data_test$ID ==i,])
      # Classification Accuracy = mean(predicted diagnosis == true diagnosis):
      ca_id[[i]] <- mean(unlist(predicted_data[[i]]) == 
                           unlist(data_test[data_test$ID == i, 'D']))  
    }
    out <- list('lm_individ' = models$lin_model,
                'rsq_indiv' = r_sq_id,
                'mean_rsq_indiv' = mean(r_sq_id),
                'pm_individ' = models,
                'ca_individ' = ca_id,
                'mean_ca_individ' = mean(unlist(ca_id)))
    
    return(out)
  }
  else if (individual == FALSE) {
    
    models_m <- data %>% 
      group_by(t) %>% # group data by time points
      do(lin_model = lm(S ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8, data = .))
    # https://stackoverflow.com/questions/1169539/linear-regression-and-group-by-in-r
    
    # R-squared for predictive model of brain region activity on symptom 
    # per time point:
    r_sq_msmt <- rep(NA, N)
    for (m in 1: N) {
      r_sq_msmt[m] <- summary(models_m$lin_model[[m]])$r.squared
    }
    
    data_train <- data %>% 
      group_by(t) %>% 
      slice_sample(prop = .8)
    
    data_test <- data %>% 
      group_by(t) %>% 
      do(subset(., !(.$rownum%in%data_train$rownum))) # select remaining observations that are NOT in training data
    
    # Train model:
    trained_models_m <- data_train %>% 
      do(prediction_mods = train(form = D ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8,
                                data = .,
                                trControl = trainControl(method = "cv", number = 5),
                                method = "glm",
                                family = "binomial")
      )
  
    
    predicted_data <- list(rep(NA, N))
    ca_msmt <- list(rep(NA, N))
    
    # Use trained model to predict a diagnosis:
    for (m in 1: N) {
      predicted_data[[m]] <- predict(trained_models_m$prediction_mods[[m]], 
                                     newdata = data_test[data_test$t ==m,])
      # CA = mean(predicted diagnosis == true diagnosis):
      ca_msmt[[m]] <- mean(unlist(predicted_data[[m]]) == 
                           unlist(data_test[data_test$t == m, 'D']))  
    }
    out <- list('lm_msmt' = models_m$lin_model,
                'rsq_msmt' = r_sq_msmt,
                'mean_rsq_msmt' = mean(r_sq_msmt),
                'pm_msmt' = models_m,
                'ca_msmt' = ca_msmt,
                'mean_ca_msmt' = mean(unlist(ca_msmt)))
    
    return(out)
    
  }

}














