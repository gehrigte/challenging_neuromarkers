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
# install.packages('foo'), and replace foo with the respective package name.

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
library(matrixStats)
library(gridExtra)
library(splitstackshape)

#---------------------
# SIMULATION FUNCTIONS
#---------------------
# RESEARCH QUESTION 1
# To define the following functions, we need to set 
# N = Number of observations per individual
# cor_bs = correlation between brain region and symptom
# cor_bb = mean correlation between brain regions
# when calling the function.

# Our simulation is based on 8 brain regions, meaning we pre-define a vector
# of brain-symptom correlation with length 8. This can be easily adapted.

make_data_mr <- function(N, cor_bs, cor_bb, var_names, low = NULL, medium = NULL, high = NULL){ 
  if(is.null(low)&is.null(medium)&is.null(high)){
    stop("Please choose to simulate data with either low, medium, or high multiple realizability (MR).")
  }
  p <- length(cor_bs)
  if (!is.null(high)) {
    # HIGH MULTIPLE REALIZABILITY (MR), i.e., activity in random brain regions can 
    # elicit the same symptom.
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
   } else if (!is.null(low)) {
     # LOW MR, i.e., activity in same brain regions of different individuals
     # elicit same symptom.
      # no randomness, i.e., the order of brain-symptom corr is fixed and we don't sample
      s <- cor_bs
      #we first create a PD matrix
      m <- matrix(cor_bb, ncol = p+1, nrow = p+1)
      diag(m) <- 1
      s1 <- c(1, s)
      m[ ,1] <- m[1, ] <- s1
      colnames(m) <- rownames(m) <-  var_names
   } else if (!is.null(medium)) {
     # MEDIUM MR, i.e., some brain regions remain the same across individuals while
     # some vary to elicit same symptom.
     s_neg <- sample(cor_bs[cor_bs < 0], replace = F)
     s_pos <- sample(cor_bs[cor_bs > 0], replace = F)
     s <- c(s_neg, s_pos)
     
     m <- matrix(cor_bb, ncol = p+1, nrow = p+1)
     diag(m) <- 1
     s1 <- c(1, s)
     m[1, ] <- m[, 1] <- s1
     colnames(m) <- rownames(m) <-  var_names
   }
  data <- rnorm_multi(n = N, vars = p+1, r = m, varnames = var_names)
  # usually the numbers are not exactly as intended, but we can decide which 
  # deviations are acceptable
  # here we make sure the first row differs from the intended variance and mean (var 0.01, mean 0.2)
  # max by .004/.005, respectively:
  while(((abs(mean(cor(data)[2:9,1]) - mean(cor_bs))) > 0.005) | 
        (abs(var(cor(data)[2:9,1]) - var(cor_bs)) > 0.004)) {
    data <- rnorm_multi(n = N, vars = p + 1, r = m)
  }
  return(data)
}


# RESEARCH QUESTION 2

# The following function 'make_data_network' allows to simulate network data, i.e., 
# activation in different brain regions are summarized to make a network. The 
# activity of (parts of) the network yields the symptom variable S.

make_data_network <- function(N, net_cor, n_var, var_names, error_mean, error_sd, 
                              low = NULL, medium = NULL, high = NULL){
  if(is.null(low)&is.null(medium)&is.null(high)){
    stop("Please choose to simulate data with either low, medium, or high multiple realizability (MR).")
  }
  
  error <- rnorm(N, error_mean, error_sd)
  #dataset for 1 individual, N measurements, n_var brain areas
  p <-  n_var
  m <- matrix(sample(net_cor, p*p, replace = TRUE), ncol = p, nrow = p)
  diag(m) <- 1
  m[((p/4)+1):p, 1:(p/4)] <-  m[1:floor(p/4), (floor(p/4)+1):p] <-  m[(floor(p/4)+1):(p/2), (floor(p/2)+1):p] <-  m[(floor(p/2)+1):p, (floor(p/4)+1):(p/2)] <-  m[(floor(3/4*p)+1):p, (floor(p/2)+1): (3/4*p)] <-  m[(floor(p/2)+1): (3/4*p),(floor(3/4*p)+1):p] <- 0
  mat <- makePsd(m, method = "correlation")
  mat[lower.tri(mat)] = t(mat)[lower.tri(mat)]
  colnames(mat) <- var_names
  
  data <-  rnorm_multi(n = N, vars = p, r = mat, varnames = var_names, empirical = FALSE)
  
  #create network variables
  data_network <- data %>%
    mutate(net1 = abs(rowProds(as.matrix(data[, 1:floor(p/4)]))),
           net2 = abs(rowProds(as.matrix(data[, (floor(p/4)+1):floor(p/2)]))),
           net3 = abs(rowProds(as.matrix(data[, (floor(p/2)+1):floor(3/4*p)]))),
           net4 = abs(rowProds(as.matrix(data[, (floor(3/4*p)+1):p])))) 
  
  if (!is.null(low)){
    data_network <- data_network %>% 
      mutate(S = rowSums(data_network[(p+1):(p+4)]) + sample(error, 1)) %>% 
      mutate_all(~(scale(.)
                   %>% as.vector))
    return(data_network)
  }  else if (!is.null(medium)) {
    net_vars <- c("net2", "net3", "net4")
    # randomly determine how many variables to include in S
    n_vars <- sample(1:3, 1)
    # randomly determine which variables to include in S
    included_vars <- sample(net_vars, n_vars, replace = FALSE)
    data_network <- data_network %>%
      mutate(S = rowSums(data_network[included_vars]) + rowSums(data_network[(floor(p)+1)]) + 
               sample(error, 1))%>%  # random sample of networks 2-4, network 1 always included + error
      mutate_all(~(scale(.) %>% 
                     as.vector))
    return(data_network)
  }  else if (!is.null(high)) {
     net_vars <- c("net1", "net2", "net3", "net4")
  
  n_vars <- sample(1:4, 1)
  included_vars <- sample(net_vars, n_vars, replace = FALSE)
  # compute S as the sum of randomly selected and included variables + error
  data_network <- data_network %>%
    mutate(S = rowSums(data_network[included_vars]) + sample(error, 1)) %>% 
    mutate_all(~(scale(.) %>% 
                   as.vector))
  return(data_network)
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



make_diagnostics <- function(data, individual = NULL, time = NULL) {
  if (is.null(individual) & is.null(time) ) {
    stop("Choose either individual-level analysis across time ('individual = T')
         or group-level analysis at one time point ('time = T'). 
         You must select one level of analysis.")
  }
    if (!is.null(individual) && individual == TRUE) {
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
      
      # Create training (80% of data) and test (20% of data) subsets:
      data_train <- data_full %>% 
        group_by(ID) %>% 
        slice_sample(prop = .8)
      
      data_test <- data_full %>% 
        group_by(ID) %>% 
        anti_join(data_train, by = 'rownum')  # select remaining observations that are NOT in training data
      
      # Train model:
      trained_models <- data_train %>% 
        do(prediction_mods = train(form = D ~ B1 + B2 + B3 + B4 + B5 + B6 + B7 + B8,
                                   data = .,
                                   trControl = trainControl(method = "cv", number = 5),
                                   method = "glm",
                                   family = "binomial"))
      
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
      out <- list('rsq_indiv' = c(unlist(r_sq_id)),
                  'mean_rsq_indiv' = mean(r_sq_id),
                  'ca_individ' = c(unlist(ca_id)),
                  'mean_ca_individ' = mean(unlist(ca_id)))
      
      return(out)
    } else if (!is.null(time) && time == T) {
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
                                   family = "binomial"))

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
      out <- list('rsq_msmt' = c(unlist(r_sq_msmt)),
                  'mean_rsq_msmt' = mean(r_sq_msmt),
                  'ca_msmt' = c(unlist(ca_msmt)),
                  'mean_ca_msmt' = mean(unlist(ca_msmt)))
      return(out)
    }
 }
  





make_diagnostics_net <- function(N_ind, ...) {
  # be sure to select low, medium, or high MR when calling make_diagnostics_net function
  input_network <- list(...)
  # Prediction using regions
  trained_models_id <- rep(NA, N_ind)
  predicted_reg_id <- rep(NA, N_ind)
  ca_reg_id <- rep(NA, N_ind)
  
  trained_models_net_id <- rep(NA, N_ind)
 
  predicted_net_id <- rep(NA, N_ind)
  ca_net_id <- rep(NA, N_ind)

  for (i in 1:N_ind) {
      # repeat process for all N_ind individuals
      data <- do.call(make_data_network, input_network) %>% 
        mutate(rownum = row_number(),
                t = rownum) %>%
        mutate(ID = i) %>% 
        mutate(D = factor(ifelse(S > median(S), 1, 0), levels = c(0, 1)))  %>% 
        select(rownum, ID, t, D, S, B1:B20, net1:net4)
        
      data_train_id <- stratified(data, group = c('D', 'ID'), .8)
      if (length(unique(unlist(data_train_id$D))) == 1) {
        # re-sample if we only have one diagnosis in our subsample
        data <- make_data_network(N, n_var = 20, net_cor = net_cor1, var_names = variables, 
                                  error_mean = 0, error_sd = 0.3, ...) %>% 
          mutate(rownum = row_number()) %>%
          mutate(ID = i) %>% 
          mutate(D = factor(ifelse(S > median(S), 1, 0), levels = c(0, 1))) %>% 
          select(rownum, ID, D, S, B1:B20, net1:net4)
        data_train_id <- stratified(data, group = c('D', 'ID'), .8)
        }
      index_train <- data_train_id$rownum
      data_test_id <- data[-index_train,] %>% 
        mutate(D = ifelse(D == 0, 'no_diagn', 'diagn'))
      
      # settings for prediction model:
      ctrlspecs <- trainControl(method="cv", 
                                  number=5, 
                                  savePredictions="all",
                                  classProbs=TRUE)
      
      trained_models_id[[i]] <- data_train_id %>% 
        select(D, B1:B20) %>% 
        mutate(D = ifelse(D == 0, 'no_diagn', 'diagn')) %>% 
        do(mods_id = train(D ~ ., data=., 
                             method="glm", 
                             family=binomial, 
                             trControl=ctrlspecs))
      trained_models_net_id[[i]] <- data_train_id %>% 
        select(D, net1:net4) %>% 
        mutate(D = ifelse(D == 0, 'no_diagn', 'diagn')) %>% 
        do(mods_id = train(D ~ ., data=., 
                           method="glm", 
                           family=binomial, 
                           trControl=ctrlspecs))
      predicted_reg_id[[i]] <- predict(trained_models_id[[i]][[1]], 
                                       newdata = data_test_id[data_test_id$ID ==i,])
      # CA = mean(predicted diagnosis == true diagnosis):
      ca_reg_id[[i]] <- mean(unlist(predicted_reg_id[[i]]) == 
                               unlist(data_test_id[data_test_id$ID == i, 'D'])) 
      predicted_net_id[[i]] <- predict(trained_models_net_id[[i]][[1]],
                                       newdata = data_test_id[data_test_id$ID == i,])
      ca_net_id[[i]] <- mean(unlist(predicted_net_id[[i]]) ==
                               unlist(data_test_id[data_test_id$ID == i, 'D']))
      
      print(paste0('Computing at ', c(i/N_ind*100), '%.'))
      } 
  data_t <- replicate(N_ind, do.call(make_data_network, input_network), simplify = F) %>%
    do.call(rbind.data.frame,.) %>%
    mutate(rownum = row_number()) %>%
    mutate(ID = rep(1:N_ind, each = N)) %>%
    mutate(D = factor(ifelse(S > median(S), 1, 0), levels = c(0, 1))) %>%
    mutate(t = factor(rep(1:N, times = N_ind))) %>%
    select(rownum, ID, t, D, S, B1:B20, net1:net4)
  
  data_train_t <- stratified(data_t, group = c('D', 't'), 0.8)
  data_test_t <- data_t%>%
    anti_join(data_train_t, by = 'rownum') %>%
    mutate(D = ifelse(D == 0, 'no_diagn', 'diagn'))
  
  # pre-define output vectors
  predicted_reg_t <- list(rep(NA, N))
  predicted_net_t <- list(rep(NA, N))
  ca_reg_t <- list(rep(NA, N))
  ca_net_t <- list(rep(NA, N))
  
  # Grouped - Region
  trained_models_t <- data_train_t %>%
    group_by(t) %>%
    select(D, B1:B20) %>%
    mutate(D = ifelse(D == 0, 'no_diagn', 'diagn')) %>%
    do(mods_reg_t = train(form = D ~ .,
                          data = .,
                          trControl = trainControl(method = "cv", number = 5),
                          method = "glm",
                          family = "binomial"))
      
  # Grouped - Network
  trained_models_net_t <- data_train_t %>%
    group_by(t) %>%
    select(D, net1:net4) %>%
    mutate(D = ifelse(D == 0, 'no_diagn', 'diagn')) %>%
    do(mods_net_t = train(D ~ ., data=.,
                          method="glm",
                          family=binomial,
                          trControl=ctrlspecs))
  for (m in 1: N) {
    predicted_reg_t[[m]] <- predict(trained_models_t$mods_reg_t[[m]],
                                    newdata = data_test_t[data_test_t$t == m,])
    # CA = mean(predicted diagnosis == true diagnosis):
    ca_reg_t[[m]] <- mean(unlist(predicted_reg_t[[m]]) ==
                            unlist(data_test_t[data_test_t$t == m, 'D']))
    
    predicted_net_t[[m]] <- predict(trained_models_net_t$mods_net_t[[m]],
                                    newdata = data_test_t[data_test_t$t == m,])
    
    ca_net_t[[m]] <- mean(unlist(predicted_net_t[[m]]) ==
                            unlist(data_test_t[data_test_t$t == m, 'D']))
    }
  return(list('ca_reg_id' = ca_reg_id,
              'mean_ca_reg_id' = mean(ca_reg_id),
              'ca_net_id' = ca_net_id,
              'mean_ca_net_id' = mean(ca_net_id),
              'ca_reg_t' = unlist(ca_reg_t),
              'mean_ca_reg_t' = mean(unlist(ca_reg_t)),
              'ca_net_t' = unlist(ca_net_t),
              'mean_ca_net_t' = mean(unlist(ca_net_t))))
  }
      



