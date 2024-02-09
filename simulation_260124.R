# Challenging the Search for Neuromarkers of Mental Disorders
# Practical Script/ Part 2

# In this script, we define our variables of interest and simulate the data. We
# call our functions that we defined in the technical script.
# This script allows the user to play around with model parameters of interest
# (e.g., different vectors of brain-symptom correlations) and see how changes
# impact model parameters such as R squared or classification accuracy, both, on 
# an individual, or group level. 
# For the individual level, set make_diagnostics(..., individual = T),
# for the group level, set make_diagnostics(..., individual = F)

#------------------
# SOURCE FUNCTIONS
#------------------
# read in functions from theoretical script/part 1
source('functions_260124.R')

#---------------
# SIMULATE DATA
#---------------
# Number of individuals we 'scan' (simulate)
N_ind <- 1000
# Number of observations per individual:
N <- 100
# Vector of brain-symptom correlations (we assume 8 involved brain regions):
vec_r_S <-  c(-0.3, -0.21, -0.12, -0.03, 0.03, 0.12, 0.21, 0.3)
# Naming our variables
variables <-  c("S", "B1","B2","B3","B4","B5","B6","B7","B8")

#-----------------------------
# HIGH MULTIPLE REALIZABILITY
#-----------------------------
# here, we call the function 'make_high_mr' N_ind times to simulate a total of 
# N*N_ind (100*1000 = 100.000) observations:
data <- replicate(N_ind, make_high_mr(N, vec_r_S, 0, variables), simplify = F)

# and add information on the individual, measurement number (time point), and the
# diagnosis, that is given, if the symptom variable exceeds the median symptom 
# across all time points and individuals:

data_full <- do.call(rbind.data.frame, data) %>% 
  mutate('ID' = factor(rep(c(1:N_ind), each = N), levels = c(1:N_ind)),
         't' = factor(rep(c(1:100), times = N_ind), levels = c(1:N)),
         'rownum' = factor(row_number()),
         'D' = ifelse(S > median(S), 1, 0), 
         'D' = as.factor(D)) %>% 
  select(rownum, ID, t, D, S, B1:B8)

# Finally, we calculate model diagnostics (this may take a while, as R calculates
# 100000 models in the background):

out_individ <- make_diagnostics(data_full, individual = T)

out_msmt <- make_diagnostics(data_full, individual = F)


#-----------------------------
# LOW MULTIPLE REALIZABILITY
#-----------------------------
data <- replicate(N_ind, make_low_mr(N, vec_r_S, 0, variables), simplify = F)

# and add information on the individual, measurement number (time point), and the
# diagnosis, that is given, if the symptom variable exceeds the median symptom 
# across all time points and individuals:

data_full <- do.call(rbind.data.frame, data) %>% 
  mutate('ID' = factor(rep(c(1:N_ind), each = N), levels = c(1:N_ind)),
         't' = factor(rep(c(1:100), times = N_ind), levels = c(1:N)),
         'rownum' = factor(row_number()),
         'D' = ifelse(S > median(S), 1, 0), 
         'D' = as.factor(D)) %>% 
  select(rownum, ID, t, D, S, B1:B8)

# Finally, we calculate model diagnostics (this may take a while, as R calculates
# 100000 models in the background):

out_individ <- make_diagnostics(data_full, individual = T)
# out_individ$mean_rsq_indiv
# [1] 0.3141727
# out_individ$mean_ca_individ
# [1] 0.6352

out_msmt <- make_diagnostics(data_full, individual = F)
# out_msmt$mean_ca_msmt
# [1] 0.67535
# out_msmt$mean_rsq_msmt
# [1] 0.2938855


