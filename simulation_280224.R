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
N_ind <- 10
# Number of observations per individual:
N <- 100
# Vector of brain-symptom correlations (we assume 8 involved brain regions):
vec_r_S <-  c(-0.3, -0.21, -0.12, -0.03, 0.03, 0.12, 0.21, 0.3)
# Naming our variables
variables <-  c("S", "B1","B2","B3","B4","B5","B6","B7","B8")

#-----------------------------
# RQ1
#-----------------------------
# HIGH MR
# here, we call the function 'make_data_mr(..., high = T)' N_ind times to simulate a total of 
# N*N_ind (100*1000 = 100.000) observations with high MR:
data <- replicate(N_ind, make_data_mr(N, cor_bs = vec_r_S, cor_bb = 0, var_names = variables, high = T), simplify = F)

# and add information on the individual, measurement number (time point), and the
# diagnosis, that is given, if the symptom variable exceeds the median symptom 
# across all time points and individuals:

data_full <- do.call(rbind.data.frame, data) %>% 
  mutate('ID' = factor(rep(c(1:N_ind), each = N), levels = c(1:N_ind)),
         't' = factor(rep(c(1:100), times = N_ind), levels = c(1:N)),
         'rownum' = factor(row_number()),
         'D' = ifelse(S > median(S), 1, 0), 
         'D' = as.factor(D)) #%>% 
 # select(rownum, ID, t, D, S, B1:B8)

# Finally, we calculate model diagnostics (this may take a while, as R calculates
# 100000 models in the background):

out_individ <- make_diagnostics(data_full, individual = T)

out_msmt <- make_diagnostics(data_full, individual = F)

# LOW MR
data <- replicate(N_ind, make_data_mr(N, vec_r_S, 0, variables, low = T), simplify = F)

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

#-------------------------------------
# RQ 2: NETWORKS AS EMERGENT REALIZERS
#-------------------------------------
#First we set values that are constant across variations:
#Brain region variables, column names:
variables <- paste("B", 1:20, sep = "")
#Average correlation 
net_cor1 <-  runif(1000, 0.29, 0.31)

N_ind <- 10
data <- replicate(N_ind, make_data_network(N, net_cor = net_cor1, var_names = variables, 
                                        error_mean = 0, error_sd = 0.3, low = T), simplify = F)

data_full <- do.call(rbind.data.frame, data) %>% 
  mutate('ID' = factor(rep(c(1:N_ind), each = N), levels = c(1:N_ind)),
         't' = factor(rep(c(1:100), times = N_ind), levels = c(1:N)),
         'rownum' = factor(row_number()),
         'D' = ifelse(S > median(S), 1, 0), 
         'D' = as.factor(D))
(out_net <- make_diagnostics(data_full, network = T))


