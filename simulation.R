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
#source('functions_260124.R')
source('functions_120624.R')
#---------------
# SIMULATE DATA
#---------------
# Number of individuals we 'scan' (simulate)
N_ind <- 1000
# Number of observations per individual:
N <- 100
# Vector of brain-symptom correlations (we assume 8 involved brain regions):
# High Measurement Error
vec_r_S <-  c(-0.3, -0.21, -0.12, -0.03, 0.03, 0.12, 0.21, 0.3)
# Naming our variables
variables <-  c("S", "B1","B2","B3","B4","B5","B6","B7","B8")

#-----------------------------
# RQ1
#-----------------------------
# HIGH MR
#--------
# here, we call the function 'make_data_mr(..., high = T)' N_ind times to simulate a total of 
# N*N_ind (100*1000 = 100.000) observations with high MR:
set.seed(1234)
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
# 100000 (N_ind * N) models in the background):

out_individ_high <- make_diagnostics(data_full, individual = T)

out_msmt_high <- make_diagnostics(data_full, time = T)
#--------
# LOW MR
#--------
set.seed(1234)
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

out_individ_low <- make_diagnostics(data_full, individual = T)

out_msmt_low <- make_diagnostics(data_full, time = T)


#---------
#MEDIUM MR
#---------
set.seed(1234)
data <- replicate(N_ind, make_data_mr(N, vec_r_S, 0, variables, medium = T), simplify = F)

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

out_individ_med <- make_diagnostics(data_full, individual = T)

out_msmt_med <- make_diagnostics(data_full, time = T)


#----------------------
# Low Measurement Error
#----------------------
vec_r_S_low <- c(-0.54, -0.38, -0.22, -0.05, 0.05, 0.22, 0.38, 0.54)
#--------
# HIGH MR
#--------
# here, we call the function 'make_data_mr(..., high = T)' N_ind times to simulate a total of 
# N*N_ind (100*1000 = 100.000) observations with high MR:
set.seed(1234)
data <- replicate(N_ind, make_data_mr(N, cor_bs = vec_r_S_low, cor_bb = 0, var_names = variables, high = T), simplify = F)

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
# 100000 (N_ind * N) models in the background):

out_individ_high2 <- make_diagnostics(data_full, individual = T)
out_msmt_high2 <- make_diagnostics(data_full, time = T)
#--------
# LOW MR
#--------
set.seed(1234)
data <- replicate(N_ind, make_data_mr(N, vec_r_S_low, 0, variables, low = T), simplify = F)

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

out_individ_low2 <- make_diagnostics(data_full, individual = T)
out_msmt_low2 <- make_diagnostics(data_full, time = T)


#---------
#MEDIUM MR
#---------
set.seed(1234)
data <- replicate(N_ind, make_data_mr(N, vec_r_S_low, 0, variables, medium = T), simplify = F)

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

out_individ_med2 <- make_diagnostics(data_full, individual = T)
out_msmt_med2 <- make_diagnostics(data_full, time = T)

#------
# PLOTS
#------
df <- data.frame(
  dataset = rep(c("High MR", "No MR", "Medium MR"), each = 2),
  variable = rep(c("Individual", "Group"), 3),
  msmt_error = rep(c('High ME', 'Low ME'), each = 6),
  mean = c(out_individ_high$mean_ca_individ, out_msmt_high$mean_ca_msmt, 
           out_individ_low$mean_ca_individ, out_msmt_low$mean_ca_msmt,
           out_individ_med$mean_ca_individ, out_msmt_med$mean_ca_msmt,
           out_individ_high2$mean_ca_individ, out_msmt_high2$mean_ca_msmt, 
           out_individ_low2$mean_ca_individ, out_msmt_low2$mean_ca_msmt,
           out_individ_med2$mean_ca_individ, out_msmt_med2$mean_ca_msmt),
  lowerCI = c(c(t.test(out_individ_high$ca_individ)$conf.int)[1], c(t.test(out_msmt_high$ca_msmt)$conf.int)[1],
              c(t.test(out_individ_low$ca_individ)$conf.int)[1], c(t.test(out_msmt_low$ca_msmt)$conf.int)[1],
              c(t.test(out_individ_med$ca_individ)$conf.int)[1], c(t.test(out_msmt_med$ca_msmt)$conf.int)[1],
              c(t.test(out_individ_high2$ca_individ)$conf.int)[1], c(t.test(out_msmt_high2$ca_msmt)$conf.int)[1],
              c(t.test(out_individ_low2$ca_individ)$conf.int)[1], c(t.test(out_msmt_low2$ca_msmt)$conf.int)[1],
              c(t.test(out_individ_med2$ca_individ)$conf.int)[1], c(t.test(out_msmt_med2$ca_msmt)$conf.int)[1]),
  upperCI = c(c(t.test(out_individ_high$ca_individ)$conf.int)[2], c(t.test(out_msmt_high$ca_msmt)$conf.int)[2],
              c(t.test(out_individ_low$ca_individ)$conf.int)[2], c(t.test(out_msmt_low$ca_msmt)$conf.int)[2],
              c(t.test(out_individ_med$ca_individ)$conf.int)[2], c(t.test(out_msmt_med$ca_msmt)$conf.int)[2],
              c(t.test(out_individ_high2$ca_individ)$conf.int)[2], c(t.test(out_msmt_high2$ca_msmt)$conf.int)[2],
              c(t.test(out_individ_low2$ca_individ)$conf.int)[2], c(t.test(out_msmt_low2$ca_msmt)$conf.int)[2],
              c(t.test(out_individ_med2$ca_individ)$conf.int)[2], c(t.test(out_msmt_med2$ca_msmt)$conf.int)[2])
)

df$dataset <- factor(df$dataset, levels = c("No MR", "Medium MR", "High MR"))
df$variable <- factor(df$variable, levels = c("Individual", "Group"))
df$msmt_error <- factor(df$msmt_error, levels = c('High ME', 'Low ME'))




# Create the plot
tiff("ca_plot.tiff", units="in", width=10, height=8, res=300)
(p <- ggplot(data = df, aes(x = variable, y = mean, fill = dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2, position = position_dodge(width = 0.5)) +
    labs(x = "Data", y = "Mean Classification Accuracy", fill = "") +
    #ggtitle("Classification Accuracy Means") +
    theme_minimal() +
    theme(legend.title = element_blank())+
    coord_cartesian(ylim = c(0.45, 1)) +
    scale_fill_viridis(discrete = TRUE) +
    geom_text(aes(label = round(mean, 2), y = mean + 0.02, group = dataset), 
              position = position_dodge(width = 0.5), vjust = -0.5, 
              size = 3, angle = 30) +
    facet_wrap(~msmt_error))
dev.off()


#-------------------------------------
# RQ 2: NETWORKS AS EMERGENT REALIZERS
#-------------------------------------
#First we set values that are constant across variations:

N <- 100
N_ind <- 1000
#Brain region variables, column names:
variables <- paste("B", 1:20, sep = "")
#Average correlation 
net_cor1 <-  runif(N_ind, 0.29, 0.31)

#----------
# LOW MR
#----------

net_low <- make_diagnostics_net(N_ind, N, n_var = 20, net_cor = net_cor1, var_names = variables,
                                error_mean = 0, error_sd = 0.3, low = T)

#----------
# MEDIUM MR
#----------
net_medium <- make_diagnostics_net(N_ind, N, n_var = 20, net_cor = net_cor1, var_names = variables,
                                   error_mean = 0, error_sd = 0.3, medium = T)

#----------
# HIGH MR
#----------
net_high <- make_diagnostics_net(N_ind, N, n_var = 20, net_cor = net_cor1, var_names = variables,
                                 error_mean = 0, error_sd = 0.3, high = T)


# PLOTS
df_net <- data.frame(
  dataset = rep(c("High MR", "No MR", "Medium MR"), each = 2),
  variable = rep(c("Individual", "Group"), 3),
  level = rep(c('Region', 'Network'), each = 6),
  mean = c(net_high$mean_ca_reg_id, net_high$mean_ca_reg_t, 
           net_low$mean_ca_reg_id, net_low$mean_ca_reg_t,
           net_medium$mean_ca_reg_id, net_medium$mean_ca_reg_t,
           net_high$mean_ca_net_id, net_high$mean_ca_net_t,
           net_low$mean_ca_net_id, net_low$mean_ca_net_t,
           net_medium$mean_ca_net_id, net_medium$mean_ca_net_t),
  lowerCI = c(c(t.test(net_high$ca_reg_id)$conf.int)[1], c(t.test(unlist(net_high$ca_reg_t))$conf.int)[1],
              c(t.test(net_low$ca_reg_id)$conf.int)[1], c(t.test(unlist(net_low$ca_reg_t))$conf.int)[1],
              c(t.test(net_medium$ca_reg_id)$conf.int)[1], c(t.test(unlist(net_medium$ca_reg_t))$conf.int)[1],
              c(t.test(net_high$ca_net_id)$conf.int)[1], c(t.test(unlist(net_high$ca_net_t))$conf.int)[1],
              c(t.test(net_low$ca_net_id)$conf.int)[1], c(t.test(unlist(net_low$ca_net_t))$conf.int)[1],
              c(t.test(net_medium$ca_net_id)$conf.int)[1], c(t.test(unlist(net_medium$ca_net_t))$conf.int)[1]),
  upperCI = c(c(t.test(net_high$ca_reg_id)$conf.int)[2], c(t.test(unlist(net_high$ca_reg_t))$conf.int)[2],
              c(t.test(net_low$ca_reg_id)$conf.int)[2], c(t.test(unlist(net_low$ca_reg_t))$conf.int)[2],
              c(t.test(net_medium$ca_reg_id)$conf.int)[2], c(t.test(unlist(net_medium$ca_reg_t))$conf.int)[2],
              c(t.test(net_high$ca_net_id)$conf.int)[2], c(t.test(unlist(net_high$ca_net_t))$conf.int)[2],
              c(t.test(net_low$ca_net_id)$conf.int)[2], c(t.test(unlist(net_low$ca_net_t))$conf.int)[2],
              c(t.test(net_medium$ca_net_id)$conf.int)[2], c(t.test(unlist(net_medium$ca_net_t))$conf.int)[2])
)

df_net$dataset <- factor(df_net$dataset, levels = c("No MR", "Medium MR", "High MR"))
df_net$variable <- factor(df_net$variable, levels = c("Individual", "Group"))
df_net$level <- factor(df_net$level, levels = c("Region", "Network"))



# Create the plot
tiff("ca_net_plot.tiff", units="in", width=8, height=5, res=300)
(p <- ggplot(data = df_net, aes(x = variable, y = mean, fill = dataset)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.5) +
    geom_errorbar(aes(ymin = lowerCI, ymax = upperCI), width = 0.2, position = position_dodge(width = 0.5)) +
    labs(x = "Data", y = "Mean Classification Accuracy", fill = "") +
    #ggtitle("Classification Accuracy Means") +
    theme_minimal() +
    theme(legend.title = element_blank())+
    coord_cartesian(ylim = c(0.45, 1)) +
    scale_fill_viridis(discrete = TRUE) +
    geom_text(aes(label = round(mean, 2), y = mean + 0.02, group = dataset), 
              position = position_dodge(width = 0.5), vjust = -0.5, 
              size = 3, angle = 30) +
    facet_wrap(~level))
dev.off()

