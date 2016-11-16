##################################################################
# State space model - non linear
##################################################################

##################################################################
# Model specification
##################################################################
# SSM params
sigma2_V <- 10
sigma2_W <- 10
theta_obs = c(sigma2_W)
theta_state = c(sigma2_V)

prior_par = c(0,sqrt(5)) # mean and sd of the prior distribution for x1

source("propagate_SSM.R")

# Problem specific functions 

nlinear_state_update = function(x_k, k, theta_state){
  error = rnorm(1, mean=0, sd=sqrt(theta_state[1]))
  x_k1 = (x_k / 2) + 25 * (x_k / (1 + x_k^2)) + 8 * cos(1.2 * (k+1)) + error
  return (x_k1)  
}


nlinear_obs_update = function (x_k, theta_obs){
  error = rnorm(1, mean=0, sd=sqrt(theta_obs[1]))
  y_k = (x_k^2 / 20) + error
  return (y_k)
}

calculate_weight_nlinear = function(observed_val, particles_vals, theta_obs){
  # This weight calculation is SSM dependant
  weight  = dnorm(observed_val, mean = particles_vals^2 / 20, sd=sqrt(theta_obs[1])) # weights here are from prob density g evaulated at y1|x_1 for all steps
}

calculate_weight_nlinear_LOG = function(observed_val, particles_vals, theta_obs){
  # This weight calculation is SSM dependant
  weight  = (-(observed_val-particles_vals)^2)/(2*theta_obs[1])-log(sqrt(2*theta_obs[1]*pi)) # LOG weights here are from prob density g evaulated at y1|x_1 for all steps
  }


##################################################################
# Data generation
##################################################################

t = 50      # number of time steps for the data generation

data=generate_data(nlinear_state_update,nlinear_obs_update, prior_par, t, plot=TRUE, theta_state, theta_obs)
latent_process_nlinear = data$latent_process
observed_process_nlinear = data$observed_process
data$plot


##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################
source("SMC.R")

N = 100  #number of particles 

#run the SMC for the non linear gaussian state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_nlinear, state_update=nlinear_state_update, observed_process=observed_process_nlinear, theta_state= c(sigma2_V), theta_obs= c(sigma2_W))

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(smc_output$particles_in_time, latent_process_nlinear)

# save plot
#plot_path = paste('plots/nonlinear_SMC_particles_in_time', '_N_', n_iter, '_' , format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), '.pdf', sep='')
#ggsave(plot_path)

##################################################################
# PIMH for a non linear gaussian model
##################################################################
source("PIMH.R")

n_iter = 100


PIMH_nonlinear = PIMH(n_iter, 
                      N, 
                      calculate_weight=calculate_weight_nlinear,  
                      state_update=nlinear_state_update, 
                      observed_process=observed_process_nlinear,
                      theta_state = c(sigma2_V),
                      theta_obs = c(sigma2_W))


plot_nonlinear=trace_plot(PIMH_nonlinear$state_values[1,], n_iter, title = "Trajectory for the firt particle")
plot_nonlinear


##################################################################
# PMMH for a non linear gaussian model
##################################################################
# source("PIMH.R")
# 
# d = 1      #update only phi
# 
# set.seed(27)
# 
# PMMH_nonlinear = PMMH_nonlinear(n_iter, 
#                       N, 
#                       d,
#                       calculate_weight=calculate_weight_nlinear,  
#                       state_update=nlinear_state_update, 
#                       observed_process=observed_process_nlinear,
#                       theta_state = c(sigma2_V),
#                       theta_obs = c(sigma2_W))
# 



##################################################################
# SMC-MCMC for a linear gaussian model - acceptence rate w/ N
# This is slow to run!
##################################################################
source("PIMH.R")
source("propagate_SSM.R")

sigma2_V <- 10
sigma2_W <- 1 #1 or 10
theta_obs = c(sigma2_W)
theta_state = c(sigma2_V)

vector_particle_numbers = c(1, 10, 25, 50, 100, 500, 1000, 2000) #the values of N to loop over
vector_times = c(10,25,50,100) #the times we want to calculate acceptance averages for
acceptance_rate_df = data.frame(t=integer(0), N=integer(0), acceptance_rate=numeric(0)) #to store the data to plot

#generate one set of data; subset for different ts later. t here must be biggest we want to compare
data=generate_data(nlinear_state_update, nlinear_obs_update, prior_par, t=100, plot=TRUE, theta_state, theta_obs)
observed_process_nlinear = data$observed_process
data$plot

# save data for plot
data_path = paste('plots/nonlinear_guassian_ssm_sigma2_W_', sigma2_W,  format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), '.rdata', sep='')
save(data, file = data_path)

# save plot
plot_path = paste('plots/nonlinear_guassian_ssm_sigma2_W_', sigma2_W,  format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), '.pdf', sep='')
ggsave(plot_path, plot = last_plot())

t_i = 0
for (t in vector_times){
  cat(t, " is t\n")
  i = 0
  for (N in vector_particle_numbers){
    cat(N, " is N\n")
    i = i + 1
    t_i = t_i + 1

    n_iter = 10 # 50,000 iterations were used in the paper, 1000 takes on the order of one hour
    
    PIMH_nonlinear = PIMH(n_iter, 
                          N, 
                          calculate_weight=calculate_weight_nlinear_LOG,  
                          state_update=nlinear_state_update, 
                          observed_process=observed_process_nlinear[1:t],
                          theta_state = c(sigma2_V),
                          theta_obs = c(sigma2_W))
    
    acceptance_rate_df[t_i, ] <- c(t, N, PIMH_nonlinear$acceptance_ratio)
  }
}

# create acceptance rate plot
acceptance_rate_df$T = as.factor(acceptance_rate_df$t)

acceptance_plot = ggplot(acceptance_rate_df, aes(x = N, y = acceptance_rate, group = T, color=T)) + 
  geom_point(aes(shape=T), size=2) + geom_line() + theme_grey(base_size = 18)

# save data for plot
data_path = paste('plots/nonlinear_guassian_acceptance_rate__sigma2_W_', sigma2_W, '_n_iter_', n_iter, '_' , format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), '.rdata', sep='')
save(acceptance_rate_df, file = data_path)

# save plot
plot_path = paste('plots/nonlinear_guassian_acceptance_rate__sigma2_W_', sigma2_W , '_n_iter_', n_iter, '_' , format(Sys.time(), "%Y_%m_%d_%H_%M_%S"), '.pdf', sep='')
ggsave(plot_path, acceptance_plot)
