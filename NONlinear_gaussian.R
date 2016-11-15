##################################################################
# State space model - non linear
##################################################################

##################################################################
# Model specification
##################################################################

# SSM params
sigma2_V <- 10
sigma2_W <- 10
prior_par = c(0,sqrt(5)) # mean and sd of the prior distribution for x1

source("propagate_SSM.R")

# Problem specific functions 

nlinear_state_update = function(x_k, k, theta_state){
  theta_state = c(sigma2_V)
  error = rnorm(1, mean=0, sd=sqrt(sigma2_V))
  x_k1 = (x_k / 2) + 25 * (x_k / (1 + x_k^2)) + 8 * cos(1.2 * (k+1)) + error
  return (x_k1)  
}


nlinear_obs_update = function (x_k, theta_obs){
  theta_obs = c(sigma2_W)
  error = rnorm(1, mean=0, sd=sqrt(sigma2_W))
  y_k = (x_k^2 / 20) + error
  return (y_k)
}

calculate_weight_nlinear = function(observed_val, particles_vals, theta_obs){
  # This weight calculation is SSM dependant
  weight  = dnorm(observed_val, mean = particles_vals^2 / 20, sd=sqrt(theta_obs)) # weights here are from prob density g evaulated at y1|x_1 for all steps
>>>>>>> 667f19c93fdaa4306fd1eff64a5a4f4cc450c72a
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


plot_nonlinear=tracePlot(PIMH_nonlinear$state_values[1,], observed_process=observed_process_nlinear)


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

vector_particle_numbers = c(1, 10, 25, 50, 100, 500, 1000, 2000) #the values of N to loop over
vector_times = c(10,25,50,100)
acceptance_ratios = c()
acceptance_rate_df = data.frame(t=integer(0), N=integer(0), acceptance_rate=numeric(0))

data=generate_data(nlinear_state_update, nlinear_obs_update, prior_par, t=100, plot=FALSE, theta_state, theta_obs)
observed_process_nlinear = data$observed_process

t_i = 0

for (t in vector_times){
  i = 0
  for (N in vector_particle_numbers){
    cat("")
    i = i + 1
    t_i = t_i + 1

    n_iter = 1000 # 50,000 iterations were used in the paper
    
    PIMH_nonlinear = PIMH(n_iter, 
                          N, 
                          calculate_weight=calculate_weight_nlinear,  
                          state_update=nlinear_state_update, 
                          observed_process=observed_process_nlinear,
                          theta_state = c(sigma2_V),
                          theta_obs = c(sigma2_W))
    
    acceptance_rate_df[t_i, ] <- c(t, N, PIMH_nonlinear$acceptance_ratio)
  }
}

acceptance_rate_df$T = as.factor(acceptance_rate_df$t)

acceptance_plot = ggplot(acceptance_rate_df, aes(x = N, y = acceptance_rate, group = T, color=T)) + 
  geom_point(aes(shape=T), size=2) + geom_line()

data_path = paste('plots/nonlinear_guassian_acceptance_rate', '_n_iter_', n_iter, '_' , Sys.time(), '.rdata', sep='')
save(acceptance_rate_df, file = data_path)

plot_path = paste('plots/nonlinear_guassian_acceptance_rate', '_n_iter_', n_iter, '_' , Sys.time(), '.pdf', sep='')

ggsave(plot_path, acceptance_plot)
