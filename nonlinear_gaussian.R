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

nlinear_state_update = function(x_k, k){
  error = rnorm(1, mean=0, sd=sqrt(sigma2_V))
  x_k1 = (x_k / 2) + 25 * (x_k / (1 + x_k^2)) + 8 * cos(1.2 * (k+1)) + error
  return (x_k1)  
}

nlinear_obs_update = function (x_k){
  error = rnorm(1, mean=0, sd=sqrt(sigma2_W))
  y_k = (x_k^2 / 20) + error
  return (y_k)
}

calculate_weight_nlinear = function(observed_val, particles_vals){
  # This weight calculation is SSM dependant
  weight  = dnorm(observed_val, mean = particles_vals^2 / 20, sd=sqrt(sigma2_W)) # weights here are from prob density g evaulated at y1|x_1 for all steps
}


##################################################################
# Data generation
##################################################################

t = 50      # number of time steps for the data generation

data=generate_data(nlinear_state_update,nlinear_obs_update, prior_par, t, plot=TRUE)
latent_process_nlinear = data$latent_process
observed_process_nlinear = data$observed_process
data$plot


##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################
source("SMC.R")

N = 100  #number of particles 

#run the SMC for the non linear gaussian state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_nlinear, state_update=nlinear_state_update, observed_process=observed_process_nlinear)

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(smc_output$particles_in_time, latent_process_nlinear)


##################################################################
# SMC-MCMC for a linear gaussian model
##################################################################
source("PIMH.R")

n_iter = 1000

PIMH_nonlinear = PIMH(n_iter, 
                   N, 
                   calculate_weight=calculate_weight_nlinear,  
                   state_update=nlinear_state_update, 
                   observed_process=observed_process_nlinear)


plot_nonlinear=tracePlot(PIMH_nonlinear$state_values[1,], 
                n_iter, 
                title = "Trajectory for the firt particle")
plot_nonlinear