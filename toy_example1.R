library(ggplot2)
library(reshape2)

##################################################################
# State space model - toy example
##################################################################
##################################################################
# Generate some observations using the SSM
##################################################################


t = 100 # num of iters

# SSM params
phi = 0.95
mu = 0.9
sigma_2 = 100
sigma = sqrt(sigma_2)
eta_2 = (sigma_2/(1-phi^2))/5
eta = sqrt(eta_2)


toy_state_update = function(x_k){
  error_value = rnorm(1, mean=0, sd=1)
  x_k1 = mu + phi*(x_k-mu) + sigma*error_value
  return (x_k1)  
}

toy_obs_update = function (x_k){
  error_value = rnorm(1, mean=0, sd=1)
  y_k = x_k + eta*error_value
  return (y_k)
}

# initialise process X (latent) and Y (observed)
latent_process_toy = vector(length = t) # to store process X in time
observed_process_toy= vector(length = t) # to store observations associated with X
x1 = rnorm(1, mean=0, sd=1) # pick starting point from rnorm 
latent_process_toy[1] = x1 # store starting point
observed_process_toy[1] = toy_obs_update(x_k = latent_process_toy[1]) # store first observation

# propogate X (latent) and Y (observed) in time by SSM
for (i in 2:t){
  latent_process_toy[i] = toy_state_update(x_k = latent_process_toy[i-1])
  observed_process_toy[i] = toy_obs_update(x_k = latent_process_toy[i])
}

# plot the evolution of the process
process_dataframe = data.frame(observed=observed_process_toy, latent=latent_process_toy)
plot_processes_in_time(process_dataframe)

# specify function to calculate weight for this SSM
calculate_weight_toy = function(observed_val, particles_vals){
  # This weight calculation is SSM dependant
  prob_argument = (observed_val-particles_vals)/eta
  weight  = dnorm(prob_argument, mean=0, sd=1) # weights here are from prob density g evaulated at y1|x_1 for all steps
}

##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################

N = 10  #number of particles 

#run the SMC for the state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_toy, state_update=toy_state_update, observed_process=observed_process_toy)
particles_in_time = smc_output$particles_in_time

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(particles_in_time, latent_process)
