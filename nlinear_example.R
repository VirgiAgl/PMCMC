##################################################################
# State space model - non linear
##################################################################
##################################################################
# Generate some observations using the SSM
##################################################################


t = 100 # num of iters

# SSM params
sigma2_V <- 1
sigma2_W <- 1


nlinear_state_update = function(x_k, k){
  error_value = rnorm(1, mean=0, sd=1)
  x_k1 = (x_k / 2) + 25 * (x_k / (1 + x_k^2)) + 8 * cos(1.2 * (k+1)) + sqrt(sigma2_V)*error_value
  return (x_k1)  
}

nlinear_obs_update = function (x_k){
  error_value = rnorm(1, mean=0, sd=1)
  y_k = (x_k^2 / 20) + sqrt(sigma2_W)*error_value
  return (y_k)
}

# initialise process X (latent) and Y (observed)
latent_process_nlinear = vector(length = t) # to store process X in time
observed_process_nlinear = vector(length = t) # to store observations associated with X
x1 = rnorm(1, mean=0, sd=sqrt(5)) # pick starting point from rnorm 
latent_process_nlinear[1] = x1 # store starting point
observed_process_nlinear[1] = nlinear_obs_update(x_k = latent_process_nlinear[1]) # store first observation

# propogate X (latent) and Y (observed) in time by SSM
for (i in 2:t){
  latent_process_nlinear[i] = nlinear_state_update(x_k = latent_process_nlinear[i-1], k = i - 1)
  observed_process_nlinear[i] = nlinear_obs_update(x_k = latent_process_nlinear[i])
}

# plot the evolution of the process
nlinear_process_dataframe = data.frame(observed=observed_process_nlinear, latent=latent_process_nlinear)
plot_processes_in_time(nlinear_process_dataframe)

# specify function to calculate weight for this SSM
calculate_weight_nlinear = function(observed_val, particles_vals){
  # This weight calculation is SSM dependant
  prob_argument = (observed_val - (particles_vals^2 / 20)) / sqrt(sigma2_W)
  weight  = dnorm(prob_argument, mean=0, sd=1) # weights here are from prob density g evaulated at y1|x_1 for all steps
}

##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################

N = 100  #number of particles 

#run the SMC for the state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_nlinear, state_update=nlinear_state_update, observed_process=observed_process_nlinear)
nlinear_particles_in_time = smc_output$particles_in_time

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(nlinear_particles_in_time, latent_process_nlinear)
