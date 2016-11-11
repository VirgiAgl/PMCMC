library(ggplot2)
library(reshape2)

##################################################################
# Some plotting functions
##################################################################
plot_processes_in_time = function(process_dataframe){
  # Takes a dataframe of 1d processes, e.g. one
  # column each for observed and latent, and plots 
  # as lines in time with labels determined by column 
  # headings in process dataframe
  process_dataframe$time = 1:nrow(process_dataframe)
  melted_data_df = melt(process_dataframe, id.vars='time')
  ggplot(melted_data_df, aes(x=time, y = value, group = variable, colour = variable)) +   geom_line()
}

plot_particles_and_latent_process = function(particles_in_time, latent_process){
  data_df = data.frame(p=particles_in_time)
  data_df$time = 1:nrow(data_df)
  melted_data_df = melt(data_df, id.vars='time')
  melted_data_df$type = 'particles'
  latent_df = data.frame(l=latent_process)
  latent_df$time = 1:nrow(latent_df)
  melted_latent_df = melt(latent_df, id.vars='time')
  melted_latent_df$type = 'latent process'
  df = rbind(melted_data_df, melted_latent_df)
  head(df)
  ggplot(df, aes(x=time, y = value, group = variable, colour = variable, linetype=type, alpha=0.8)) +   geom_line() +guides(colour=FALSE, alpha=FALSE)
  
}


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

SMC = function(N, calculate_weight, state_update, observed_process){
  # Function to perform SMC for a state space model. 
  # N is desired number of particles
  # calculate_weight is a function to calculate the weight at each timestep
  # state_update is a function specifying the SSM
  t = length(observed_process)
  particles_in_time = matrix(NA, ncol=N, nrow= t) # N particles at t timesteps
  weights_in_time = matrix(NA, ncol=N, nrow= t) # N weights at t timesteps
  particles_in_time[1,] = rnorm(N, mean=0, sd=1) # Initialise with proposal density
  
  resample_count = 0
  
  for (i in 1:t){
    if (i >= 2) { # All steps other than 1st
      particles_in_time[i,] = sapply(X = particles_in_time[i-1,], FUN = state_update)#phi=phi, mu=mu, sigma=sigma)
    }
    
    weight = calculate_weight(observed_val = observed_process[i], particles_vals = particles_in_time[i,])
    weight_norm = weight/sum(weight)
    weights_in_time[i,] = weight_norm
    
    ESS = sum((weights_in_time[i,])^2)^-1 
    if (ESS<N/2){
      #only resample if weights are very variable such that ESS is small
      resample_index = sample(1:N, replace=TRUE, prob=weights_in_time[1,])
      particles_in_time[,1:N] = particles_in_time[,resample_index]
      resample_count = resample_count + 1
    }
  }
  
  cat(100.0*(resample_count / t), "% of timesteps we resample")
  
  out = list(particles_in_time=particles_in_time)
  return(out)
}

N = 10  #number of particles 

#run the SMC for the state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_toy, state_update=toy_state_update, observed_process=observed_process_toy)
particles_in_time = smc_output$particles_in_time

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(particles_in_time, latent_process)
