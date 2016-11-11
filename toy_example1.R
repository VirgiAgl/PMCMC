library(ggplot2)
library(reshape2)

##################################################################
# State space model - toy example
##################################################################


state_update = function(x_k, phi, mu, sigma){
  error_value = rnorm(1, mean=0, sd=1)
  x_k1 = mu + phi*(x_k-mu) + sigma*error_value
  return (x_k1)  
}

obs_update = function (x_k, phi, eta){
  error_value = rnorm(1, mean=0, sd=1)
  y_k = x_k + eta*error_value
  return (y_k)
}


##################################################################
# Some plotting functions
##################################################################
plot_processes_in_time = function(process_dataframe){
  # Takes a dataframe of 1d processes, e.g. one
  # column each for observed and latent, and plots 
  # as lines in time with labels determined by column 
  # headings in process dataframe
  data_df$time = 1:nrow(data_df)
  melted_data_df = melt(data_df, id.vars='time')
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
  melted_latent_df$type = 'latent variables'
  df = rbind(melted_data_df, melted_latent_df)
  head(df)
  ggplot(df, aes(x=time, y = value, group = variable, colour = variable, linetype=type, alpha=0.8)) +   geom_line() +guides(colour=FALSE, alpha=FALSE)
  
}


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

# initialise process X (latent) and Y (observed)
latent_process = vector(length = t) # to store process X in time
observed_process= vector(length = t) # to store observations associated with X
x1 = rnorm(1, mean=0, sd=1) # pick starting point from rnorm 
latent_process[1] = x1 # store starting point
observed_process[1] = obs_update(x_k = latent_process[1], phi=phi, eta=eta) # store first observation

# propogate X (latent) and Y (observed) in time by SSM
for (i in 2:t){
  latent_process[i] = state_update(x_k = latent_process[i-1], phi=phi, mu=mu, sigma=sigma)
  observed_process[i] = obs_update(x_k = latent_process[i], phi=phi, eta=eta)
}

# plot the evolution of the process
process_dataframe = data.frame(observed=observed_process, latent=latent_process)
plot_processes_in_time(process_dataframe)


##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################

# Initialise N particles
N = 10  #number of particles 
particles_in_time = matrix(NA, ncol=N, nrow= t) # N particles at t timesteps
weights_in_time = matrix(NA, ncol=N, nrow= t) # N weights at t timesteps
particles_in_time[1,] = rnorm(N, mean=0, sd=1) # Initialise with proposal density

resample_count = 0

for (i in 1:t){
  if (i >= 2) { # All steps other than 1st
    particles_in_time[i,] = sapply(X = particles_in_time[i-1,], FUN = state_update, phi=phi, mu=mu, sigma=sigma)
  }
  
  # This weight calculation is SSM dependant
  prob_argument = (observed_process[i]-particles_in_time[i,])/eta
  weight  = dnorm(prob_argument, mean=0, sd=1) # weights here are from prob density g evaulated at y1|x_1 for all steps
  
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


##################################################################
# Plot the particles trajectories in path space
##################################################################

# plot for the particles trajectories over the state space
plot_particles_and_latent_process(particles_in_time, latent_process)


##################################################################
# Run Kalman filter for comparison of linear guassian SSMs
##################################################################

# kalman filter
#fkf(a0=latent_process, as.matrix(1), as.matrix(mu-phi*mu), as.matrix(0), Tt=as.matrix(phi), Zt=as.matrix(1), HHt=as.matrix(sigma), GGt=as.matrix(eta), yt=as.matrix(observed_process) )
