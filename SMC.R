SMC = function(N, calculate_weight, state_update, observed_process){
  # Function to perform SMC for a state space model. 
  # N is desired number of particles
  # calculate_weight is a function to calculate the weight at each timestep
  # state_update is a function specifying the SSM
  t = length(observed_process)
  particles_in_time = matrix(NA, ncol=N, nrow= t) # N particles at t timesteps
  particle_mean_in_time = rep(NA, t) # N particles at t timesteps
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
      resample_index = sample(1:N, replace=TRUE, prob=weights_in_time[i,])
      particles_in_time[,1:N] = particles_in_time[,resample_index]
      resample_count = resample_count + 1
    }
    
    particle_mean_in_time[i] = sum(particles_in_time[i,])/N
  }
  
  cat(100.0*(resample_count / t), "% of timesteps we resample")
  
  out = list(particles_in_time=particles_in_time, particle_mean_in_time=particle_mean_in_time)
  return(out)
}