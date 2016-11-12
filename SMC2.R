SMC = function(N, calculate_weight, state_update, observed_process){
  # Function to perform SMC for a state space model. 
  # N is desired number of particles
  # calculate_weight is a function to calculate the weight at each timestep
  # state_update is a function specifying the SSM
  
  t = length(observed_process)
  d=1         ##dimension of the state space (in this case is one because we don't have parameters to update)
  
  particles_in_time = matrix(NA, ncol=N, nrow= t) # N particles at t timesteps
  particle_mean_in_time=matrix(0,nrow=t,ncol=d)   # N particles at t timesteps
  particle_variance_in_time=matrix(0,nrow=t,ncol=d)
  lik_in_time = rep(NA, t)                        # Log likelihood value at eah t timestep 
  weights_in_time = matrix(NA, ncol=N, nrow= t)   # N weights at t timesteps
  particles_in_time[1,] = rnorm(N, mean=0, sd=1)  # Initialise with proposal density
  
  
  for (i in 1:t){
    if (i >= 2) { # All steps other than 1st
      particles_in_time[i,] = sapply(X = particles_in_time[i-1,], FUN = state_update, k=i)#phi=phi, mu=mu, sigma=sigma)
    }
    
    weight = calculate_weight(observed_val = observed_process[i], particles_vals = particles_in_time[i,])
    
    if (i = 1){
      logl=logl+log(mean(w))
    }
    logl=logl+log(mean(w)) 
    logl=log(mean(weight)) #log of estimate of likelihood 
    lik_in_time[1]=logl 
    
    weight_norm = weight/sum(weight)
    weights_in_time[i,] = weight_norm
    
    ##STORAGE OF FEATURES OF FILTER -- MEAN AND VAR 
    particle_mean_in_time=matrix(0,nrow=t,ncol=d) 
    particle_variance_in_time=matrix(0,nrow=t,ncol=d) 
    particle_mean_in_time[i,]=sum(weight_norm*particles_in_time[i,]) # expected value of a random var
    particle_variance_in_time[i,]=sum(weight_norm*particles_in_time[i,]^2)-particle_mean_in_time[i,]^2  # variance of a random var
    
    resample_index = sample(1:N, replace=TRUE, prob=weights_in_time[i,])
    particles_in_time[,1:N] = particles_in_time[,resample_index]
  }
  
  cat(100.0*(resample_count / t), "% of timesteps we resample")
  
  out = list(particles_in_time=particles_in_time, particle_mean_in_time=particle_mean_in_time)
  return(out)
}