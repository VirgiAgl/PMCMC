
# This function generates the data independently on the model you want to analyse

source("plotting_functions.R")

generate_data = function (state_update, obs_update, prior_par, t, plot=TRUE, theta_state, theta_obs){
  latent_process = vector(length = t)   # to store process X in time
  observed_process = vector(length = t) # to store observations associated with X
  
  # initialise process X (latent) and Y (observed)
  latent_process[1] = rnorm(1, mean=prior_par[1], sd=prior_par[2]) # store the first set of particles
  observed_process[1] = obs_update(x_k = latent_process[1], theta_obs)
  
  # propagate the state space model
  for (i in 2:t){
    latent_process[i] = state_update(x_k = latent_process[i-1], k = i - 1, theta_state)
    observed_process[i] = obs_update(x_k = latent_process[i], theta_obs)
  }
  
  if (plot==TRUE){
    process_dataframe = data.frame(observed=observed_process, latent=latent_process)
    plot = plot_processes_in_time(process_dataframe)  
  }
  
  return(list(latent_process=latent_process, observed_process=observed_process, plot=plot))
}
