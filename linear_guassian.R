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
sigma_2 = 0.01
sigma = sqrt(sigma_2)
eta_2 = (sigma_2/(1-phi^2))/5
eta = sqrt(eta_2)

toy_state_update = function(x_k, k){
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
  weight  = dnorm(observed_val, mean = particles_vals, sd=eta) # weights here are from prob density g evaulated at y1|x_1 for all steps
}

##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################

N = 100  #number of particles 

#run the SMC for the state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_toy, state_update=toy_state_update, observed_process=observed_process_toy)
particles_in_time = smc_output$particles_in_time
particle_mean_in_time = smc_output$particle_mean_in_time

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(particles_in_time, latent_process_toy)

##################################################################
# Kalman filter for predictive distribution
##################################################################
library(FKF)

toy_model_state_space <- function(phi, mu, sigma, y_values, x_values) {
  Tt = matrix(phi)
  Zt = matrix(1)
  dt = as.matrix(mu-phi*mu)
  ct = as.matrix(0)
  HHt = matrix(sigma)
  GGt = matrix(mu)
  yt = as.matrix(y_values)
  a0 <- x_values[1]
  P0 <- matrix(1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt, HHt=HHt, yt=rbind(y_values)))
}

# Compute the Kalman filter       
sp = toy_model_state_space(phi=phi, mu=mu, sigma=sigma ,y_values=observed_process_toy[1:20], x_values=latent_process_toy[1:20])
ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = sp$yt)


# Store the values for the plot the predictive distributions
predicted_exp_value = t(cbind(ans$at[1], ans$att))
predicted_var = as.matrix(c(ans$Pt[1], as.vector(ans$Ptt)))
x = seq(-8,8, by=0.01)
values_to_plot = matrix(NA, nrow = length(x),ncol=length(predicted_exp_value))
for (i in 1:(20+1)){
  values_to_plot[,i] = dnorm(x, mean=predicted_exp_value[i,],sd=sqrt(predicted_var[i,]))
}
data_frame = as.data.frame(values_to_plot)
data_frame$x_values = x
data_frame_melt = melt(data_frame, id=c("x_values"))
n_values=as.matrix(table(data_frame_melt$variable))

particle_means_df = data.frame(means = particle_mean_in_time[1:21])
particle_means_df$variable = colnames(data_frame)[1:21]

data_frame_melt$variable = factor(data_frame_melt$variable, levels=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V20','V21'))
particle_means_df$variable = factor(particle_means_df$variable, levels = c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V20','V21'))

# Plot the sequence of predictive distributions
plot = ggplot() +  
  geom_point(data=data_frame_melt, aes(x=value, y = x_values, group = variable, colour = "black"), size=0.05, colour="#000099") + 
  geom_point(data = particle_means_df, aes(x=0, y = means, colour = "red", group = variable, size=20)) +
  ggtitle("Predictive densities") +
  labs(x="Time index",y="State values") + 
  facet_wrap(~variable, ncol = n_values[1] )+
  theme(strip.background = element_blank(), strip.text.x = element_blank())
plot##################################################################
# State space model - toy example
##################################################################
##################################################################
# Generate some observations using the SSM
##################################################################

t = 100 # num of iters

# SSM params
phi = 0.95
mu = 0.9
sigma_2 = 0.01
sigma = sqrt(sigma_2)
eta_2 = (sigma_2/(1-phi^2))/5
eta = sqrt(eta_2)

toy_state_update = function(x_k, k){
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

N = 100  #number of particles 

#run the SMC for the state space model
smc_output = SMC(N=N, calculate_weight=calculate_weight_toy, state_update=toy_state_update, observed_process=observed_process_toy)
particles_in_time = smc_output$particles_in_time
particle_mean_in_time = smc_output$particle_mean_in_time
final_weights = smc_output$final_weights
particle_index = sample(1:N, size = 1, prob = final_weights)
sample_particle = particles_in_time[, particle_index]

# plot for the particles trajectories over the state space, along with the actual latent process used to generate the data we train on
plot_particles_and_latent_process(particles_in_time, latent_process_toy)

##################################################################
# Kalman filter for predictive distribution
##################################################################

toy_model_state_space <- function(phi, mu, sigma, y_values, x_values) {
  Tt = matrix(phi)
  Zt = matrix(1)
  dt = as.matrix(mu-phi*mu)
  ct = as.matrix(0)
  HHt = matrix(sigma)
  GGt = matrix(mu)
  yt = as.matrix(y_values)
  a0 <- x_values[1]
  P0 <- matrix(1)
  return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt, HHt=HHt, yt=rbind(y_values)))
}

# Compute the Kalman filter       
sp = toy_model_state_space(phi=phi, mu=mu, sigma=sigma ,y_values=observed_process_toy[1:20], x_values=latent_process_toy[1:20])
ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = sp$yt)


# Store the values for the plot the predictive distributions
predicted_exp_value = t(cbind(ans$at[1], ans$att))
predicted_var = as.matrix(c(ans$Pt[1], as.vector(ans$Ptt)))
x = seq(-8,8, by=0.01)
values_to_plot = matrix(NA, nrow = length(x),ncol=length(predicted_exp_value))
for (i in 1:(20+1)){
  values_to_plot[,i] = dnorm(x, mean=predicted_exp_value[i,],sd=sqrt(predicted_var[i,]))
}
data_frame = as.data.frame(values_to_plot)
data_frame$x_values = x
data_frame_melt = melt(data_frame, id=c("x_values"))
n_values=as.matrix(table(data_frame_melt$variable))

particle_means_df = data.frame(means = particle_mean_in_time[1:21])
particle_means_df$variable = colnames(data_frame)[1:21]

data_frame_melt$variable = factor(data_frame_melt$variable, levels=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V20','V21'))
particle_means_df$variable = factor(particle_means_df$variable, levels = c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14','V15','V16','V17','V18','V19','V20','V21'))

# Plot the sequence of predictive distributions
plot = ggplot() +  
  geom_point(data=data_frame_melt, aes(x=value, y = x_values, group = variable, colour = "black"), size=0.05, colour="#000099") + 
  geom_point(data = particle_means_df, aes(x=0, y = means, colour = "red", group = variable, size=20)) +
  ggtitle("Predictive densities") +
  labs(x="Time index",y="State values") + 
  facet_wrap(~variable, ncol = n_values[1] )+
  theme(strip.background = element_blank(), strip.text.x = element_blank())
plot
