
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

# initialise process X (unobserved) and Y (observed)
x_values = vector(length = t) # to store process X in time
y_values= vector(length = t) # to store observations associated with X
x1 = rnorm(1, mean=0, sd=1) # pick starting point from rnorm 
x_values[1] = x1 # store starting point
y_values[1] = obs_update(x_k = x_values[1], phi=phi, eta=eta) # store first observation

# propogate X  and Y in time by SSM
for (i in 2:t){
  x_values[i] = state_update(x_k = x_values[i-1], phi=phi, mu=mu, sigma=sigma)
  y_values[i] = obs_update(x_k = x_values[i], phi=phi, eta=eta)
}

# plot the evolution of the process
plot(y_values, type= "l")
lines(x_values, col="red")


##################################################################
# Use SMC/particle filter to sample from underlying process X
##################################################################

# Initialise N particles
N = 10  #number of particles 
matrix_X = matrix(NA, ncol=N, nrow= t) # N particles at t timesteps
matrix_W = matrix(NA, ncol=N, nrow= t) # N weights at t timesteps
matrix_X[1,] = rnorm(N, mean=0, sd=1) # Initialise with proposal density

resample_count = 0

for (i in 1:t){
  if (i >= 2) { # All steps other than 1st
    matrix_X[i,] = sapply(X = matrix_X[i-1,], FUN = state_update, phi=phi, mu=mu, sigma=sigma)
  }
  
  prob_argument = (y_values[i]-matrix_X[i,])/eta
  weight  = dnorm(prob_argument, mean=0, sd=1) # weights here are from prob density g evaulated at y1|x_1 for all steps
  weight_norm = weight/sum(weight)
  matrix_W[i,] = weight_norm
  
  ESS = sum((matrix_W[i,])^2)^-1 
  
  if (ESS<N/2){
    resample_index = sample(1:N, replace=TRUE, prob=matrix_W[1,])
    matrix_X[,1:N] = matrix_X[,resample_index]
    resample_count = resample_count + 1
  }
}

cat(100.0*(resample_count / t), "% of timesteps we resample")


##################################################################
# Plot the particles trajectories in path space
##################################################################

# plot for the particles trajectories over the state space
plot(matrix_X[,1], type="l", ylim = c(min(matrix_X),max(matrix_X)))
for (i in 2:N){
  lines(matrix_X[,i], type="l", col=i+10)
}
lines(x_values, col="red", lty=2)


##################################################################
# Run Kalman filter for comparison of linear guassian SSMs
##################################################################

# kalman filter
#fkf(a0=x_values, as.matrix(1), as.matrix(mu-phi*mu), as.matrix(0), Tt=as.matrix(phi), Zt=as.matrix(1), HHt=as.matrix(sigma), GGt=as.matrix(eta), yt=as.matrix(y_values) )
