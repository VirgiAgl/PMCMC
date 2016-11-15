source("plotting_functions.R")
source("SMC.R")

# we assume the prior distribution for the parameter to be a standard normal (p)
# we assume the updating distribution for the parameter to be a normal distribution 
# centered on the previous state and with sd=1 (q)

# SSM params
mu = 0.9
sigma_2 = 0.01
sigma = sqrt(sigma_2)
eta_2 = 0.02
eta = sqrt(eta_2)
prior_par = c(0,1) # mean and sd of the prior distribution for x1
t = 100

data=generate_data(linear_state_update,linear_obs_update, prior_par, t, plot=TRUE, phi=0.95)
observed_process_linear = data$observed_process



n_acceptance=0
phi_values = rep(NA, n_iter+1)                        # initialise vector to store the values of the unknown parameter
state_values = matrix(NA, nrow = t, ncol = n_iter+1 ) #store the state values at each step
lik_values = vector(length = n_iter+1 )               #store the lik value at each step


phi_values[1] = rnorm(1,mean=1, sd=1)                 # to be chosen arbitrarly 
phi = phi_values[1]

# run an SMC for iteration 0
smc_output = SMC(N = N, 
                 calculate_weight = calculate_weight, 
                 state_update = state_update, 
                 observed_process = observed_process_linear,
                 theta_state = c(mu, phi, sigma),
                 theta_obs = c(eta))

index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ]) #to check
proposed_x = smc_output$particles_in_time[,index]
proposed_lik = smc_output$lik_in_time[t]

# store the first two values
lik_values[1]=proposed_lik
state_values[,1]=proposed_x


for (i in 1:n_iter){
  
  phi = rnorm(1,mean=phi_values[i], sd=1) #proposal updating distribution

  # run an SMC for each iteration from 1 to n_iter
  smc_output = SMC(N=N, 
                   calculate_weight=calculate_weight, 
                   state_update=state_update, 
                   observed_process=observed_process_linear, 
                   theta_state = c(mu, phi, sigma),
                   theta_obs = c(eta))
  
  # sample the path x1:xT to consider
  index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])
  proposed_x = smc_output$particles_in_time[,index]
  proposed_lik = smc_output$lik_in_time[t]
  
  # compute the accepatance probability
  first_term = exp(proposed_lik-lik_values[i])
  second_term = dnorm(phi, mean=1, sd=1)/dnorm(phi_values[i],mean=1, sd=1)
  third_term = dnorm(phi_values[i],mean = phi,sd=1)/dnorm(phi,mean = phi_values[i],sd=1)
  acc_prob = min(c(1, first_term*second_term*third_term))
  
  # accept or reject the new value
  if(runif(1) < acc_prob){
    state_values[, i+1] = proposed_x
    lik_values[i+1] = proposed_lik
    phi_values[i+1] = phi
    n_acceptance = n_acceptance + 1
  } else {
    state_values[, i+1] = state_values[,i] 
    lik_values[i+1] = lik_values[i]
    phi_values[i+1] = phi[i]
  }
  
  # compute the ratio of accepted values
  acceptance_ratio = n_acceptance/n_iter
}

out=list(state_values = state_values,
         lik_values = lik_values,
         phi_values = phi_values,
         acceptance_ratio = acceptance_ratio
)

return(out)
