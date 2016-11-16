source("plotting_functions.R")
source("SMC.R")

set.seed(27)

PMMH_linear= function(n_iter, 
                N, 
                d,              # dimension of the unknown parameter space
                calculate_weight, 
                state_update, 
                observed_process, 
                theta_state,     # for the linear model c(mu, phi, sigma)
                theta_obs) {     # for the linear model c(eta)
  
  t=length(observed_process)     #number of observed step in times
  n_acceptance=0
  theta_unknown = matrix(NA, nrow=d, ncol=n_iter+1)     #initialise vector to store the values of the unknown parameters
  state_values = matrix(NA, nrow = t, ncol = n_iter+1 ) #store the state values at each step
  lik_values = vector(length = n_iter+1 )               #store the lik value at each step
  
  for (i in 1:d){                                       #generate the first value for the unknown pars
    theta_unknown[d,1] = rnorm(1,mean=0.95, sd=0)          #to be chosen arbitrarly (assume all unknown pars have the same prior distrubution)
  }
  
  theta_state[2]=theta_unknown[1,1]                     #phi is the second element in the theta_state vector and on the first row in the matrix of unknown pars 
  
  # run an SMC for iteration 0
  smc_output = SMC(N = N, 
                   calculate_weight = calculate_weight, 
                   state_update = state_update, 
                   observed_process = observed_process,
                   theta_state,
                   theta_obs)
  
  index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])  # to check
  proposed_x = smc_output$particles_in_time[,index]
  proposed_lik = smc_output$lik_in_time[t]
  
  # store the first two values
  lik_values[1]=proposed_lik
  state_values[,1]=proposed_x
  
  
  for (i in 1:n_iter){
    cat("iteration in the MCMC = ", i, "\n")
    
                                                    #generate new values for the unknow pars
    theta_unknown[1,i+1] = rnorm(1,mean=theta_unknown[1,i], sd=1)  #proposal updating distribution

    theta_state[2]=theta_unknown[1,i+1]                              #phi is the second element in the theta_state vector and on the first row in the matrix of unknown pars 
    
    
    # run an SMC for each iteration from 1 to n_iter
    smc_output = SMC(N=N, 
                     calculate_weight=calculate_weight, 
                     state_update=state_update, 
                     observed_process=observed_process, 
                     theta_state,
                     theta_obs)
    
    # sample the path x1:xT to consider
    index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])
    proposed_x = smc_output$particles_in_time[,index]
    proposed_lik = smc_output$lik_in_time[t]
    
    # compute the acceptance probability
    first_term = exp(proposed_lik-lik_values[i])
    second_term = dnorm(theta_unknown[1,i+1], mean=1, sd=1)/dnorm(theta_unknown[1,i],mean=1, sd=1)
    third_term = dnorm(theta_unknown[1,i],mean = theta_unknown[1,i+1],sd=1)/dnorm(theta_unknown[1,i+1],mean = theta_unknown[1,i],sd=1)
    acc_prob = min(c(1, first_term*second_term*third_term))
    
    # accept or reject the new value
    if(runif(1) < acc_prob){
      state_values[, i+1] = proposed_x
      lik_values[i+1] = proposed_lik
      n_acceptance = n_acceptance + 1
    } else {
      state_values[, i+1] = state_values[,i] 
      lik_values[i+1] = lik_values[i]
      theta_unknown[1,i+1] = theta_unknown[1,i]
    }
    
    # compute the ratio of accepted values
    acceptance_ratio = n_acceptance/n_iter
  }
  
  out=list(state_values = state_values,
           lik_values = lik_values,
           theta_unknown = theta_unknown,
           acceptance_ratio = acceptance_ratio
  )
  
  return(out)
  
}



set.seed(27)
PMMH_nonlinear = function(n_iter, 
                       N, 
                       d,               # dimension of the unknown parameter space
                       calculate_weight, 
                       state_update, 
                       observed_process, 
                       theta_state,     # for the linear model c(mu, phi, sigma)
                       theta_obs) {     # for the linear model c(eta)
  
  t=length(observed_process)                            #number of observed step in times
  n_acceptance=0
  theta_unknown = matrix(NA, nrow=d, ncol=n_iter+1)     #initialise vector to store the values of the unknown parameters
  state_values = matrix(NA, nrow = t, ncol=n_iter+1 )   #store the state values at each step
  lik_values = vector(length = n_iter+1 )               #store the lik value at each step
  
  for (i in 1:d){                                       #generate the first value for the unknown pars
    theta_unknown[d,1] = rnorm(1,mean=1, sd=1)          #to be chosen arbitrarly (assume all unknown pars have the same prior distrubution)
  }
  
  theta_state[1]=theta_unknown[1,1]                     #phi is the second element in the theta_state vector and on the first row in the matrix of unknown pars 
  
  # run an SMC for iteration 0
  smc_output = SMC(N = N, 
                   calculate_weight = calculate_weight, 
                   state_update = state_update, 
                   observed_process = observed_process_linear,
                   theta_state,
                   theta_obs)
  
  index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])  # to check
  proposed_x = smc_output$particles_in_time[,index]
  proposed_lik = smc_output$lik_in_time[t]
  
  # store the first two values
  lik_values[1]=proposed_lik
  state_values[,1]=proposed_x
  
  
  for (i in 1:n_iter){
    cat("iteration in the MCMC = ", i, "\n")
    
    for (i in 1:d){                                                  #generate new values for the unknow pars
      theta_unknown[d,i+1] = rnorm(1,mean=theta_unknown[d,i], sd=1)  #proposal updating distribution
    }
    
    theta_state[1]=theta_unknown[1,i+1]                              #phi is the second element in the theta_state vector and on the first row in the matrix of unknown pars 
    
    
    # run an SMC for each iteration from 1 to n_iter
    smc_output = SMC(N=N, 
                     calculate_weight=calculate_weight, 
                     state_update=state_update, 
                     observed_process=observed_process_linear, 
                     theta_state,
                     theta_obs)
    
    # sample the path x1:xT to consider
    index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])
    proposed_x = smc_output$particles_in_time[,index]
    proposed_lik = smc_output$lik_in_time[t]
    
    # compute the acceptance probability
    first_term = exp(proposed_lik-lik_values[i])
    second_term = dnorm(theta_unknown[1,i+1], mean=1, sd=1)/dnorm(theta_unknown[1,i],mean=1, sd=1)
    third_term = dnorm(theta_unknown[1,i],mean = theta_unknown[1,i+1],sd=1)/dnorm(theta_unknown[1,i+1],mean = theta_unknown[1,i],sd=1)
    acc_prob = min(c(1, first_term*second_term*third_term))
    
    # accept or reject the new value
    if(runif(1) < acc_prob){
      state_values[, i+1] = proposed_x
      lik_values[i+1] = proposed_lik
      n_acceptance = n_acceptance + 1
    } else {
      state_values[, i+1] = state_values[,i] 
      lik_values[i+1] = lik_values[i]
      theta_unknown[1,i+1] = theta_unknown[1,i]
    }
    
    # compute the ratio of accepted values
    acceptance_ratio = n_acceptance/n_iter
  }
  
  out=list(state_values = state_values,
           lik_values = lik_values,
           phi_values = phi_values,
           theta_unknown = theta_unknown,
           acceptance_ratio = acceptance_ratio
  )
  
  return(out)
  
}
