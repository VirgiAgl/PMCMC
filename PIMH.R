source("plotting_functions.R")

PIMH = function(n_iter, 
                N, 
                calculate_weight, 
                state_update, 
                observed_process){
  
  #initialisation
  t=length(observed_process)     #number of observed step in times
  n_acceptance=0
  state_values = matrix(NA, nrow = N, ncol = n_iter+1 ) #store the state values at each step
  lik_values = vector(length = n_iter+1 )               #store the lik value at each step
  
  # run an SMC for iteration 0
  smc_output = SMC(N = N, 
                    calculate_weight = calculate_weight, 
                    state_update = state_update, 
                    observed_process = observed_process
                    )
  
  index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ]) #to check
  proposed_x = smc_output$particles_in_time[,index]
  proposed_lik = smc_output$lik_in_time[t]
  
  # store the first two values
  lik_values[1]=proposed_lik
  state_values[,1]=proposed_x
  
  
  for (i in 1:n_iter){
    if (i %% 1000== 0){
      cat("Step ", i, " of ", n_iter, '\n') 
    }
    # run an SMC for each iteration from 1 to n_iter
    smc_output = SMC(N=N, 
                      calculate_weight=calculate_weight, 
                      state_update=state_update, 
                      observed_process=observed_process
                      )
    
    # sample the path x1:xT to consider
    index = sample(1:N, 1, prob = smc_output$weights_in_time[t, ])
    proposed_x = smc_output$particles_in_time[,index]
    proposed_lik = smc_output$lik_in_time[t]
    
    # compute the accepatance probability
    acc_prob = min(c(1, proposed_lik/(lik_values[i]) ))
    
    # accept or reject the new value
    if(runif(1) < acc_prob){
      state_values[, i+1] = proposed_x
      lik_values[i+1] = proposed_lik
      n_acceptance = n_acceptance + 1
    } else {
      state_values[, i+1] = state_values[,i] 
      lik_values[i+1] = lik_values[i]
    }
    
    # compute the ratio of accepted values
    acceptance_ratio = n_acceptance/n_iter
  }
  
  out=list(state_values = state_values,
           lik_values = lik_values,
           acceptance_ratio = acceptance_ratio
           )
  
  return(out)
}

