
n_iter=100                         #number of MCMC iterations 
N=100                               #number of particles 

PIMH = function(n_iter, N, y_values, calculate_weight, state_update, observed_process){
  t=length(y_values)                                    #number of observed step in times
  n_acceptance=0                                        #initialise number of accepted values
  state_values = matrix(NA, nrow = N, ncol = n_iter+1 ) #initialise a matrix to store the state values at each iteration
  lik_values = vector(length = n_iter+1 )               #initialise a matrix to store the lik values at each iteration 
  
  # run an SMC for iteration 0
  smc_output = SMC2(N=N, calculate_weight=calculate_weight, state_update=state_update, observed_process=observed_process)
  index = sample(1:t, 1)
  proposed_x = smc_output$particles_in_time[,index]
  proposed_lik = smc_output$lik_in_time[index]
  
  # store the first two values
  lik_values[1]=proposed_lik
  state_values[,1]=proposed_x
  
  
  for (i in 1:n_iter){

    # run an SMC for each iteration from 1 to n_iter
    smc_output = SMC2(N=N, calculate_weight=calculate_weight, state_update=state_update, observed_process=observed_process)
    # sample the path x1:xT to consider
    index = sample(1:t, 1)
    proposed_x = smc_output$particles_in_time[,index]
    proposed_lik = smc_output$lik_in_time[index]
    lik_values[i+1]=proposed_lik
    
    # compute the accepatance probability
    acc_prob = min(c(1,lik_values[i]/(lik_values[i+1]) ))
    
    # accept or reject the new value
    if(runif(1)<acc_prob){
      state_values[,i+1] = proposed_x
      lik_values[i+1] = proposed_lik
      n_acceptance = n_acceptance + 1
    } else {
      state_values[,i+1] = state_values[,i] 
      lik_values[i+1] = lik_values[i]
    }
    
    # compute the ratio of accepted values
    acceptance_ratio = n_acceptance/n_iter
  }
  
  return(out=list(state_values=state_values,lik_values=lik_values,acceptance_ratio=acceptance_ratio))
}

PIMH_toy = PIMH(n_iter, N, observed_process_toy, calculate_weight=calculate_weight_toy,  state_update=toy_state_update, observed_process=observed_process_toy )

plot1=tracePlot(PIMH_toy$state_values[1,], n_iter, title="Markov Chain for the first particle")
plot1
