
# prior need to store the prior parameters for both the state prior and the parameters prior 
# (mean, var) of the state and then (mean, var) of the par

par=c(phi,mu,sigma,eta) # define a parameter vector
prior=c(0,1.5, 1,1)     # prior means and sds for X_1,logit(phi)  
h=0
LW=F

# phi is unknow, eta sigma mu are known 
sigma_2 = 0.01
par[2]= 0.9           #mu
par[3]= sqrt(sigma_2) #sigma
#eta_2 = (sigma_2/(1-par[1]^2))/5
eta_2 = 0.02
par[4] = sqrt(eta_2)  #eta


SMCpar=function(N,y,prior,par,h=0,LW=F,PLOT=T) { 
  n=length(y)  # number of observations 
  #d=2         # dimension of the state space -- STATE=(X_t,logit(phi)) 
  par_dim=1
  
  # initialise a matrix to store the simulate particles Xt and the unknown parameters
  particles_in_time=matrix(NA, ncol=N, nrow= t)          # N particles at t timesteps
  parameters_in_time=matrix(NA, ncol=par_dim, nrow= t)   # 1 unknow parameter at t timesteps
  
  particle_mean_in_time=matrix(0,nrow=t,ncol=1)          # Mean at each t timestep
  particle_variance_in_time=matrix(0,nrow=t,ncol=1)
  parameter_mean_in_time=matrix(0,nrow=t,ncol=1)          # Mean at each t timestep
  parameter_variance_in_time=matrix(0,nrow=t,ncol=1)
  
  lik_in_time = rep(NA, t)                               # Log likelihood value at each t timestep 
  weights_in_time = matrix(NA, ncol=N, nrow= t)          # N weights at t timesteps
  logl = 0                                               #initialize variable to store likelihood values
  
  # generate the first set of particles 
  particles_in_time[1,] = rnorm(N, mean=prior[1], sd=prior[3])     # Initialise with proposal density
  parameters_in_time[1,] = rnorm(1, mean=prior[2], sd=prior[4])    #for(j in 1:d) particles[,j]=rnorm(N,prior[j],prior[j+d]) 
  
  # calculate weights. This depedends on eta which depends on phi which is unknown and updated at each step
  #eta_2 = (sigma_2/(1-parameters_in_time[1,]^2))/5
  #eta = sqrt(eta_2)
  eta = par[4]
  weight = calculate_weight(observed_val = observed_process[i], particles_vals = particles_in_time[i,])
  
  
  # log of estimate of likelihood 
  logl=logl+log(mean(weight))   #log of the estimate of the likelihood
  lik_in_time[1]=logl           #store the likelihood value
  
  weight_norm = weight/sum(weight)  #normalize the weights
  weights_in_time[1,] = weight_norm #store the weights
  
  # Store the value of the mean and the variance (before resampling) 
  # compute the mean and the variance both for the states and for the parameters
  particle_mean_in_time[1,]=sum(weight_norm*particles_in_time[1,])                                    # expected value of a random var
  particle_variance_in_time[1,]=sum(weight_norm*particles_in_time[1,]^2)-particle_mean_in_time[1,]^2  # variance of a random var
  parameter_mean_in_time[1,] = sum(weight_norm*parameters_in_time[1,])  
  parameter_variance_in_time[1,] = sum(weight_norm*parameters_in_time[1,]^2)-parameter_mean_in_time[1,]^2  # variance of a random var

  
  for(i in 2:n){

    ##(1) Resample (both the states and the parameters) 
    # resample (not for the last step)
    resample_index = sample(1:N, replace=TRUE, prob=weights_in_time[i-1,])
    particles_in_time[,1:N] = particles_in_time[,resample_index]
    parameter_in_time[,1:d] = parameter_in_time[,resample_index]

    
    ##(2) Propagate particles 
    
    #First update parameters 
    if(h>0 && V.st[i-1,2]>0){ ##JITTER 
      if(LW){ #Liu-West -- shrinkage 
        lambda=sqrt(1-h) 
        parameter_in_time[i,]=lambda*parameter_in_time[i-1,]+(1-lambda)*parameter_mean_in_time[i-1,] 
      } 
      #Jitter -- for both LW and non LW 
      parameter_in_time[i,]=parameter_in_time[i-1,]+rnorm(N,0,sqrt(h*parameter_variance_in_time[i-1,])) 
    } 
    
    # update the particles considering that
    # with these updated parameters the prior for the state will change 
    phi=exp(parameter_in_time[i,])/(1+exp(parameter_in_time[i,])) #you update the value of phi to them update the state
    #sigma=gamma*sqrt(1-phi^2)  #rule to update the sigma
    particles_in_time[i,]=(particles_in_time[i-1,]*phi+(par[2]-(phi*par[2])))+rnorm(N,0,par[3]) 
    
    
    ##(3) Weight  
    # check if also the parameter for the update of the state change
    eta = par[4]
    weight = calculate_weight(observed_val = observed_process[i], particles_vals = particles_in_time[i,]) 
    
    ##(4) update log of estimate of likelihood 
    logl=logl+log(mean(weight))   #log of the estimate of the likelihood
    lik_in_time[i,]=logl           #store the likelihood value
    
    
    weight_norm = weight/sum(weight)  #normalize the weights
    weights_in_time[i,] = weight_norm #store the weights
    
    # Store the value of the mean and the variance (before resampling) 
    # compute the mean and the variance both for the states and for the parameters
    particle_mean_in_time[i,]=sum(weight_norm*particles_in_time[i,])                                    # expected value of a random var
    particle_variance_in_time[i,]=sum(weight_norm*particles_in_time[i,]^2)-particle_mean_in_time[i,]^2  # variance of a random var
    parameter_mean_in_time[i,] = sum(weight_norm*parameters_in_time[i,])  
    parameter_variance_in_time[i,] = sum(weight_norm*parameters_in_time[i,]^2)-parameter_mean_in_time[i,]^2  # variance of a random var

  }
  out = list(particles_in_time=particles_in_time, particle_mean_in_time=particle_mean_in_time, particle_variance_in_time=particle_variance_in_time, weights_in_time=weights_in_time, parameters_in_time=parameters_in_time,  parameter_mean_in_time=parameter_mean_in_time, parameter_variance_in_time=parameter_variance_in_time, lik_in_time=lik_in_time)
  return(out)
} 