
# define the noise processes
U_dist = function(x){
  return(dnorm(x, mean= 0, sd=1))
} 

V_dist = function(x){
  return(dnorm(x, mean= 0, sd=1))
} 

# define the state equation
state_update = function(x_k, phi=0.95, mu=0.9, sigma_2= 100){
  sigma = sqrt(sigma_2)
  eta_2 = (sigma_2/(1-phi^2))/5
  eta = sqrt(eta_2)
  error_value = rnorm(1, mean=0, sd=1)
  x_k1 = mu + phi*(x_k-mu) + sigma*error_value
  return (x_k1)  
}

obs_update = function (x_k, N, phi=0.95, mu=0.9, sigma_2= 100){
  sigma = sqrt(sigma_2)
  eta_2 = (sigma_2/(1-phi^2))/5
  eta = sqrt(eta_2)
  error_value = rnorm(N, mean=0, sd=1)
  y_k = x_k + eta*error_value
  return (y_k)
}

# parameters for the algorithm
t = 100 #n_iter 
N = 100  #number of particles 

# get the state X
x_values= vector(length = t)
x1 = rnorm(1, mean=0, sd=1)
x_values[1] = x1

for (i in 2:t){
  x_values[i] = state_update(x_k = x_values[i-1])
}

# get the data Y
y_values= vector(length = t)
for (i in 1:t){
  y_values[i] = obs_update(x_k = x_values[i], N=1) #evaluate the obs eq. once
}

plot(y_values, type= "l")
lines(x_values, col="red")

# define the matrix to store the N particles overtime
matrix_X = matrix(NA, ncol=N, nrow= t)
matrix_W = matrix(NA, ncol=N, nrow= t)


phi=0.95 
mu=0.9
sigma_2= 100
sigma = sqrt(sigma_2)
eta_2 = (sigma_2/(1-phi^2))/5
eta = sqrt(eta_2)

# first step

number_resample = 0

# generate N particles from the chosen initial proposal
x1 = rnorm(N, mean=0, sd=1)
# store the first particles as the first row of matrix x
matrix_X[1,] = x1

# compute and normalize the weights
prob_argument = (y_values[1]-matrix_X[1,])/eta #values at which we want to evaluate the probability function g
weight1  = dnorm(prob_argument, mean=0, sd=1) #in this case the weights are given by the probability density g evaulated at y1|x_1
weight1_norm = weight1/sum(weight1)
matrix_W[1,] = weight1_norm #the weight sum up to one

ESS = sum((matrix_W[1,])^2)^-1 

if (ESS<N/2){
  resample_index = sample(1:N, replace=TRUE, prob=matrix_W[1,])
  matrix_X[,1:N] = matrix_X[,resample_index]
  number_resample = number_resample + 1
}

# second step 
for (i in 2:t){
  matrix_X[i,] = sapply(X =matrix_X[i-1,], FUN = state_update)
  prob_argument = (y_values[i]-matrix_X[i,])/eta
  weight1  = dnorm(prob_argument, mean=0, sd=1)
  weight1_norm = weight1/sum(weight1)
  matrix_W[i,] = weight1_norm
  
  ESS = sum((matrix_W[i,])^2)^-1 
  
  if (ESS<N/2){
    resample_index = sample(1:N, replace=TRUE, prob=matrix_W[1,])
    matrix_X[,1:N] = matrix_X[,resample_index]
    number_resample = number_resample + 1
  }
}


# plot for the particles trajectories over the state space
plot(matrix_X[,1], type="l", ylim = c(min(matrix_X),max(matrix_X)))
for (i in 2:N){
  lines(matrix_X[,i], type="l", col=i+10)
}
lines(x_values, col="red", lty=2)


# kalman filter
#fkf(a0=x_values, as.matrix(1), as.matrix(mu-phi*mu), as.matrix(0), Tt=as.matrix(phi), Zt=as.matrix(1), HHt=as.matrix(sigma), GGt=as.matrix(eta), yt=as.matrix(y_values) )
