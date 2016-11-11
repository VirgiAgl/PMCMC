
# define the state equation
state_update = function(x_k, phi=0.95, mu=0.9, sigma_2= 0.5){
  sigma = sqrt(sigma_2)
  eta_2 = (sigma_2/(1-phi^2))/5
  eta = sqrt(eta_2)
  error_value = rnorm(1, mean=0, sd=1)
  x_k1 = mu + phi*(x_k-mu) + sigma*error_value
  return (x_k1)  
}

# define the obs equation
obs_update = function (x_k, N, phi=0.95, mu=0.9, sigma_2= 0.5){
  sigma = sqrt(sigma_2)
  eta_2 = (sigma_2/(1-phi^2))/5
  eta = sqrt(eta_2)
  error_value = rnorm(N, mean=0, sd=1)
  y_k = x_k + eta*error_value
  return (y_k)
}Q

# iteration steps of the state space model 
t = 20 

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


# Parameters of the Kalman Filter
phi=0.95 
mu=0.9
sigma_2= 0.5
sigma = sqrt(sigma_2)
eta_2 = (sigma_2/(1-phi^2))/5
eta = sqrt(eta_2)




# Traslate the model in state space form
toy_model_state_space <- function(phi, mu, sigma,y_values, x_values) {
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
sp = toy_model_state_space(phi=phi, mu=mu, sigma=sigma ,y_values=y_values, x_values=x_values)
ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = sp$yt)

# Store the values for the plot the predictive distributions
predicted_exp_value = t(cbind(ans$at[1], ans$att))
predicted_var = as.matrix(c(ans$Pt[1], as.vector(ans$Ptt)))
x = seq(-6,6, by=0.01)
values_to_plot = matrix(NA, nrow = length(x),ncol=length(predicted_exp_value))
for (i in 1:(t+1)){
  values_to_plot[,i] = dnorm(x, mean=predicted_exp_value[i,],sd=sqrt(predicted_var[i,]))
}
data_frame = as.data.frame(values_to_plot)
data_frame$x_values = x
data_frame_melt = melt(data_frame, id=c("x_values"))
n_values=as.matrix(table(data_frame_melt$variable))


# Plot the sequence of predictive distributions
ggplot(data_frame_melt, aes(x=value, y = x_values, group = variable, colour = "black")) +  geom_point(size=0.05, colour="#000099") + ggtitle("Predictive densities") +
  labs(x="Time index",y="State values") + facet_wrap(~variable, ncol = n_values[1] )  + theme(strip.background = element_blank(), strip.text.x = element_blank())

