library(ggplot2)
library(reshape2)

##################################################################
# Some plotting functions
##################################################################
plot_processes_in_time = function(process_dataframe){
  # Takes a dataframe of 1d processes, e.g. one
  # column each for observed and latent, and plots 
  # as lines in time with labels determined by column 
  # headings in process dataframe
  process_dataframe$time = 1:nrow(process_dataframe)
  melted_data_df = melt(process_dataframe, id.vars='time')
  plot = ggplot(melted_data_df, aes(x=time, y = value, group = variable, colour = variable)) +   geom_line()
  
  return(plot)
  }

plot_particles_and_latent_process = function(particles_in_time, latent_process){
  data_df = data.frame(p=particles_in_time)
  data_df$time = 1:nrow(data_df)
  melted_data_df = melt(data_df, id.vars='time')
  melted_data_df$type = 'particles'
  latent_df = data.frame(l=latent_process)
  latent_df$time = 1:nrow(latent_df)
  melted_latent_df = melt(latent_df, id.vars='time')
  melted_latent_df$type = 'latent process'
  df = rbind(melted_data_df, melted_latent_df)
  head(df)
  ggplot(df, aes(x=time, y = value, group = variable, colour = variable, alpha=0.8)) +   
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(aes(linetype=type, alpha=0.2)) +
    guides(colour=FALSE, alpha=FALSE)
}


trace_plot = function(X, n_iter, title){
    plot = qplot(seq(from= 1, to=n_iter+1, by=1), X ,geom="line", main="",
    xlab="Iteration index", ylab=title, xlim=c(0,n_iter+1), alpha = I(1/1000), size=I(2.5)) +
    theme(axis.title.x = element_text(size = 12), title = element_text(colour='black'),
    axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
    axis.line = element_line(colour="black", size = 1, linetype = "solid"))+
    geom_point(size=2, shape=19, alpha=0.1)
    return(plot)
}


plot_particles_and_latent_process_NOLEGEND = function(particles_in_time, latent_process){
  data_df = data.frame(p=particles_in_time)
  data_df$time = 1:nrow(data_df)
  melted_data_df = melt(data_df, id.vars='time')
  melted_data_df$type = 'particles'
  latent_df = data.frame(l=latent_process)
  latent_df$time = 1:nrow(latent_df)
  melted_latent_df = melt(latent_df, id.vars='time')
  melted_latent_df$type = 'latent process'
  df = rbind(melted_data_df, melted_latent_df)
  head(df)
  ggplot(df, aes(x=time, y = value, group = variable, colour = variable, alpha=0.8)) +   
    scale_linetype_manual(values = c("dashed", "solid")) +
    geom_line(aes(linetype=type, alpha=0.2)) +
    guides(colour=FALSE, alpha=FALSE)+
    theme(legend.position="none", plot.title = element_text(hjust = 0.5))
}


moving_average = function(n_moving, n_iter,X){
  moving_par = rep(NA, n_iter/n_moving)
  moving_par[1]= X$theta_unknown[n_moving]/n_moving
  for (i in seq(2*n_moving, n_iter, by=n_moving)) {
    moving_par[i/n_moving]=(X$moving_par[i]-X$moving_par[i-n_moving])/n_moving
  }
  return(moving_alpha)
}
