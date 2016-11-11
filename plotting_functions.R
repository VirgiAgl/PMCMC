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
  ggplot(melted_data_df, aes(x=time, y = value, group = variable, colour = variable)) +   geom_line()
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
  ggplot(df, aes(x=time, y = value, group = variable, colour = variable, linetype=type, alpha=0.8)) +   geom_line() +guides(colour=FALSE, alpha=FALSE)
  
}