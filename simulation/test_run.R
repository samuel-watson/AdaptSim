mod_addint <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001),
                          par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 10000)

# IV

mod_addint <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(0,0,1,0,0,10), 
                          par_upper = c(0.5,0.8,20,1,1,1000), par_discrete = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE), n = 10000)



mod_addint$sim(cl = cl )

sim_data <- mod_addint$last_sim_output
par_vals <- mod_addint$par_vals

mod_addint$last_sim_output <- sim_data
mod_addint$par_vals <- par_vals

mod_addint$adapt(stat="beta", # model for the type 1 error
               model="linear",
               type = "addint", # Use appropriate model type
               sampling = 150,
               L=1.1,
               d = 1,
               m = 10)

mod_addint$last_stan_fit$profiles()

cl_plot_points <- as.matrix(expand.grid(gamma = 0.1,rho = 0.5,n_z = 1:20, b = 0.15, lambda = 0.5, n = 300))
cl_plot_points_data <- mod_addint$predict(cl_plot_points)
dfp_cl <- data.frame(n_z = 1:20,
                     mean = apply(cl_plot_points_data,1,mean),
                     lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                     uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))


ggplot(data = dfp_cl, aes(x = n_z, y = mean)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "Number of instruments", y = "Type I error", fill = "Iterations") +
  theme_minimal()

cl_plot_points <- as.matrix(expand.grid(gamma = 0.1,rho = 0.5,n_z = 2, b = seq(0,1,length.out = 20), lambda = 0.5, n = 300))
cl_plot_points_data <- mod_addint$predict(cl_plot_points)
dfp_cl <- data.frame(b = seq(0,1,length.out = 20),
                     mean = apply(cl_plot_points_data,1,mean),
                     lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                     uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))


ggplot(data = dfp_cl, aes(x = b, y = mean)) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "b", y = "Type I error", fill = "Iterations") +
  theme_minimal()

cl_plot_points <- as.matrix(expand.grid(gamma = 0.01,rho = 0.8,n_z = 1:20, b =0.05, lambda = seq(0,0.5,length.out = 20), n = 300))
cl_plot_points_data <- mod_addint$predict(cl_plot_points)
dfp_cl <- data.frame(l = cl_plot_points[,5],
                     n_z = cl_plot_points[,3],
                     mean = apply(cl_plot_points_data,1,mean),
                     lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                     uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))


ggplot(data = dfp_cl, aes(x = l, y = n_z, fill = mean)) +
  geom_tile() +
  labs(x = "lambda", y = "Type I error", fill = "Iterations") +
  theme_minimal()+
  scale_fill_viridis_c()

mod_addint$sample(10000, type = "var")

ggplot(data = as.data.frame(mod_addint$par_vals), aes(x = V3, y = V5))+
  geom_density_2d_filled()
