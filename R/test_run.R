mod_addint <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001),
                          par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 10000)


mod_addint$sim(cl = cl )

sim_data <- mod_addint$last_sim_output
par_vals <- mod_addint$par_vals

mod_addint$last_sim_output <- sim_data
mod_addint$par_vals <- par_vals

mod_addint$adapt(stat="p1", # model for the type 1 error
               model="binomial",
               type = "as", # Use appropriate model type
               sampling = 150,
               L=1.1,
               d = 1,
               m = 10)

mod_addint$last_stan_fit$profiles()

cl_plot_points <- as.matrix(expand.grid(cl = 6:60, m = 10, cv = 0.2, icc = 0.05))
cl_plot_points_data <- mod_addint$predict(cl_plot_points)
dfp_cl <- data.frame(cl = 6:60,
                     mean = apply(cl_plot_points_data,1,mean),
                     lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                     uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))


ggplot(data = dfp_cl, aes(x = cl, y = mean)) +
  geom_hline(yintercept = 0.05, lty = 1) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "Number of clusters", y = "Type I error", fill = "Iterations") +
  theme_minimal()

mod_addint$sample(1000, type = "entr")

ggplot(data = as.data.frame(mod_addint$par_vals), aes(x = V1, y = V4))+
  geom_density_2d_filled()
