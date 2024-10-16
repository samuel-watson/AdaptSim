# WEAK INSTRUMENTS VERSION

require(glmmrBase)
require(pbapply)
if(!require(AdaptSim)){
  devtools::install_github("samuel-watson/AdaptSim") # note you will need the package devtools for this
  require(AdaptSim)
}

data_sim <- function(x){
  Z <- matrix(rnorm(x[6]*x[3]),nrow = x[6], ncol = x[3])
  L <- matrix(c(1,x[2],0,sqrt(1 - x[2]^2)),nrow = 2)
  u <- t(L%*%t(matrix(rnorm(x[6]*2),ncol = 2)))
  pi <- x[4]*exp(-c(0:(x[3]-1))*x[5])
  D <- Z%*%pi + u[,2]
  gamma <- rnorm(x[3],0,x[1])
  Y <- Z%*%gamma + u[,1]
  df <- data.frame(y = Y, d = D)
  df <- cbind(df, as.data.frame(Z))
  colnames(df)[3:ncol(df)] <- paste0("z",1:x[3])
  return(df)
}



fit_model <- function(data=data){
  fit1 <- ivreg::ivreg(formula = paste0("y ~ 1 | d | ", paste0(paste0("z",1:(ncol(data)-2)),collapse = " + ")), data = data)
  
  if(is(fit1,"ivreg")){
    s1 <- summary(fit1)
    tstat <- s1$coefficients[2,3]
    pval <- 2*(1-pnorm(abs(tstat)))
    
    return(c(beta = s1$coefficients[2,1],
             se = s1$coefficients[2,2],
             p = pval,
             p1 = I(pval < 0.05)*1,
             cover = I((s1$coefficients[2,1] + qnorm(0.975)*s1$coefficients[2,2]) > 0 &&
                         (s1$coefficients[2,1] - qnorm(0.975)*s1$coefficients[2,2]) < 0)*1))
  } else {
    # return NA if the model fails to fit
    return(c(beta = NA,
             se = NA,
             p = NA,
             p1 = NA,
             cover = NA))
  }
}

x <- c(0.1,0.5,4,0.5,0.02,100)
df <- data_sim(x)
fit_model(df)

# Define the 15 points provided by the instructor

standard_sim <- function(x)fit_model(data_sim(x))

if(file.exists(paste0(getwd(),"/results/iv_sim_example.RDS"))){
  # deliberately chosen values with maximal spread
  
  # sample 100 points
  vals <- matrix(
    c(runif(100,0,0.5),
      runif(100,0,0.8),
      sample(1:20,100,replace = TRUE),
      runif(100,0,1),
      runif(100,0,0.5),
      sample(10:1000,100,replace = TRUE)),
    ncol = 6
  )
  
  cl <- parallel::makeCluster(parallel::detectCores()-2)
  parallel::clusterEvalQ(cl,.libPaths("C:/R/lib"))
  parallel::clusterEvalQ(cl,library(ivreg))
  parallel::clusterExport(cl, list("data_sim", "fit_model", "standard_sim", "vals"))
  
  # Function to run simulations for a given point
  run_simulation_for_point <- function(index) {
    point <- vals[index, ]
    pbapply::pbreplicate(10000, standard_sim(point), cl = cl)
  }
  
  # Run simulations for the provided points
  results <- list()
  for (i in 1:nrow(vals)) {
    results[[i]] <- run_simulation_for_point(i)
  }
  
  saveRDS(vals, paste0(getwd(),"/results/iv_vals.RDS"))
  saveRDS(results, paste0(getwd(),"/results/iv_sim_example.RDS")) # change the path as required
  parallel::stopCluster(cl)
} else {
  results <- readRDS(paste0(getwd(),"/results/iv_sim_example.RDS"))
  vals <- readRDS(paste0(getwd(),"/results/iv_vals.RDS"))
}

# Compute true values for the provided points
true_values <- sapply(results, function(res) mean(res["beta", ]))

mod_additive <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(0,0,1,0,0,10), par_upper = c(0.5,0.8,20,1,0.5,1000), par_discrete = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE), n = 10000)
mod_addint <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(0,0,1,0,0,10), par_upper = c(0.5,0.8,20,1,0.5,1000), par_discrete = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE), n = 1000)
mod_linear <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(0,0,1,0,0,10), par_upper = c(0.5,0.8,20,1,0.5,1000), par_discrete = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE), n = 10000)
mod_as <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(0,0,1,0,0,10), par_upper = c(0.5,0.8,20,1,0.5,1000), par_discrete = c(FALSE,FALSE,TRUE,FALSE,FALSE,TRUE), n = 10000)

cl <- parallel::makeCluster(4) # you can change the number of cores
parallel::clusterEvalQ(cl,.libPaths("C:/R/lib"))
parallel::clusterEvalQ(cl,library(ivreg))
parallel::clusterExport(cl,c("data_sim","fit_model"))

#########################################
# code block
#########################################

block_size <- 10000
total_iterations <- 100000

results_list <- list()  # Only include the models you want to run
data_list <- list() ## SAM: we can save the simulated data in a list to recall it for each model
par_list <- list() ## SAM: we also need to save the parameter values
blocks <- seq(block_size, total_iterations, by = block_size) ## SAM: I've also moved this outside the loop to help with indexing the data
gen_new_points <- TRUE ## SAM: Flag for if we need new points

for (mod_type_name in c("mod_additive","mod_addint","mod_linear","mod_as")) {  # Only include the models you want to run
  mod_type <- get(mod_type_name)
  
  for (block in blocks) {
    cat("Running block of", block, "iterations for", mod_type_name, "\n")
    
    if(length(data_list) < which(blocks == block)){
      tryCatch({
        mod_type$sim(cl = cl)
      }, error = function(e) {
        message("Error during parallel execution: ", e)
        stop(e)
      })
      data_list[[length(data_list)+1]] <- mod_type$last_sim_output
      par_list[[length(par_list)+1]] <- mod_type$par_vals
      gen_new_points <- TRUE
    } else {
      mod_type$last_sim_output <- data_list[[which(blocks == block)]]
      mod_type$par_vals <- par_list[[which(blocks == block)]]
      gen_new_points <- FALSE
    }
    
    # Setting adaptation model parameters
    model_type <- switch(mod_type_name,
                         mod_additive = "additive",
                         mod_addint = "addint",
                         mod_linear = "as",
                         mod_as = "as")
    
    mod_type$adapt(stat="p1", # model for the type 1 error
                   model="linear",
                   type = model_type, # Use appropriate model type
                   sampling = 150,
                   L=1.1,
                   d = ifelse(mod_type_name == "mod_linear",1,2),
                   m = 10)
    
    pred_points <- vals
    pred_points_data <- mod_type$predict(pred_points)
    
    rmse <- c()
    for(i in 1:nrow(pred_points_data)) rmse <- c(rmse, sqrt(mean((pred_points_data[i, ] - true_values[i])^2)))
    
    cl_plot_points <- as.matrix(expand.grid(gamma = 0.1,rho = 0.5,n_z = 1:20, b = 0.5, lambda = 0.2, n = 300))
    cv_plot_points <- as.matrix(expand.grid(gamma = 0.1,rho = 0.5,n_z = 2, b = seq(0,1,length.out = 20), lambda = 0.2, n = 300))
    icc_plot_points <- as.matrix(expand.grid(gamma = seq(0,0.5,length.out = 20),rho = 0.5,n_z = 2, b = 0.5, lambda = 0.2, n = 300))
    
    cl_plot_points_data <- mod_type$predict(cl_plot_points)
    cv_plot_points_data <- mod_type$predict(cv_plot_points)
    icc_plot_points_data <- mod_type$predict(icc_plot_points)
    
    dfp_cl <- data.frame(n_z = 6:60,
                         mean = apply(cl_plot_points_data,1,mean),
                         lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                         uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))
    
    dfp_cv <- data.frame(b = seq(0,1,length.out = 20),
                         mean = apply(cv_plot_points_data,1,mean),
                         lci = apply(cv_plot_points_data,1,function(i)quantile(i,0.025)),
                         uci = apply(cv_plot_points_data,1, function(i)quantile(i,0.975)))
    
    dfp_icc <- data.frame(gamma = seq(0,0.5,length.out = 20),
                          mean = apply(icc_plot_points_data,1,mean),
                          lci = apply(icc_plot_points_data,1,function(i)quantile(i,0.025)),
                          uci = apply(icc_plot_points_data,1, function(i)quantile(i,0.975)))
    
    results_list[[mod_type_name]][[which(blocks == block)]] <- list(pred_points_data = pred_points_data,
                                                                    dfp_cl = dfp_cl,
                                                                    dfp_cv = dfp_cv,
                                                                    dfp_icc = dfp_icc,
                                                                    rmse = rmses)
    if(gen_new_points){
      tryCatch({
        mod_type$sample(block_size, type = "none")
      }, error = function(e) {
        message("Error during model sampling: ", e)
        next
      })
    }
    
    saveRDS(results_list, paste0(getwd(),"/results/iv_results_list.RDS"))
    saveRDS(data_list, paste0(getwd(),"/results/iv_data_list.RDS"))
    saveRDS(par_list, paste0(getwd(),"/results/iv_par_list.RDS"))
  }
}


#########################
# SUMMARIZE THE OUTPUT  #
#########################

require(ggplot2)
require(gridExtra)

# Create a list to store plots
plot_list <- list()

# Only include the models you want to run
for (mod_type_name in c("mod_additive","mod_linear")) {
  dfp_all_cl <- data.frame()
  dfp_all_cv <- data.frame()
  dfp_all_icc <- data.frame()
  
  for (block in blocks) {
    dfp_cl <- results_list[[mod_type_name]][[which(blocks == block)]]$dfp_cl
    dfp_cl$iterations <- block
    dfp_cl$model <- mod_type_name
    dfp_all_cl <- rbind(dfp_all_cl, dfp_cl)
    
    dfp_cv <- results_list[[mod_type_name]][[which(blocks == block)]]$dfp_cv
    dfp_cv$iterations <- block
    dfp_cv$model <- mod_type_name
    dfp_all_cv <- rbind(dfp_all_cv, dfp_cv)
    
    dfp_icc <- results_list[[mod_type_name]][[which(blocks == block)]]$dfp_icc
    dfp_icc$iterations <- block
    dfp_icc$model <- mod_type_name
    dfp_all_icc <- rbind(dfp_all_icc, dfp_icc)
  }
  
  # Create the plot for the current model (Number of clusters)
  p_cl <- ggplot(data = dfp_all_cl, aes(x = cl, y = mean, fill = as.factor(iterations))) +
    geom_hline(yintercept = 0.05, lty = 1) +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_line() +
    labs(x = "Number of clusters", y = "Type I error", fill = "Iterations") +
    ggtitle(paste("Model Type:", mod_type_name, "- Clusters")) +
    theme_minimal()
  
  # Create the plot for the current model (CV)
  p_cv <- ggplot(data = dfp_all_cv, aes(x = cv, y = mean, fill = as.factor(iterations))) +
    geom_hline(yintercept = 0.05, lty = 1) +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_line() +
    labs(x = "Coefficient of Variation", y = "Type I error", fill = "Iterations") +
    ggtitle(paste("Model Type:", mod_type_name, "- CV")) +
    theme_minimal()
  
  # Create the plot for the current model (ICC)
  p_icc <- ggplot(data = dfp_all_icc, aes(x = icc, y = mean, fill = as.factor(iterations))) +
    geom_hline(yintercept = 0.05, lty = 1) +
    geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
    geom_line() +
    labs(x = "Intraclass Correlation Coefficient", y = "Type I error", fill = "Iterations") +
    ggtitle(paste("Model Type:", mod_type_name, "- ICC")) +
    theme_minimal()
  
  # Add the plots to the list
  plot_list[[paste(mod_type_name, "clusters")]] <- p_cl
  plot_list[[paste(mod_type_name, "cv")]] <- p_cv
  plot_list[[paste(mod_type_name, "icc")]] <- p_icc
}

# Display all plots together
do.call(grid.arrange, c(plot_list, ncol = 2))

## as an alternative
require(patchwork)

(plot_list$`mod_additive clusters` + plot_list$`mod_linear clusters`) / (plot_list$`mod_additive icc` + plot_list$`mod_linear icc`) / (plot_list$`mod_additive cv` + plot_list$`mod_linear cv`) + plot_layout(guides = "collect")


# Optionally, save the plots
for (plot_name in names(plot_list)) {
  ggsave(filename = paste0("plot_", plot_name, ".png"), plot = plot_list[[plot_name]], width = 8, height = 6)
}

# RMSE Plot
rmse_df <- data.frame()

for (mod_type_name in c("mod_additive", "mod_linear")) {
  for (block in blocks) {
    rmse_data <- results_list[[mod_type_name]][[which(blocks == block)]]$rmse
    rmse_df <- rbind(rmse_df, data.frame(iterations = block, rmse = rmse_data, model = mod_type_name))
  }
}

rmse_plot <- ggplot(rmse_df, aes(x = iterations, y = rmse, color = model, group = interaction(model, rownames(rmse_df)))) +
  geom_line() +
  geom_point() +
  labs(x = "Iterations", y = "RMSE") +
  theme_minimal()

ggsave(filename = "plot_rmse.png", plot = rmse_plot, width = 8, height = 6)




p_cl <- ggplot(data = results_list$mod_additive[[10]]$dfp_cl, aes(x = cl, y = mean, fill = as.factor(iterations))) +
  geom_hline(yintercept = 0.05, lty = 1) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "Number of clusters", y = "Type I error", fill = "Iterations") +
  ggtitle(paste("Model Type:", mod_type_name, "- Clusters")) +
  theme_minimal()

p_cl <- ggplot(data = dfp_all_cl[dfp_all_cl$iterations > 20000,], aes(x = cl, y = mean, fill = as.factor(iterations))) +
  geom_hline(yintercept = 0.05, lty = 1) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "Number of clusters", y = "Type I error", fill = "Iterations") +
  ggtitle(paste("Model Type:", mod_type_name, "- Clusters")) +
  theme_minimal()
