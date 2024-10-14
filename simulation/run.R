require(glmmrBase)
require(pbapply)
if(!require(AdaptSim)){
  devtools::install_github("samuel-watson/AdaptSim") # note you will need the package devtools for this
  require(AdaptSim)
}

data_sim <- function(x){
  df <- data.frame(cl = 1:x[1])

  sizes <- rep(0,x[1])
  while(any(sizes<=0)){
    sizes[sizes <= 0] <- round(rnorm(length(sizes[sizes <= 0]), x[2], x[3]*x[2]),0)
  }

  df <- df[rep(1:x[1],sizes),,drop=FALSE]
  df$int <- 0
  df[df$cl > x[1]/2,'int'] <- 1

  mod <- Model$new(
    formula = ~ int + (1|gr(cl)),
    covariance = c(sqrt(x[4])),
    mean = rep(0,2),
    data = df,
    family = gaussian(),
    var_par = sqrt(1-x[4])
  )

  df$y <- mod$sim_data()
  return(df)
}

fit_model <- function(data=data){
  fit1 <- lme4::lmer(y~int+(1|cl),data=data, REML= FALSE)#tryCatch(suppressWarnings(lme4::lmer(y~int+(1|cl),data=data)),error=function(i)NA)

  if(is(fit1,"lmerMod")){
    s1 <- summary(fit1)
    tstat <- s1$coefficients[2,3]
    pval <- 2*(1-pnorm(abs(tstat)))

    return(c(conv = fit1@optinfo$conv$opt,
             beta = s1$coefficients[2,1],
             se = s1$coefficients[2,2],
             p = pval,
             p1 = I(pval < 0.05)*1,
             cover = I((s1$coefficients[2,1] + qnorm(0.975)*s1$coefficients[2,2]) > 0 &&
                         (s1$coefficients[2,1] - qnorm(0.975)*s1$coefficients[2,2]) < 0)*1))
  } else {
    # return NA if the model fails to fit
    return(c(conv = 1,
             beta = NA,
             se = NA,
             p = NA,
             p1 = NA,
             cover = NA))
  }
}

df <- data_sim(c(10,10,0.4,0.05))
fit_model(df)

# Define the 15 points provided by the instructor
vals <- matrix(c(
  34, 199, 0.23637690, 0.09193269,
  45, 130, 0.11546262, 0.19733145,
  13,  38, 0.00000000, 0.23485321,
  17,  28, 0.25326985, 0.03161119,
  6, 176, 0.90833488, 0.00100000,
  20,   5, 1.00000000, 0.40000000,
  8, 200, 0.23894541, 0.14370243,
  20,  33, 1.00000000, 0.10095784,
  60,   5, 0.98865012, 0.39488326,
  57, 189, 0.94585527, 0.37815695,
  46,   26, 0.05879151, 0.31761570,
  44, 142, 0.69983646, 0.27835633,
  31, 146, 0.62735035, 0.01168396,
  57,  53, 0.07166746, 0.02867074,
  60, 200, 1.00000000, 0.01887049
), ncol = 4, byrow = TRUE)

standard_sim <- function(x)fit_model(data_sim(x))

if(file.exists(paste0(getwd(),"/results/sim_example.RDS"))){
  cl <- parallel::makeCluster(parallel::detectCores()-2)
  parallel::clusterEvalQ(cl,library(glmmrBase))
  parallel::clusterEvalQ(cl,library(lme4))
  parallel::clusterExport(cl, list("data_sim", "fit_model", "standard_sim", "vals"))

  # Function to run simulations for a given point
  run_simulation_for_point <- function(index) {
    point <- vals[index, ]
    pbapply::pbreplicate(1000, standard_sim(point), cl = cl)
  }

  # Run simulations for the provided points
  results <- list()
  for (i in 1:nrow(vals)) {
    results[[i]] <- run_simulation_for_point(i)
  }

  saveRDS(results, paste0(getwd(),"/results/sim_example.RDS")) # change the path as required
  parallel::stopCluster(cl)
} else {
  results <- readRDS(paste0(getwd(),"/results/sim_example.RDS"))
}





# Compute true values for the provided points
true_values <- sapply(results, function(res) mean(res["p1", ]))
t1e_values <- sapply(results, function(res) quantile(res["p1", ],0.025) < 0.05 & quantile(res["p1", ],0.975) > 0.05)

mod_additive <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001), par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 10000)
mod_addint <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001), par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 1000)
mod_linear <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001), par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 10000)
mod_as <- adapt$new(data_fn = "data_sim", fit_fn = "fit_model", par_lower = c(6,5,0,0.001), par_upper = c(60,200,1,0.4), par_discrete = c(TRUE,FALSE,FALSE,FALSE), n = 10000)


cl <- parallel::makeCluster(4) # you can change the number of cores
parallel::clusterEvalQ(cl,library(glmmrBase))
parallel::clusterEvalQ(cl,library(lme4))
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

data_list <- readRDS("D:/data_list_2.RDS")
par_list <- readRDS("D:/par_list_2.RDS")

for (mod_type_name in c("mod_as")) {  # Only include the models you want to run
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
                         mod_linear = "linear",
                         mod_as = "as")

    mod_type$adapt(stat="p1", # model for the type 1 error
                   model="binomial",
                   type = model_type, # Use appropriate model type
                   sampling = 150,
                   L=1.1,
                   d = 2,
                   m = 10)

    pred_points <- vals
    pred_points_data <- mod_type$predict(pred_points)

    rmse <- c()
    for(i in 1:nrow(pred_points_data)) rmse <- c(rmse, sqrt(mean((pred_points_data[i, ] - true_values[i])^2)))

    pred_t1e <- c()
    for(i in 1:nrow(pred_points_data)){
     pred_t1e <- c(pred_t1e, quantile(pred_points_data[i, ],0.025) < 0.05 & quantile(pred_points_data[i, ],0.975) > 0.05)
    }
    #sum(pred_t1e*t1e_values)/length(t1e_values)

    cl_plot_points <- as.matrix(expand.grid(cl = 6:60, m = 10, cv = 0.2, icc = 0.05))
    cv_plot_points <- as.matrix(expand.grid(cl = 30, m = 10, cv = seq(0, 1, length.out = 100), icc = 0.05))
    icc_plot_points <- as.matrix(expand.grid(cl = 30, m = 10, cv = 0.2, icc = seq(0.001, 0.4, length.out = 100)))

    cl_plot_points_data <- mod_type$predict(cl_plot_points)
    cv_plot_points_data <- mod_type$predict(cv_plot_points)
    icc_plot_points_data <- mod_type$predict(icc_plot_points)

    dfp_cl <- data.frame(cl = 6:60,
                         mean = apply(cl_plot_points_data,1,mean),
                         lci = apply(cl_plot_points_data,1,function(i)quantile(i,0.025)),
                         uci = apply(cl_plot_points_data,1, function(i)quantile(i,0.975)))

    dfp_cv <- data.frame(cv = seq(0, 1, length.out = 100),
                         mean = apply(cv_plot_points_data,1,mean),
                         lci = apply(cv_plot_points_data,1,function(i)quantile(i,0.025)),
                         uci = apply(cv_plot_points_data,1, function(i)quantile(i,0.975)))

    dfp_icc <- data.frame(icc = seq(0.001, 0.4, length.out = 100),
                          mean = apply(icc_plot_points_data,1,mean),
                          lci = apply(icc_plot_points_data,1,function(i)quantile(i,0.025)),
                          uci = apply(icc_plot_points_data,1, function(i)quantile(i,0.975)))

    results_list[[mod_type_name]][[which(blocks == block)]] <- list(pred_points_data = pred_points_data,
                                                                    dfp_cl = dfp_cl,
                                                                    dfp_cv = dfp_cv,
                                                                    dfp_icc = dfp_icc,
                                                                    rmse = rmse,
                                                                    pred_t1e = pred_t1e)
    if(gen_new_points){
      tryCatch({
        mod_type$sample(block_size, type = "none")
      }, error = function(e) {
        message("Error during model sampling: ", e)
        next
      })
    }

    saveRDS(results_list, paste0(getwd(),"/results/results_list.RDS"))
    saveRDS(data_list, paste0(getwd(),"/results/data_list.RDS"))
    saveRDS(par_list, paste0(getwd(),"/results/par_list_3.RDS"))
  }
}

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
