require(glmmrBase)
require(pbapply)
if(!require(AdaptSim)){
  devtools::install_github("samuel-watson/AdaptSim") # note you will need the package devtools for this
  require(AdaptSim)
}


cond_number <- function(t,k,m,icc,cac,type){
  df <- data_trial(t,k,m,type)
  nT <- length(unique(df$t))
  model <- Model$new(
    formula = ~ int + factor(t) - 1 + (1|gr(cl)) + (1|gr(cl,t)),
    covariance = c(icc*cac,icc*(1-cac)),
    mean = rep(0,nT+1),
    data = df,
    family = gaussian(),
    var_par = 1-icc
  )
  df$y <- model$sim_data()
  M <- model$information_matrix()
  K <<- model$small_sample_correction("KR")
  e <- eigen(M)$values
  CAB <- M[1,2:nT,drop=FALSE] %*% solve(M[2:nT, 2:nT]) %*% t(M[1,2:nT,drop=FALSE])
  asd <- drop(solve(M)[1,1])
  
  return(list(df, max(e)/min(e), max(e), min(e), prod(e), 
              sum(e), length(e), drop(CAB), drop(M[1,1]), asd))
}

fit_model <- function(data=data){
  fit1 <- lme4::lmer(y~int+factor(t) + (1|cl/t),data=data[[1]], REML= FALSE)#tryCatch(suppressWarnings(lme4::lmer(y~int+(1|cl),data=data)),error=function(i)NA)
  
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
                         (s1$coefficients[2,1] - qnorm(0.975)*s1$coefficients[2,2]) < 0)*1,
             e1 = log(data[[2]]),
             e2 = log(data[[3]]),
             e3 = log(data[[4]]),
             e4 = log(data[[5]]),
             e5 = log(data[[6]]),
             elen = data[[7]],
             cab = data[[8]],
             m1 = data[[9]],
             a = data[[10]]))
  } else {
    # return NA if the model fails to fit
    return(c(conv = 1,
             beta = NA,
             se = NA,
             p = NA,
             p1 = NA,
             cover = NA,
             e1 = log(data[[2]]),
             e2 = log(data[[3]]),
             e3 = log(data[[4]]),
             e4 = log(data[[5]]),
             e5 = log(data[[6]]),
             elen = data[[7]],
             cab = data[[8]],
             m1 = data[[9]],
             a = data[[10]]))
  }
}

tmp <- fit_model(cond_number(10,2,20,0.001,0.1,"sw"))

# generate comparisons
block_size <- 20 

vals_e <- data.frame(
  t = sample(4:10,block_size,replace=T),
  k =  sample(1:2,block_size,replace = T),
  m=  sample(10:100,block_size,replace = T),
  icc =  runif(block_size,0.001,0.2),
  cac =  runif(block_size,0.01,0.99),
  type = sample(c("sw","inc","st","par","par_co"),block_size,replace=T)
)
vals_e[vals_e$type%in%c("par","par_co"),"k"] <- sample(4:40,sum(vals_e$type%in%c("par","par_co")),replace=T)

cl <- parallel::makeCluster(4)
parallel::clusterEvalQ(cl,library(glmmrBase))
parallel::clusterEvalQ(cl,library(lme4))
parallel::clusterExport(cl,c("cond_number","fit_model","data_trial","vals_e"))

# Function to run simulations for a given point

# Run simulations for the provided points
results_e <- list()
for (i in 18:nrow(vals_e)) {
  parallel::clusterExport(cl,c("i"))
  results_e[[i]] <- pbapply::pbreplicate(1000, fit_model(cond_number(vals_e$t[i],vals_e$k[i],
                                                                   vals_e$m[i],vals_e$icc[i],
                                                                   vals_e$cac[i],vals_e$type[i])), cl = cl)
}




saveRDS(vals_e,    paste0(getwd(),"/results/vals_eigen.RDS"))
saveRDS(results_e, paste0(getwd(),"/results/sim_example_eigen.RDS"))

vals_e <- readRDS(paste0(getwd(),"/results/vals_eigen.RDS"))
results_e <- readRDS(paste0(getwd(),"/results/sim_example_eigen.RDS"))

true_values <- sapply(results_e, function(res) mean(res["p1", ]))
t1e_values <- sapply(true_values, function(x) x - qnorm(0.975)*sqrt(x*(1-x)/1000) < 0.05 & x + qnorm(0.975)*sqrt(x*(1-x)/1000) > 0.05)

results_edf1 <- data.frame(t(results_e[[1]][7:14,1]))
for(i in 2:length(results_e)){
  results_edf1 <- rbind(results_edf1, data.frame(t(results_e[[i]][7:14,1])))
}

mod_additive <- adapt$new(data_fn = "data_trial", 
                          fit_fn = "fit_model", 
                          par_lower = c(1,2), 
                          par_upper = c(6,9), 
                          par_discrete = c(FALSE,FALSE), 
                          n = 20000)

cl <- parallel::makeCluster(4) # you can change the number of cores
parallel::clusterEvalQ(cl,.libPaths("C:/R/lib"))
parallel::clusterEvalQ(cl,library(glmmrBase))
parallel::clusterEvalQ(cl,library(lme4))
parallel::clusterExport(cl,c("cond_number","fit_model","data_trial"))

#########################################
# code block
#########################################
block_size <- 20000

vals <- data.frame(
  t = sample(4:10,block_size,replace=T),
  k =  sample(1:2,block_size,replace = T),
  m=  sample(10:30,block_size,replace = T),
  icc =  runif(block_size,0.001,0.2),
  cac =  runif(block_size,0.01,0.99),
  type = sample(c("sw","inc","st","par","par_co"),block_size,replace=T)
)
vals[vals$type%in%c("par","par_co"),"k"] <- sample(4:30,sum(vals$type%in%c("par","par_co")),replace=T)

parallel::clusterExport(cl,c("vals"))

res <- pbapply::pbsapply(1:nrow(vals),
                  function(i)fit_model(cond_number(vals$t[i],vals$k[i],
                                                   vals$m[i],vals$icc[i],
                                                   vals$cac[i],vals$type[i])),
                  cl = cl)

res2 <- res
res3 <- res
res4 <- res

saveRDS(list(res2,res3,res4), paste0(getwd(),"/results/sim_outputs.RDS"))

pars4 <- t(res3[c(7:14),])
pars4[,8] <- log(pars4[,8])

mod_additive$last_sim_output <- res3[1:6,]
mod_additive$par_vals <- pars4[,c(1,8)]

parallel::stopCluster(cl)


mod_additive$adapt(stat="p1", # model for the type 1 error
               model="binomial",
               type = "full", # Use appropriate model type
               sampling = 150,
               L=1.1,
               d = 1,
               m = 10)

results_edf1$m1 <- log(results_edf1$m1)
pred_e <- mod_additive$predict(results_edf1[,c(1,8)])
dfpr <- cbind(data.frame(pred = rowMeans(pred_e), true = true_values, 
                         lci = apply(pred_e,1,function(i)quantile(i,0.025)),
                         uci = apply(pred_e,1,function(i)quantile(i,0.975)),
                         prob = apply(pred_e,1,function(i)length(i[i>0.05])/length(i))),
              vals_e[1:17,], results_edf1)

rowMeans(pred_e)

true_values

dfpr$type_label <- factor(dfpr$type,
                          levels = unique(dfpr$type),
                          labels = c(
                            "Stepped-wedge",
                            "Parallel",
                            "Parallel XO",
                            "Incomplete S-W",
                            "Staircase"
                          ))

ggplot(data = dfpr, aes(x = log(true), y = log(pred), color = type_label))+
  geom_abline(intercept = 0, slope = 1, lty = 2)+
  geom_hline(yintercept = log(0.05),lty =3)+
  geom_vline(xintercept = log(0.05),lty =3)+
  geom_point()+
  geom_errorbar(aes(ymin = log(lci), ymax= log(uci)))+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_color_discrete(name = "Trial type")

ggplot(data = dfpr, aes(x = true, y = 2*(1-prob), color = type_label))+
  geom_vline(xintercept = 0.05+3*sqrt(0.05*0.95/10000),lty =3)+
  geom_vline(xintercept = 0.05-3*sqrt(0.05*0.95/10000),lty =3)+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x= "Log type I error", y = "'p-value'")+
  scale_color_discrete(name = "Trial type")+
  scale_x_log10()

dfv <- cbind(vals_e[1:17,], results_edf1, true_values)
dfv$isin <- I(true_values < 0.07)

ggplot(data = dfv, aes( y = e1, x = m1, 
                        color = true_values,
                        shape = factor(isin)))+
  geom_point()+
  scale_colour_viridis_c()

ggplot(data = dfv, aes( y = true_values, x = m1, color = type))+
  geom_point()


dfpv <- cbind(vals, mod_additive$par_vals)
qplot(x = true_values, y = rowMeans(pred_e))+
  geom_abline(intercept = 0, slope = 1)

epoints1 <- expand.grid(l = seq(1,6,length.out = 20), s = 4)
epoints2 <- expand.grid(l = 4, s = seq(3,9,length.out = 20))
epoints3 <- expand.grid(l = seq(1,6,length.out = 20), s = seq(3,9,length.out = 20))

plot_points_data1 <- mod_additive$predict(epoints1)
plot_points_data2 <- mod_additive$predict(epoints2)
plot_points_data3 <- mod_additive$predict(epoints3)

dfp_cl <- data.frame(l = epoints1$l,
                     s = epoints1$s,
                     mean = apply(plot_points_data1,1,mean),
                     lci = apply(plot_points_data1,1,function(i)quantile(i,0.025)),
                     uci = apply(plot_points_data1,1, function(i)quantile(i,0.975)))


p_cl <- ggplot(data = dfp_cl, aes(x = l, y = mean)) +
  geom_hline(yintercept = 0.05, lty = 2) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "cn", y = "Type I error", fill = "Iterations") +
  theme_minimal(); p_cl

dfp_cl <- data.frame(l = epoints2$l,
                     s = epoints2$s,
                     mean = apply(plot_points_data2,1,mean),
                     lci = apply(plot_points_data2,1,function(i)quantile(i,0.025)),
                     uci = apply(plot_points_data2,1, function(i)quantile(i,0.975)))


p_cl <- ggplot(data = dfp_cl, aes(x = s, y = mean)) +
  geom_hline(yintercept = 0.05, lty = 2) +
  geom_ribbon(aes(ymin = lci, ymax = uci), alpha = 0.2) +
  geom_line() +
  labs(x = "l", y = "Type I error", fill = "Iterations") +
  theme_minimal(); p_cl

####################################

dfp_cl <- data.frame(l = epoints3$l,
                     s = epoints3$s,
                     mean = apply(plot_points_data3,1,mean),
                     lci = apply(plot_points_data3,1,function(i)quantile(i,0.025)),
                     uci = apply(plot_points_data3,1, function(i)quantile(i,0.975)))


p_cl <- ggplot(data = dfp_cl, aes(x = s, y = l, fill = mean)) +
  geom_tile()+
  scale_fill_viridis_c(name = "Type I Error")+
  theme_minimal()+
  labs(x = "Log information", y = "Log condition number"); p_cl



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
                   model="binomial",
                   type = model_type, # Use appropriate model type
                   sampling = 150,
                   L=1.1,
                   d = ifelse(mod_type_name == "mod_linear",1,2),
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
    saveRDS(par_list, paste0(getwd(),"/results/par_list.RDS"))
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
