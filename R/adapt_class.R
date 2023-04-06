
## NOTE: I'VE CREATED A NEW CLASS USING R6 THAT HOLDS ALL OF THE FUNCTIONS
## AND OUTPUT SO THAT WE CAN KEEP TRACK OF PREVIOUS VALUES AND THE
## NUMBER OF SIMULATIONS. THIS SHOULD ALSO MAKE SAVING AND RUNNING EASIER.
## YOU CAN SEE EXAMPLE OF HOW TO RUN IT IN THE OTHER FILE.

adapt <- R6::R6Class("adapt",
                     public = list(
                       #name of data simulation function
                       data_fn = NULL,
                       #name of model fitting function
                       fit_fn = NULL,
                       # lower bounds
                       par_lower = NULL,
                       # upper bounds
                       par_upper = NULL,
                       # whether the parameter is discrete
                       par_discrete = NULL,
                       # basis functions per dimension
                       m = NULL,
                       # simulation values of parameters
                       par_vals = NULL,
                       # adaption data from last adapt run,
                       adapt_output = NULL,
                       # save last stan fit in case
                       last_stan_fit = NULL,
                       # the output of the lat sim function
                       last_sim_output = NULL,
                       # last kappa value
                       kappa = NULL,
                       # last lambda values
                       lambda = NULL,
                       # initialisation
                       # args are self-explanatory
                       # n is initial starting number of simulations
                       initialize = function(data_fn,
                                             fit_fn,
                                             par_lower,
                                             par_upper,
                                             n,
                                             m,
                                             par_discrete = rep(FALSE,length(par_upper))){
                         ## add checks of the functions here

                         self$par_lower = par_lower
                         self$par_upper = par_upper
                         self$par_discrete = par_discrete
                         self$data_fn = get(data_fn)
                         self$fit_fn = get(fit_fn)
                         self$m = m
                         args <- list(1:m)
                         if(length(par_upper)>1){
                           for(i in 2:length(par_upper)){
                             args[[i]] <- 1:m
                           }
                         }
                         private$ind <- as.matrix(do.call(expand.grid,args))
                         private$priors_m = rep(0,nrow(private$ind))
                         private$priors_sd = rep(1,nrow(private$ind))
                         private$prior_intercept = c(0,1)
                         private$prior_lengthscale = as.matrix(cbind(rep(0,length(par_upper)),rep(0.5,length(par_upper))))
                         private$prior_fscale = c(0,0.5)
                         private$prior_varpar = c(0,1)
                         self$par_vals <- matrix(0,nrow=n,ncol=length(par_upper))
                         # simulate starting values
                         for(i in 1:length(par_upper)){
                           self$par_vals[,i] <- runif(n,par_lower[i],par_upper[i])
                           if(self$par_discrete[i])self$par_vals[,i] <- round(self$par_vals[,i],0)
                         }
                         private$all_vals = self$par_vals
                         message("Stan files will be compiled if this is the first time executing this class.")
                         bin_file <- system.file("stan","approx_gp_binom.stan",package = "AdaptSim",mustWork = TRUE)
                         lin_file <- system.file("stan","approx_gp_lin.stan",package = "AdaptSim",mustWork = TRUE)

                         private$mod_bin = cmdstanr::cmdstan_model(bin_file)
                         private$mod_lin = cmdstanr::cmdstan_model(lin_file)
                       },
                       # model fitting simulation
                       sim = function(cl = NULL){
                         if(is.null(cl)){
                           out <- pbapply::pbsapply(1:nrow(self$par_vals),function(i){
                             self$fit_fn(self$data_fn(self$par_vals[i,]))
                           })
                         } else {
                           out <- pbapply::pbsapply(1:nrow(self$par_vals),function(i){
                             self$fit_fn(self$data_fn(self$par_vals[i,]))
                           },cl=cl)
                         }
                         self$last_sim_output = out
                       },
                       adapt = function(stat,
                                        model= "binomial",
                                        alpha = 0.05,
                                        kappa = 1,
                                        samp_n = 20,
                                        n_new = 100,
                                        L = 1.2,
                                        append = FALSE,
                                        chains = 3,
                                        parallel_chains = 6,
                                        warmup = 500,
                                        sampling = 500){

                         if(any(is.na(self$last_sim_output[2,]) | is.na(self$last_sim_output[stat,]))){
                           idxna <- which((is.na(self$last_sim_output[2,]) | is.na(self$last_sim_output[stat,])))
                           self$par_vals <- self$par_vals[-idxna,]
                           self$last_sim_output <- self$last_sim_output[,-idxna]
                         }

                         if(is(samp_n,"data.frame")){
                           if(ncol(samp_n)!=length(self$par_upper))stop("Wrong number of columns")
                           for(i in 1:length(self$par_upper)){
                             if(any(samp_n[,i]>self$par_upper[i])|any(samp_n[,i]<self$par_lower[i]))warning("Prediction values out of range")
                           }
                           xs <- samp_n
                           colnames(xs) <- paste0("Var",1:ncol(xs))
                           x <- rbind(self$par_vals,as.matrix(xs))
                         } else if(is(samp_n,"double")){
                           q1 <- self$par_upper[1] - self$par_lower[1]
                           l1 <- self$par_lower[1] + q1/(2*samp_n)
                           u1 <- self$par_upper[1] - q1/(2*samp_n)
                           xs <- matrix(seq(l1,u1,length.out=samp_n),ncol=1)
                           if(length(self$par_upper)>1){
                             for(i in 2:length(self$par_upper)){
                               q1 <- self$par_upper[i] - self$par_lower[i]
                               l1 <- self$par_lower[i] + q1/(2*samp_n)
                               u1 <- self$par_upper[i] - q1/(2*samp_n)
                               xs <- expand.grid(xs, matrix(seq(l1,u1,length.out=samp_n),ncol=1))
                             }
                           }
                           x <- rbind(self$par_vals,as.matrix(xs))

                         } else {
                           stop("samp_n wrong type")
                         }
                         nsamp <- nrow(xs)

                         x_grid <- x
                         for(i in 1:ncol(x_grid)){
                           x_grid[,i] <- 2*(x_grid[,i]-min(x_grid[,i]))/(max(x_grid[,i])-min(x_grid[,i])) - 1
                         }

                         dat <- list(
                           D = ncol(self$par_vals),
                           L = rep(L,length(self$par_upper)),
                           M = self$m,
                           M_nD = nrow(private$ind),
                           Nsample = nrow(self$par_vals),
                           Npred = nsamp,
                           y = self$last_sim_output[stat,],
                           x_grid = x_grid,
                           indices = private$ind,
                           intercept_prior = private$prior_intercept,
                           beta_prior = cbind(private$priors_m,private$priors_sd),
                           lengthscale_prior = private$prior_lengthscale,
                           fscale_prior = private$prior_fscale
                         )

                         if(private$prior_intercept[1] == 0) {
                          if(model == "binomial"){
                            dat$intercept_prior[1] <- log(alpha/(1-alpha))
                          } else if(model == "linear"){
                            dat$intercept_prior[1] <- alpha
                          }
                         }

                         if(model == "binomial"){
                           message("Using binomial model")
                           fit <- private$mod_bin$sample(data = dat,
                                                         chains = chains,
                                                         parallel_chains = parallel_chains,
                                                         iter_warmup = warmup,
                                                         iter_sampling = sampling,
                                                         refresh = 100)
                         } else if(model == "linear"){
                           message("Using linear model")
                           dat$
                           fit <- private$mod_lin$sample(data = dat,
                                                         chains = chains,
                                                         parallel_chains = parallel_chains,
                                                         iter_warmup = warmup,
                                                         iter_sampling = sampling,
                                                         refresh = 100)
                         } else {
                           stop("Model incorrectly specified")
                         }

                         self$last_stan_fit = fit
                         ypred <- fit$draws("y_grid_predict",format = "matrix")
                         sige <- fit$draws("sigma_e",format = "matrix")
                         dfp <- cbind(as.data.frame(xs),data.frame(mean = colMeans(ypred)))
                         dfp$lci <- apply(ypred,2,function(i)quantile(i,0.025))
                         dfp$uci <- apply(ypred,2,function(i)quantile(i,0.975))

                         dfp$prob <- NA
                         dfp$entr <- NA

                         ## average correlation kappa estimates commented out below.
                         ## see code below that for updated mean CrI width approach

                         lambda <- fit$draws("phi",format = "matrix")
                         lambda <- colMeans(lambda)
                         self$lambda <- lambda

                         # if(kappa == 1){
                         #   k <- private$kappa_1(alpha,lambda)
                         #   self$kappa <- k
                         #   for(i in 1:nrow(dfp)){
                         #     dfp$prob[i] <- sum(abs(ypred[,i]-alpha) < k)/nrow(ypred)
                         #     dfp$entr[i] <- -dfp$prob[i]*log(dfp$prob[i],2) - (1-dfp$prob[i])*log((1-dfp$prob[i]),2)
                         #     if(dfp$prob[i]==0 || dfp$prob[i]==1)dfp$entr[i] <- 0
                         #   }
                         # } else if(kappa== 2){
                         #   k <- rep(NA,nrow(dfp))
                         #   for(i in 1:nrow(dfp)){
                         #     k[i] <- private$kappa_2(theta = unlist(unname(dfp[i,1:ncol(self$par_vals)])),alpha,lambda)
                         #
                         #     dfp$prob[i] <- sum(abs(ypred[,i]-alpha) < k[i])/nrow(ypred)
                         #     dfp$entr[i] <- -dfp$prob[i]*log(dfp$prob[i],2) - (1-dfp$prob[i])*log((1-dfp$prob[i]),2)
                         #     if(dfp$prob[i]==0 || dfp$prob[i]==1)dfp$entr[i] <- 0
                         #   }
                         #   self$kappa <- cbind(dfp[,1:ncol(self$par_vals)],k=k)
                         # }

                         if(kappa==1){
                           self$kappa <- mean(abs(dfp$uci - dfp$lci)/2)
                           for(i in 1:nrow(dfp)){
                                 dfp$prob[i] <- sum(abs(ypred[,i]-alpha) < self$kappa)/nrow(ypred)
                                 dfp$entr[i] <- -dfp$prob[i]*log(dfp$prob[i],2) - (1-dfp$prob[i])*log((1-dfp$prob[i]),2)
                                 if(dfp$prob[i]==0 || dfp$prob[i]==1)dfp$entr[i] <- 0
                               }
                         } else if(kappa == 2) {
                           self$kappa <- abs(dfp$uci - dfp$lci)/2
                           for(i in 1:nrow(dfp)){
                             dfp$prob[i] <- sum(abs(ypred[,i]-alpha) < self$kappa[i])/nrow(ypred)
                             dfp$entr[i] <- -dfp$prob[i]*log(dfp$prob[i],2) - (1-dfp$prob[i])*log((1-dfp$prob[i]),2)
                             if(dfp$prob[i]==0 || dfp$prob[i]==1)dfp$entr[i] <- 0
                           }
                         }


                         # get posterior values
                         beta <- fit$draws("beta",format = "matrix")
                         private$priors_m <- colMeans(beta)
                         private$priors_sd <- apply(beta,2,sd)
                         intercept <- fit$draws("intercept",format = "matrix")
                         private$prior_intercept <- c(mean(intercept),sd(intercept))
                         phi <- fit$draws("phi",format = "matrix")
                         for(i in 1:length(self$par_upper)){
                           private$prior_lengthscale[,1] <- colMeans(phi)
                           private$prior_lengthscale[,2] <- apply(phi,2,sd)
                         }
                         sigmae <- fit$draws("sigma_e",format = "matrix")
                         private$prior_fscale <- c(mean(sigmae),sd(sigmae))
                         if(model%in%c("linear","beta")){
                           sigma <- fit$draws("sigma",format = "matrix")
                           private$prior_varpar <- c(mean(sigma),sd(sigma))
                         }

                         ## draw a new sample of parameters
                         newsamp <- sample(1:nrow(dfp),n_new,replace = TRUE,prob = dfp$entr)

                         d1 <- diff(dfp$Var1)[diff(dfp$Var1)>0][1]
                         dfs <- matrix(dfp$Var1[newsamp]+runif(length(newsamp),0,d1))

                         if(length(self$par_upper)>1){
                           for(i in 2:length(self$par_upper)){
                             d1 <- diff(dfp[,paste0("Var",i)])[diff(dfp[,paste0("Var",i)])>0][1]
                             dfs <- cbind(dfs,matrix(dfp[,paste0("Var",i)][newsamp]+runif(length(newsamp),0,d1)))
                           }
                         }

                         for(i in 1:length(self$par_upper)){
                           if(self$par_discrete[i])dfs[,i] <- round(dfs[,i],0)
                         }

                         if(append){
                           self$par_vals <- rbind(self$par_vals,as.matrix(dfs))
                         } else {
                           self$par_vals <- as.matrix(dfs)
                         }

                         private$all_vals <- rbind(private$all_vals,as.matrix(dfs))
                         self$adapt_output = dfp
                       },
                       n_sims = function(){
                         nrow(private$all_vals)
                       },
                       priors = function(){
                         dfp <- as.data.frame(private$ind)
                         dfp$mean <- private$priors_m
                         dfp$sd <- private$priors_sd
                         return(list(beta = dfp, intercept = private$prior_intercept))
                       }
                     ),
                     private = list(
                       priors_m = NULL,
                       priors_sd = NULL,
                       prior_intercept = NULL,
                       prior_lengthscale = NULL,
                       prior_fscale = NULL,
                       prior_varpar = NULL,
                       ind = NULL,
                       all_vals = NULL,
                       mod_lin = NULL,
                       mod_bin = NULL,
                       kappa_1 = function(alpha, lambda, sigma_e = NULL){
                         N <- nrow(private$all_vals)
                         rho <- rep(NA,length(self$par_upper))
                         for(i in 1:length(self$par_upper)){
                           rho[i] <- lambda[i] * sqrt(pi) * pracma::erf(self$par_upper[i]/lambda[i]) -
                             lambda[i]^2/6 * (1- exp(-lambda[i]^2/self$par_upper[i]^2))
                         }
                         n <- 1 + (N - 1)*prod(rho)
                         if(is.null(sigma_e)){
                           return(1.96 * sqrt(alpha*(1-alpha)/n))
                         } else {
                           ## add return for linear model...
                         }
                       },
                       kappa_2 = function(theta, alpha, lambda, sigma_e = NULL){
                         N <- nrow(private$all_vals)
                         b <- self$par_upper - self$par_lower
                         rho <- rep(NA,length(b))
                         for(i in 1:length(b)){
                           l <- 0.5/b[i] * sqrt(pi) * lambda[i]
                           if(theta[i] == self$par_lower[i] || theta[i] == self$par_upper[i]){
                             rho[i] <- l * pracma::erf(b[i]/lambda[i])
                           } else if(theta[i] < (self$par_lower[i]+self$par_upper[i])/2){
                             rho[i] <- l * (pracma::erf((theta[i]-self$par_lower[i])/lambda[i]) + pracma::erf((self$par_upper[i] - theta[i])/lambda[i]))
                           } else if(theta[i] == (self$par_lower[i]+self$par_upper[i])/2){
                             rho[i] <- 2*l * pracma::erf(b[i]/2*lambda[i])
                           } else if(theta[i] > (self$par_lower[i]+self$par_upper[i])/2){
                             rho[i] <- l * (pracma::erf((self$par_upper[i] - theta[i])/lambda[i]) + pracma::erf((theta[i]-self$par_lower[i])/lambda[i]))
                           }
                         }
                         n <- 1 + (N - 1)*prod(rho)
                         return(1.96 * sqrt(alpha*(1-alpha)/n))
                       }
                     ))
