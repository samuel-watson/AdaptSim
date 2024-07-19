
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
                                             par_discrete = rep(FALSE,length(par_upper))){
                         ## add checks of the functions here

                         self$par_lower = par_lower
                         self$par_upper = par_upper
                         self$par_discrete = par_discrete
                         self$data_fn = get(data_fn)
                         self$fit_fn = get(fit_fn)
                         private$prior_intercept = c(0,1)
                         private$prior_fscale = c(0,0.5)
                         private$prior_varpar = c(1,1)
                         self$par_vals <- matrix(0,nrow=n,ncol=length(par_upper))
                         # simulate starting values
                         for(i in 1:n)self$par_vals[i,] <- private$sample_values()
                         private$all_vals = self$par_vals
                         message("Stan files will be compiled if this is the first time executing this class.")
                         bin_file <- system.file("stan","approx_gp_binom.stan",package = "AdaptSim",mustWork = TRUE)
                         lin_file <- system.file("stan","approx_gp_lin.stan",package = "AdaptSim",mustWork = TRUE)
                         bin_embed_file <- system.file("stan","approx_gp_binom_linembed.stan",package = "AdaptSim",mustWork = TRUE)
                         bin_additive_file <- system.file("stan","approx_gp_binom_additive.stan",package = "AdaptSim",mustWork = TRUE)
                         bin_as_file <- system.file("stan","approx_gp_binom_as.stan",package = "AdaptSim",mustWork = TRUE)
                         lin_embed_file <- system.file("stan","approx_gp_lin_linembed.stan",package = "AdaptSim",mustWork = TRUE)
                         lin_additive_file <- system.file("stan","approx_gp_lin_additive.stan",package = "AdaptSim",mustWork = TRUE)
                         lin_as_file <- system.file("stan","approx_gp_lin_as.stan",package = "AdaptSim",mustWork = TRUE)
                         private$mod_bin = cmdstanr::cmdstan_model(bin_file)
                         private$mod_bin_embed = cmdstanr::cmdstan_model(bin_embed_file)
                         private$mod_bin_additive = cmdstanr::cmdstan_model(bin_additive_file)
                         private$mod_bin_as = cmdstanr::cmdstan_model(bin_as_file)
                         private$mod_lin = cmdstanr::cmdstan_model(lin_file)
                         private$mod_lin_embed = cmdstanr::cmdstan_model(lin_embed_file)
                         private$mod_lin_additive = cmdstanr::cmdstan_model(lin_additive_file)
                         private$mod_lin_as = cmdstanr::cmdstan_model(lin_as_file)
                       },
                       # model fitting simulation
                       sim = function(cl = NULL){
                         if(is.null(cl)){
                           out <- pbapply::pbsapply(1:nrow(self$par_vals),function(i){
                             self$fit_fn(self$data_fn(self$par_vals[i,]))
                           })
                         } else {
                           message(paste0("Parallel execution with ",length(cl)," cores"))
                           # I don't know why this is slower than executing in the global environment
                           # can be improved
                           out <- pbapply::pbsapply(1:nrow(self$par_vals),function(i){
                             self$fit_fn(self$data_fn(self$par_vals[i,]))
                           },cl=cl)
                         }
                         self$last_sim_output = out
                       },
                       # type can be "linear" for linear embedding, "additive" for additive, "as" for active subspace, and "full" for complete D-dimensional model
                       adapt = function(stat,
                                        model= "binomial",
                                        type = "linear",
                                        alpha = 0.05,
                                        L = 1.2,
                                        m = 5,
                                        d = 2,
                                        append = FALSE,
                                        chains = 3,
                                        parallel_chains = 6,
                                        warmup = 500,
                                        sampling = 500){

                         if(!type%in%c("linear","additive","as","full"))stop("Type must be linear, additive, as, or full")
                         if(!model%in%c("binomial","linear"))stop("Model must be binomial or linear")
                         if(length(m)==1) {
                           self$m = rep(m,ifelse(type == "as", d, length(self$par_upper)))
                         } else {
                           if(((length(m)!=length(par_upper))&type!="as")|(((length(m)!=d&type=="as"))))stop("m wrong length, should be either length 1 or D/d")
                           self$m = m
                         }
                         if(length(self$par_upper)<=2 & type == "as")stop("No active subspace method for D < 3")
                         args <- list(1:self$m[1])
                         if(type %in% c("as","full")){
                           upper_dim <- ifelse(type == "as",d,length(self$par_upper))
                           for(i in 2:upper_dim){
                             args[[i]] <- 1:self$m[i]
                           }
                           private$ind <- as.matrix(do.call(expand.grid,args))
                         } else if(type == "linear") {
                           private$ind <- matrix(1:self$m[1],ncol=1)
                         } else {
                           total_fn <- sum(self$m)
                         }

                         if(is.null(private$priors_m))
                         {
                           private$priors_m = rep(0,ifelse(type == "additive", sum(self$m) ,nrow(private$ind)))
                           private$priors_sd = rep(1,ifelse(type == "additive", sum(self$m) ,nrow(private$ind)))
                           private$prior_lengthscale = as.matrix(cbind(rep(0,ifelse(type == "as", d, length(self$par_upper))),
                                                                       rep(0.5,ifelse(type == "as", d, length(self$par_upper)))))
                         }
                         if(type == "as"&is.null(private$a_mat_prior)){
                           private$a_mat_prior <- matrix(NA,2,d*length(self$par_upper))
                           private$a_mat_prior[1,] <- 0
                           private$a_mat_prior[2,] <- 1
                         }

                         if(any(is.na(self$last_sim_output[2,]) | is.na(self$last_sim_output[stat,])))
                         {
                           idxna <- which((is.na(self$last_sim_output[2,]) | is.na(self$last_sim_output[stat,])))
                           self$par_vals <- self$par_vals[-idxna,]
                           self$last_sim_output <- self$last_sim_output[,-idxna]
                         }

                         x_grid <- private$transform(as.matrix(self$par_vals))

                         if(length(L)==1)
                         {
                           if(type == "linear"){
                             private$lvec <- c(L*length(self$par_upper))
                           } else if(type == "as") {
                             private$lvec <- rep(L*d,d)
                           } else {
                             private$lvec <- rep(L,length(self$par_upper))
                           }
                         } else {
                           if(((length(L)!=length(self$par_upper))&type!="as")|((length(L)!=d)&type!="as"))stop("L wrong length")
                           private$lvec = L
                         }

                         if(type != "additive"){
                           LAMBDA <- private$lambda_nD(private$lvec,private$ind,ifelse(type == "linear",1,
                                                                                       ifelse(type == "as", d, length(self$par_upper))))
                         } else {
                           LAMBDA <- rep(NA, total_fn)
                           for(i in 1:ncol(x_grid)){
                             for(j in 1:self$m[i]){
                               mcol <- ifelse(i == 1, 0, sum(self$m[1:(i-1)]))
                               LAMBDA[mcol + j] <- private$lambda_nD(private$lvec[i],j,1)
                             }
                           }
                         }


                         if(type == "linear" & is.null(private$avec_prior)){
                             private$avec_prior <- rep(0,length(self$par_upper))
                             private$avec_conc <- 0
                         }
                         
                         if(type == "linear" & nrow(private$prior_lengthscale) != 1){
                           warning("Length scale prior wrong dimension for linear model, taking only the first entry")
                           private$prior_lengthscale <- private$prior_lengthscale[1,,drop=FALSE]
                         } 

                         dat <- list(
                           D = ncol(self$par_vals),
                           L = private$lvec,
                           Nsample = nrow(self$par_vals),
                           y = self$last_sim_output[stat,],
                           x_grid = x_grid,
                           lambda = LAMBDA,
                           intercept_prior = private$prior_intercept,
                           beta_prior = cbind(private$priors_m,private$priors_sd),
                           lengthscale_prior = private$prior_lengthscale,
                           fscale_prior = private$prior_fscale
                         )

                         if(type != "additive") {
                           dat <- append(dat,list(indices = private$ind,
                                                  M_nD = nrow(private$ind)))
                         } else {
                           dat <- append(dat,list(total_fn = sum(self$m),
                                                  M_nD = self$m))
                         }

                         if(type == "linear") dat <- append(dat,list(avec_prior = private$avec_prior,
                                                                     avec_conc = private$avec_conc))

                         if(type == "as"){
                           private$as_d <- d
                           dat <- append(dat,list(d = d,
                                                  a_mat_prior_mean = private$a_mat_prior[1,],
                                                  a_mat_prior_sd = private$a_mat_prior[2,]))
                         }

                         if(private$prior_intercept[1] == 0) {
                          if(model == "binomial"){
                            dat$intercept_prior[1] <- log(alpha/(1-alpha))
                          } else if(model == "linear"){
                            dat$intercept_prior[1] <- alpha
                          }
                         }

                         if(model=="linear"){
                           dat <- append(dat,list(sigma_prior = private$prior_varpar))
                         }

                         if(model == "binomial"){
                           message("Using binomial model")
                           if(type == "linear"){
                             fit <- private$mod_bin_embed$sample(data = dat,
                                                           chains = chains,
                                                           parallel_chains = parallel_chains,
                                                           iter_warmup = warmup,
                                                           iter_sampling = sampling,
                                                           refresh = 100)
                           } else if (type == "additive") {
                             fit <- private$mod_bin_additive$sample(data = dat,
                                                                 chains = chains,
                                                                 parallel_chains = parallel_chains,
                                                                 iter_warmup = warmup,
                                                                 iter_sampling = sampling,
                                                                 refresh = 100)
                           } else if (type == "as") {
                             fit <- private$mod_bin_as$sample(data = dat,
                                                                    chains = chains,
                                                                    parallel_chains = parallel_chains,
                                                                    iter_warmup = warmup,
                                                                    iter_sampling = sampling,
                                                                    refresh = 100)
                           } else {
                             fit <- private$mod_bin$sample(data = dat,
                                                           chains = chains,
                                                           parallel_chains = parallel_chains,
                                                           iter_warmup = warmup,
                                                           iter_sampling = sampling,
                                                           refresh = 100)
                           }
                         } else if(model == "linear"){
                           message("Using linear model")
                           dat$sigma_prior <- private$prior_varpar
                           fit <- private$mod_lin$sample(data = dat,
                                                         chains = chains,
                                                         parallel_chains = parallel_chains,
                                                         iter_warmup = warmup,
                                                         iter_sampling = sampling,
                                                         refresh = 100)
                         } else {
                           stop("Model incorrectly specified")
                         }

                         private$model <- model
                         private$type <- type
                         self$last_stan_fit = fit

                         # # get posterior values
                         beta <- fit$draws("beta",format = "matrix")
                         private$weights <- beta
                         private$priors_m <- colMeans(beta)
                         private$priors_sd <- apply(beta,2,sd)
                         private$intercept <- fit$draws("intercept",format = "matrix")
                         private$prior_intercept <- c(mean(private$intercept),sd(private$intercept))
                         private$phi <- fit$draws("phi",format = "matrix")
                         if(type == "linear")private$a_vector <- fit$draws("a",format = "matrix")
                         for(i in 1:length(self$par_upper)){
                           private$prior_lengthscale[,1] <- colMeans(private$phi)
                           private$prior_lengthscale[,2] <- apply(private$phi,2,sd)
                         }
                         private$sigmae <- fit$draws("sigma_e",format = "matrix")
                         private$prior_fscale <- c(mean(private$sigmae),sd(private$sigmae))
                         if(model%in%c("linear","beta")){
                           sigma <- fit$draws("sigma",format = "matrix")
                           m_sigma <- mean(sigma)
                           v_sigma <- var(sigma)
                           b <- m_sigma/v_sigma
                           a <- m_sigma*b
                           private$prior_varpar <- c(a,b)
                         }
                         if(type == "linear"){
                           private$avec_prior <- colMeans(private$a_vector)
                           private$avec_conc <- 1/var(drop(matrix(private$a_vector)))
                         }
                         if(type == "as"){
                           private$a_mat <- fit$draws("A",format = "matrix")
                           aa <- mod$last_stan_fit$draws("Amat",format="matrix")
                           private$a_mat_prior[1,] <- colMeans(private$a_mat)
                           private$a_mat_prior[2,] <- apply(private$a_mat,2,sd)
                         }
                       },
                       sample = function(n, type = "var", append = FALSE, kappa = NULL, alpha = 0.05){
                         if(is.null(private$intercept))stop("No MCMC samples")
                         x <- private$sample_values()
                         x <- matrix(x,nrow=1)
                         p2 <- self$predict(x)
                         if(type == "var"){
                            e1 <- var(p2[1,])
                         } else if(type == "entr") {
                           if(private$model == "binomial"){
                             if(is.null(kappa)){
                               prob1 <- mean(p2[1,] > alpha)
                             } else {
                               prob1 <- mean(p2[1,] - alpha < kappa)
                             }
                             e1 <- -prob1*log(prob1,2) - (1-prob1)*log((1-prob1),2)
                           } else {
                             stop("Entropy not compatible with linear model")
                           }
                         } else if (type == "none") {
                           x <- matrix(NA,nrow=n, ncol= length(self$par_upper))
                           for(i in 1:n)x[i,] <- private$sample_values()
                         } else {
                           stop("Type not recognized")
                         }

                         if(type != "none"){
                           for(i in 1:(n-1)){
                             x2 <- private$sample_values()
                             x2 <- matrix(x2,nrow=1)
                             p2 <- self$predict(x2)
                             if(type == "var"){
                               e2 <- var(p2[1,])
                             } else {
                               if(private$model == "binomial"){
                                 if(is.null(kappa)){
                                   prob1 <- mean(p2[1,] > alpha)
                                 } else {
                                   prob1 <- mean(p2[1,] - alpha < kappa)
                                 }
                                 e2 <- ifelse(prob1 == 1, 0, ifelse(prob1==0, 1, -prob1*log(prob1,2) - (1-prob1)*log((1-prob1),2)))
                               }
                             }
                             u1 <- runif(1)
                             if( e2/e1 > u1) {
                               x <- rbind(x,x2)
                               e1 <- e2
                             } else {
                               x <- rbind(x,x[nrow(x),,drop=FALSE])
                             }
                             cat("\rIter: ",i," of ",n)
                           }
                         }


                         private$all_vals <- rbind(private$all_vals,x)
                         if(append){
                           self$par_vals <- rbind(self$par_vals,x)
                         } else {
                           self$par_vals <- x
                         }
                       },
                       predict = function(x){
                         message(paste0("Predicting from model type: ",private$type))
                         if(is.null(private$intercept))stop("No MCMC samples")
                         x_grid <- private$transform(x)
                         iter <- nrow(private$weights)
                         dim <- ifelse(private$type == "additive", sum(self$m), nrow(private$ind))
                         PHI <- matrix(0,nrow(x),dim)
                         vals <- matrix(NA,nrow=nrow(x),ncol=iter)
                         n_d <- ifelse(private$type == "as",private$as_d,length(self$par_upper))

                         if(private$type != "additive"){
                           LAMBDA <- private$lambda_nD(private$lvec,private$ind,ifelse(private$type=="linear",1,n_d))
                         } else {
                           LAMBDA <- rep(NA, dim)
                           for(i in 1:length(self$par_upper)){
                             for(j in 1:self$m[i]){
                               mcol <- ifelse(i == 1, 0, sum(self$m[1:(i-1)]))
                               LAMBDA[mcol + j] <- private$lambda_nD(private$lvec[i],j,1)
                             }
                           }
                         }

                         if(private$type=="linear"){
                           xa <- x_grid %*% t(private$a_vector)
                           diagSPD <- private$spd_nD(private$sigmae,private$phi,LAMBDA,1)
                           diagSPD <- sqrt(diagSPD)
                           spd_beta <- t(diagSPD) * t(private$weights)
                           for(j in 1:iter){
                             for(i in 1:dim) PHI[,i] <- private$phi_nD(private$lvec,unlist(private$ind[i,]),xa[,j])
                             vals[,j] <- c(private$intercept[j,]) + PHI %*% c(spd_beta[,j])
                             cat("\rIter: ",j," of ",iter)
                           }
                         } else if(private$type == "additive"){
                           diagSPD <- matrix(NA, nrow = iter, ncol = dim)
                           for(i in 1:ncol(x_grid)){
                             mcol <- ifelse(i == 1, 0, sum(self$m[1:(i-1)]))
                             diagSPD[,(mcol+1):(mcol+self$m[i])] <- private$spd_nD(private$sigmae,private$phi[,i],LAMBDA[(mcol+1):(mcol+self$m[i])],1)
                             for(j in 1:self$m[i]){
                               PHI[,mcol + j] <- private$phi_nD(private$lvec[i],j,x_grid[,i,drop=FALSE])
                             }
                           }
                           diagSPD <- sqrt(diagSPD)
                           vals <- c(private$intercept) + PHI %*% (t(diagSPD) * t(private$weights))
                         } else if(private$type=="as"){
                           diagSPD <- private$spd_nD(private$sigmae,private$phi,LAMBDA,private$as_d)
                           diagSPD <- sqrt(diagSPD)
                           spd_beta <- t(diagSPD) * t(private$weights)
                           for(j in 1:iter){
                             A <- matrix(c(private$a_mat[j,]),nrow=length(self$par_upper))
                             xa <- x_grid %*% A
                             for(i in 1:dim) PHI[,i] <- private$phi_nD(private$lvec,unlist(private$ind[i,]),xa)
                             vals[,j] <- c(private$intercept[j,]) + PHI %*% c(spd_beta[,j])
                             cat("\rIter: ",j," of ",iter)
                           }
                         } else {
                           for(i in 1:dim) PHI[,i] <- private$phi_nD(private$lvec,unlist(private$ind[i,]),x_grid)
                           diagSPD <- private$spd_nD(private$sigmae,private$phi,LAMBDA,length(self$par_upper))
                           diagSPD <- sqrt(diagSPD)
                           vals <- c(private$intercept) + PHI %*% (t(diagSPD) * t(private$weights))
                         }

                         if(private$model == "binomial") vals <- exp(vals)/(1+exp(vals))
                         return(vals)
                       },
                       n_sims = function(){
                         nrow(private$all_vals)
                       },
                       priors = function(){
                         dfp <- as.data.frame(private$ind)
                         dfp$mean <- private$priors_m
                         dfp$sd <- private$priors_sd
                         out <- list(beta = dfp,
                                     intercept = private$prior_intercept,
                                     lengthscale = private$prior_lengthscale,
                                     fvar = private$prior_fscale,
                                     variance = private$prior_varpar)
                         if(!is.null(private$type)){
                           if(private$type == "as"){
                             out <- append(out, list(a_mat = private$a_mat_priors))
                           } else if(private$type == "linear"){
                             out <- append(out, list(avec_prior = private$avec_prior,
                                                     avec_conc = private$avec_conc))
                           }
                         }
                         return(out)
                       },
                       set_priors = function(beta_mean = NULL,
                                             beta_sd = NULL,
                                             lengthscale = NULL,
                                             fvar = NULL,
                                             variance = NULL){
                         if(!is.null(beta_mean)){
                           if(nrow(private$ind)!=length(beta_mean))stop("beta mean wrong size")
                           private$priors_m = beta_mean
                         }
                         if(!is.null(beta_sd)){
                           if(nrow(private$ind)!=length(beta_sd))stop("beta sd wrong size")
                           private$priors_sd = beta_sd
                         }
                         if(!is.null(lengthscale)){
                           if(length(self$par_upper)!=nrow(lengthscale) & nrow(lengthscale) > 1)stop("lengthscale wrong size")
                           if(ncol(lengthscale)!=2)stop("lengthscale should have 2 columns for mean and sd")
                           private$prior_lengthscale = lengthscale
                         }
                         if(!is.null(beta_sd)){
                           if(length(fvar)!=2)stop("beta sd wrong size")
                           private$prior_fscale = fvar
                         }
                         if(!is.null(variance)){
                           if(length(fvar)!=2)stop("beta sd wrong size")
                           private$prior_varpar = variance
                         }
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
                       mod_bin_embed = NULL,
                       mod_bin_additive = NULL,
                       mod_bin_as = NULL,
                       mod_lin_embed = NULL,
                       mod_lin_additive = NULL,
                       mod_lin_as = NULL,
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
                       },
                       weights = NULL,
                       lambda_nD = function(L, m, D){
                         lam <- t(m) * (pi/(2*L))
                         return(lam^2)
                       },
                       spd_nD = function(sigma, phi, w, D){
                         phisq <- phi^2
                         S = c(sigma)^2 * sqrt(2*pi)^D * (apply(phi,1,prod) * exp(-0.5*(phisq%*%w)))
                         return(S)
                       },
                       phi_nD = function(L, m, x){
                         fi <- t(1/(sqrt(L))* sin(m * (pi*(t(x)+L) * (1/(2*L))) ))
                         fi1 <- apply(fi,1,prod)
                         return(fi1)
                       },
                       transform = function(x){
                         x_grid <- x
                         for(i in 1:ncol(x_grid)){
                           x_grid[,i] <- 2*(x_grid[,i]-self$par_lower[i])/(self$par_upper[i]-self$par_lower[i]) - 1
                         }
                         return(x_grid)
                       },
                       back_transform = function(x){
                         x_grid <- x
                         for(i in 1:ncol(x_grid)){
                           x_grid[,i] <- (self$par_upper[i]-self$par_lower[i])*(x_grid[,i] + 1)/2 + self$par_lower[i]
                         }
                         return(x_grid)
                       },
                       sample_values = function(){
                         x <- rep(NA, length(self$par_upper))
                         for(i in 1:length(x)){
                           x[i] <- runif(1,self$par_lower[i],self$par_upper[i])
                           if(self$par_discrete[i])x[i] <- round(x[i],0)
                         }
                         return(x)
                       },
                       sigmae = NULL,
                       phi = NULL,
                       intercept = NULL,
                       lvec = NULL,
                       model = NULL,
                       type = NULL,
                       a_vector = NULL,
                       a_mat = NULL,
                       a_mat_prior = NULL,
                       avec_prior = NULL,
                       avec_conc = NULL,
                       as_d = NULL
                     ))
