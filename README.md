# AdaptSim

This is a prototype R package for running "adaptive simulations". The package defines the R6 Adapt class, which takes data simulation and model fitting functions as 
arguments. It then has two main functions, one for running a block of simulations, and another for running an adaptive step to estimate a model and sample new parameter 
values for the next block of simulations. Note, the documentation is not written yet for this package. 

## Example: Type I error with a linear mixed model
We provide an example of estimating the type I error for a linear mixed model for a cluster randomised trial, varying the number of clusters and the intraclass 
correlation coefficient. We aim to replicate some of the results in [Leyrat et al (2018)](https://doi.org/10.1093/ije/dyx169). 

### Data simulation and model fitting
First we define a function that simulates data sets, which takes as an argument a vector with the values of the iteration of the simulation. 
We use the `glmmrBase` package for data simulation of 
GLMMs:
```
data_sim <- function(x){
  df <- nelder(as.formula(paste0("~cl(",x[1],") > ind(10)")))
  df$int <- 0
  df[df$cl > x[1]/2,'int'] <- 1

  mod <- Model$new(
    covariance = list(
      formula = ~(1|gr(cl)),

      parameters = c(sqrt(x[2]))
    ),
    mean = list(
      formula = ~ int,
      parameters = rep(0,2) # change the parameters as required
    ),
    data = df,
    family = gaussian(),
    verbose = FALSE,
    var_par = sqrt(1-x[2])
  )
  data <- mod$sim_data(type = "data")
  return(data)
}
```
then we define a function that estimates the model of interest and retruns a named vector of the statistics of interest. We use `lme4` for model fitting:
```
fit_model <- function(data=data){
  fit1 <- lme4::lmer(y~int+(1|cl),data=data)#tryCatch(suppressWarnings(lme4::lmer(y~int+(1|cl),data=data)),error=function(i)NA)

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
```

### Create a new Adapt class object
Creation of a new Adapt class object is relatively simple. We provide the names of the two functions defined above, the upper and lower limits of the 
ranges of the varying parameters in the simulation (here, the number of clusters and the ICC), the inital number of simulations (we choose 10,000), and 
the number of basis functions for the adaptive step:
```
mod <- adapt$new(data_fn = "data_sim",
                 fit_fn = "fit_model",
                 par_lower = c(4,0),
                 par_upper = c(40,0.1),
                 par_discrete = c(TRUE,FALSE),
                 n = 10000,
                 m = c(5,10))
```
We also indicate that the first parameter (the number of clusters) is discrete, which forces the class to round parameter values to the nearest integer. We set the 
priors for the adaptive Bayesian model (although, the defaults usually work well):
```
mod$set_priors(lengthscale = matrix(c(0,0,0.5,0.5),ncol=2),
               fvar = c(0,0.5))
```

### Run the simulations
We can run the simulations in parallel by providing a cluster to the simulation function. We must export the relevant libraries to the cluster as well when creating
it. To create the cluster:
```
cl <- parallel::makeCluster(parallel::detectCores()-1)
parallel::clusterEvalQ(cl,library(glmmrBase))
parallel::clusterEvalQ(cl,library(lme4))
```
Then, to run the simualtions we simply call:
```
mod$sim(cl=cl)
```
The results are stored in the object and can be accessed if needed using `mod$sim_output`. 

### Run the adaptive step
We aim to estimate the type I error. So, our model is a binomial model where the outcome is whether the p-value was above or below the threshold 0.05 here, which 
is returned in the `fit_model` function as `p1`. We also need to indicate how many new simulation parameter values to draw, we choose 10,000 again. By default,
the function will generate a lattice over the domain of the simulation parameters at which to predict the value of the statistic of interest. Then new 
values are drawn by sampling uniformly within these cells, with the number of draws proportional to the adaptive statistic. If we instead want to specify a lattice
so that predictions are made at desirable values, then we can provide these values to the function. To replicate the results of Leyrat et al, we want to extract 
predictions at an ICC of 0.001 and 0.05 and for all values of the number of clusters. So we first define the data frame `samp_df` below, which is passed to the 
adaptive function in this and subsequent steps:
```
samp_df <- expand.grid(n = 4:40, icc = c(seq(0.001,0.05,length.out=8),seq(0.05,0.1,length.out=8)))
mod$adapt(stat="p1",
          n_new = 10000,
          kappa = 1,
          model="binomial",
          samp_n = samp_df,
          sampling = 150,
          L=3.2)
```
We can then repeat the procedure by calling `mod$sim` and `mod$adapt` until we reach our goal. 

### Extracting predictions
We can extract the predictions after each adaptive step as they are stored as `mod$adapt_output`. We have generated plots to compare against Leyrat et al using the 
following `ggplot2` command
```
ggplot(data=mod$adapt_output[mod$adapt_output$Var2%in%c(0.001,0.05),],
       aes(x=Var1,y=mean))+
  geom_hline(yintercept=0.05,lty=1)+
  geom_hline(yintercept = max(0,0.05-mod2$kappa), lty=3)+
  geom_hline(yintercept = min(1,0.05+mod2$kappa), lty=3)+
  geom_ribbon(aes(ymin=lci,ymax=uci),fill="#268bd2",alpha=0.2)+
  geom_line(color="#dc322f")+
  facet_wrap(~factor(Var2))+
  ggthemes::theme_solarized_2()+
  labs(x="Number of clusters",y="Type I error")
```
