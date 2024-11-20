data_trial <- function(t,k,m,type){
  if(type %in% c("st","st_m","inc","sw")){
    df <- nelder(formula(paste0("~cl(",t-1,")*t(",t,")")))
    df$int <- I(df$t > df$cl)*2
    if(type %in% c("st","st_m")){
      #df <- df[df$t > 1 & df$t < t,]
      df <- df[df$cl == df$t | df$cl == (df$t - 1),]
      #df$t <- df$t - 1
    }
    if(type == "inc"){
      df <- df[df$t > (df$cl - 2) & df$t < (df$cl + 3),]
    }
    if(k > 1){
      df0 <- df
      for(i in 1:k){
        df1 <- df0
        df1$cl <- df1$cl + max(df$cl)
        df <- rbind(df,df1)
      }
    }
    df$cl <- as.numeric(as.factor(df$cl))
  } else {
    df <- nelder(formula(paste0("~cl(",k*2,")*t(",t,")")))
    if(type == "par"){
      df$int <- I(df$cl > k)*1
    } else {
      df$int <- 0
      for(i in 1:t){
        if(i%%2 == 0){
          df[df$cl <= k & df$t == i, 'int'] <- 1
        } else {
          df[df$cl > k & df$t == i, 'int'] <- 1
        }
      }
    }
  }
  
  df <- df[rep(1:nrow(df),each = m),] 
  
  return(df)
}

df1 <- data_trial(9,1,8,"st")
df2 <- data_trial(9,1,10,"st_m")
df3 <- data_trial(9,30,10,"par")

model1 <- Model$new(
  formula = ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
  data = df1,
  family= gaussian()
)

model2 <- Model$new(
  formula = ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
  data = df2,
  family= gaussian()
)

model3 <- Model$new(
  formula = ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
  data = df3,
  family= gaussian()
)

M1 <- model1$information_matrix()
M2 <- model2$information_matrix()
M3 <- model3$information_matrix()
M1k <- model1$small_sample_correction("KR")
M2k <- model2$small_sample_correction("KR")
M3k <- model3$small_sample_correction("KR")

e1 <- eigen(M1)
e2 <- eigen(M2)
e3 <- eigen(M3)

y2 <- model2$sim_data()
model2$update_y(y2)
M2b <- model2$information_matrix(oim = TRUE)

eigen(solve(M3k$vcov_beta))$values




data_sim <- function(t,k,m,icc,cac,type){
  # df <- data.frame(cl = 1:x[1])
  # 
  # if(x[3]==0){
  #   sizes <- rep(x[2],x[1])
  # } else {
  #   sizes <- rep(0,x[1])
  #   while(any(sizes<=0)){
  #     sizes[sizes <= 0] <- round(rnorm(length(sizes[sizes <= 0]), x[2], x[3]*x[2]),0)
  #   }
  # }
  # 
  # 
  # df <- df[rep(1:x[1],sizes),,drop=FALSE]
  # df$int <- 0
  # df[df$cl > x[1]/2,'int'] <- 1
  # 
  df <- data_trial(t,k,m,type)
  mod <- Model$new(
    formula = ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
    covariance = c(icc*cac,icc*(1-cac)),
    mean = rep(0,2),
    data = df,
    family = gaussian(),
    var_par = 1-icc
  )
  
  df$y <- mod$sim_data()
  return(df)
}

cond_number <- function(t,k,m,icc,cac,type){
  df <- data_trial(t,k,m,type)
  model <- Model$new(
    formula = ~ int + factor(t) + (1|gr(cl)) + (1|gr(cl,t)),
    covariance = c(icc*cac,icc*(1-cac)),
    mean = rep(0,length(unique(df$t))+1),
    data = df,
    family = gaussian(),
    var_par = 1-icc
  )
  M <- model$information_matrix()
  e <- eigen(M)$values
  return(c(max(e)/min(e), max(e), min(e), prod(e), sum(e)))
}

cond_number(6,2,20,0.1,0.8,"sw")
#seq(0.001,0.1,length.out = 10)
test <- expand.grid(t = 5:10,k = 1:2,icc = c(0.001,0.05), cac = c(0.2,0.8), type = c("sw","st"), cn = NA, l = NA, s = NA, prod = NA, sum = NA)

for(i in 1:nrow(test)){
  res <- cond_number(test$t[i],test$k[i],20,test$icc[i],test$cac[i],test$type[i])
  test[i,6:10] <- res
}

ggplot(data = test, aes(x = t, y = s, color = factor(icc), lty = factor(cac)))+
  geom_line()+
  facet_grid(k~type)
