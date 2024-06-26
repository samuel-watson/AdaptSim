functions {
  vector spd_nD(real sigma, real phi, row_vector w, int D){
    row_vector[size(w)] S;
    S = sigma^2 * sqrt(2*pi())^D * phi * exp(-0.5*((phi * phi) * w));
    return sqrt(S');
  }
  vector phi_nD(real L, int m, vector x) {
    int r = rows(x);
    vector[r] fi1;
    fi1 = 1/sqrt(L)*sin(m*pi()*(x+L)/(2*L));
    return fi1;
  }
  real partial_sum_lpdf(array[] real y,int start, int end, vector mu,real sigma){
    return normal_lpdf(y[start:end]|mu[start:end], sigma);
  }
  real partial_gauss_lpdf(array[] real y,int start, int end, array[] real mu, array[] real sigma){
    return normal_lpdf(y[start:end]|mu[start:end],sigma[start:end]);
  }
}
data {
  int<lower=1> D; //number of dimensions
  array[D] real L;
  array[D] int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  int<lower=1> total_fn;
  int<lower=1> Nsample; //number of observations per time period
  array[Nsample] int y;
  matrix[Nsample,D] x_grid; //prediction grid and observations
  vector[total_fn] lambda;
  array[2] real intercept_prior;
  array[total_fn,2] real beta_prior;
  array[D,2] real lengthscale_prior;
  array[2] real fscale_prior;
  array[2] real sigma_prior;
}
transformed data {
  matrix[Nsample,total_fn] PHI;
  for(i in 1:D){
    int m = i == 1 ? 0 : sum(M_nD[1:(i-1)]);
    for(j in 1:M_nD[i]){
      PHI[,m + j] = phi_nD(L[i], j, x_grid[,i]);
    }
  }
}
parameters {
  vector[total_fn] beta;
  row_vector<lower=1e-05>[D] phi; //length scale
  real<lower=1e-6> sigma_e;
  real intercept;
  real<lower=1e-6> sigma;
}

transformed parameters{
  vector[Nsample] f;
  vector[total_fn] diagSPD;
  for(i in 1:D){
    int m = i == 1 ? 0 : sum(M_nD[1:(i-1)]);
    diagSPD[(m+1):(m+M_nD[i])] = spd_nD(sigma_e, phi[i], lambda[(m+1):(m+M_nD[i])]', 1);
  }
  f = intercept + PHI * (diagSPD .* beta);
}
model{
  int grainsize = 1;
  for(d in 1:D)phi[d] ~ normal(lengthscale_prior[d,1],lengthscale_prior[d,2]);
  sigma_e ~ normal(fscale_prior[1],fscale_prior[2]);
  intercept ~ normal(intercept_prior[1],intercept_prior[2]);
  sigma ~ gamma(sigma_prior[1], sigma_prior[2]);
  target += reduce_sum(partial_gauss_lpdf,to_array_1d(beta),grainsize,beta_prior[,1],beta_prior[,2]);
  target += reduce_sum(partial_sum_lpdf,y,grainsize,f,sigma);
}
