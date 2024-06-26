functions {
  vector spd_nD(real sigma, row_vector phi, matrix w, int D){
    row_vector[cols(w)] S;
    S = sigma^2 * sqrt(2*pi())^D * prod(phi) * exp(-0.5*((phi .* phi) * w));
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
  array[1] real L;
  int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  int<lower=1> Nsample; //number of observations per time period
  array[Nsample] int y;
  matrix[Nsample,D] x_grid; //prediction grid and observations
  matrix[1,M_nD] lambda;
  array[M_nD,1] int indices;
  array[2] real intercept_prior;
  array[M_nD,2] real beta_prior;
  array[1,2] real lengthscale_prior;
  array[2] real fscale_prior;
  array[D] real avec_prior;
  real<lower=0> avec_conc;
  array[2] real sigma_prior;
}
parameters {
  vector[M_nD] beta;
  row_vector<lower=1e-05>[1] phi; //length scale
  real<lower=1e-6> sigma_e;
  real intercept;
  unit_vector[D] a;
  real<lower=1e-6> sigma;
}

transformed parameters{
  vector[Nsample] f;
  vector[M_nD] diagSPD;
  vector[Nsample] xa;
  matrix[Nsample,M_nD] PHI;

  xa = x_grid * a;
  for (m in 1:M_nD){
    PHI[,m] = phi_nD(L[1], indices[m,1], xa);
  }
  diagSPD = spd_nD(sigma_e, phi, lambda, 1);
  f = intercept + PHI * (diagSPD .* beta);

}
model{
  int grainsize = 1;
  phi ~ normal(lengthscale_prior[1,1],lengthscale_prior[1,2]);
  sigma_e ~ normal(fscale_prior[1],fscale_prior[2]);
  intercept ~ normal(intercept_prior[1],intercept_prior[2]);
  a ~ von_mises(avec_prior,avec_conc);
  sigma ~ gamma(sigma_prior[1],sigma_prior[2]);
  target += reduce_sum(partial_gauss_lpdf,to_array_1d(beta),grainsize,beta_prior[,1],beta_prior[,2]);
  target += reduce_sum(partial_sum_lpdf,y,grainsize,f,sigma);
}

