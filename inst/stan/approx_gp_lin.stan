functions {
  vector spd_nD(real sigma, row_vector phi, matrix w, int D){
    row_vector[cols(w)] S;
    S = sigma^2 * sqrt(2*pi())^D * prod(phi) * exp(-0.5*((phi .* phi) * w));
    return sqrt(S');
  }
  vector phi_nD(array[] real L, array[] int m, matrix x) {
    int c = cols(x);
    int r = rows(x);

    matrix[r,c] fi;
    vector[r] fi1;
    for (i in 1:c){
      fi[,i] = 1/sqrt(L[i])*sin(m[i]*pi()*(x[,i]+L[i])/(2*L[i]));
    }
    fi1 = fi[,1];
    for (i in 2:c){
      fi1 = fi1 .* fi[,i];
    }
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
  int<lower=1> M_nD; //total basis functions m1*m2*...*mD
  int<lower=1> Nsample; //number of observations per time period
  int<lower=1> Npred;
  array[Nsample] real y;
  matrix[Nsample+Npred,D] x_grid; //prediction grid and observations
  array[M_nD,D] int indices;
  matrix[D,M_nD] lambda;
  array[2] real intercept_prior;
  array[M_nD,2] real beta_prior;
  array[D,2] real lengthscale_prior;
  array[2] real fscale_prior;
  array[2] real sigma_prior;
}
transformed data {
  matrix[Nsample,M_nD] PHI;

  for (m in 1:M_nD){
    PHI[,m] = phi_nD(L, indices[m,], x_grid);
  }
}

parameters {
  vector[M_nD] beta;
  row_vector<lower=1e-05>[D] phi; //length scale
  real<lower=1e-05> sigma;
  real<lower=0> sigma_e;
  real intercept;
}

transformed parameters{
  vector[Nsample] f;
  vector[M_nD] diagSPD;
  diagSPD = spd_nD(sigma_e, phi, lambda, D);
  f = intercept + PHI * (diagSPD .* beta);

}
model{
  int grainsize = 1;
  // for(i in 1:M_nD)beta[i] ~ normal(beta_prior[i,1],beta_prior[i,2]);
  for(d in 1:D)phi[d] ~ normal(lengthscale_prior[d,1],lengthscale_prior[d,2]);
  sigma_e ~ normal(fscale_prior[1],fscale_prior[2]);
  intercept ~ normal(intercept_prior[1],intercept_prior[2]);
  sigma ~ normal(sigma_prior[1],sigma_prior[2]);
  //y ~ normal(f,sigma);
  target += reduce_sum(partial_gauss_lpdf,to_array_1d(beta),grainsize,beta_prior[,1],beta_prior[,2]);
  target += reduce_sum(partial_sum_lpdf,y,grainsize,f[1:Nsample],sigma);
}

// generated quantities{
//   vector[Npred] y_grid_predict;
// 
//   for(i in (Nsample+1):(Nsample+Npred)){
//     y_grid_predict[i-Nsample] = intercept + f[i];
//   }
// }
