functions {
  vector spd_nD(real sigma, row_vector phi, matrix w, int D){
    row_vector[cols(w)] S;
    S = sigma^2 * sqrt(2*pi())^D * prod(phi) * exp(-0.5*((phi .* phi) * w));
    return sqrt(S');
  }
 matrix phi_nD(data array[] real L, data array[,] int m, matrix x, data int M_nD) {
    int c = cols(x);
    int r = rows(x);
    vector[r] xcol;
    matrix[r,M_nD] phi = rep_matrix(1.0,r,M_nD);
    real sqrtl;
    for(i in 1:c){
      xcol = (x[,i]+L[i])*pi()/(2*L[i]);
      sqrtl = 1/sqrt(L[i]);
      phi *= pow(sqrtl,M_nD);
      for(j in 1:M_nD){
        phi[,j] = phi[,j] .* sin(m[j,i]*xcol);
      }
    }
    return phi;
  }
  real partial_sum_lpdf(array[] real y,int start, int end, vector mu,real sigma){
    return normal_lpdf(y[start:end]|mu[start:end], sigma);
  }
  real partial_gauss_lpdf(array[] real y,int start, int end, array[] real mu, array[] real sigma){
    return normal_lpdf(y[start:end]|mu[start:end],sigma[start:end]);
  }
}
data {
  int<lower=2> D; //number of dimensions
  int<lower=2,upper=D> d; // subspace number of dimensions
  array[d] real L;
  int<lower=1> M_nD; //total basis functions m1*m2*...*md
  int<lower=1> Nsample; //number of observations
  array[Nsample] real y;
  matrix[Nsample,D] x_grid; //prediction grid and observations
  matrix[d,M_nD] lambda;
  array[M_nD,d] int indices;
  array[2] real intercept_prior;
  array[M_nD,2] real beta_prior;
  array[d,2] real lengthscale_prior;
  array[2] real fscale_prior;
  array[d*D] real a_mat_prior_mean;
  array[d*D] real a_mat_prior_sd;
  array[2] real sigma_prior;
}
parameters {
  vector[M_nD] beta;
  row_vector<lower=1e-05>[d] phi; //length scale
  real<lower=1e-6> sigma_e;
  real intercept;
  matrix[D,d] Amat;
  real<lower=1e-6> sigma;
}

transformed parameters{
  vector[Nsample] f;
  vector[M_nD] diagSPD;
  matrix[D,d] A = qr_thin_Q(Amat);
  matrix[Nsample,M_nD] PHI;
  matrix[Nsample,d] XA;

  A = qr_thin_Q(Amat);
  XA = x_grid * A;
  PHI = phi_nD(L, indices, XA, M_nD);
  diagSPD = spd_nD(sigma_e, phi, lambda, 1);
  f = intercept + PHI * (diagSPD .* beta);

}
model{
  int grainsize = 1;
  for(i in 1:d)phi[i] ~ normal(lengthscale_prior[i,1],lengthscale_prior[i,2]);
  sigma_e ~ normal(fscale_prior[1],fscale_prior[2]);
  intercept ~ normal(intercept_prior[1],intercept_prior[2]);
  for(i in 1:(d*D))to_vector(Amat)[i] ~ normal(a_mat_prior_mean[i],a_mat_prior_sd[i]);
  sigma ~ gamma(sigma_prior[1], sigma_prior[2]);
  target += reduce_sum(partial_gauss_lpdf,to_array_1d(beta),grainsize,beta_prior[,1],beta_prior[,2]);
  target += reduce_sum(partial_sum_lpdf,y,grainsize,f,sigma);
}

