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
  real partial_sum_lpmf(array[] int y,int start, int end, vector mu){
    return bernoulli_logit_lpmf(y[start:end]|mu[start:end]);
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
  array[Nsample] int y;
  matrix[Nsample,D] x_grid; //prediction grid and observations
  matrix[d,M_nD] lambda;
  array[M_nD,d] int indices;
  array[2] real intercept_prior;
  array[M_nD,2] real beta_prior;
  array[d,2] real lengthscale_prior;
  array[2] real fscale_prior;
  array[d*D] real a_mat_prior_mean;
  array[d*D] real a_mat_prior_sd;
}
parameters {
  vector[M_nD] beta;
  row_vector<lower=1e-05>[d] phi; //length scale
  real<lower=1e-6> sigma_e;
  real intercept;
  matrix[D,d] Amat;
}

transformed parameters{
  vector[Nsample] f;
  vector[M_nD] diagSPD;
  matrix[D,d] A = qr_thin_Q(Amat);
  matrix[Nsample,M_nD] PHI;

  for (m in 1:M_nD){
    PHI[,m] = phi_nD(L, indices[m,], x_grid * A);
  }
  diagSPD = spd_nD(sigma_e, phi, lambda, 1);
  f = intercept + PHI * (diagSPD .* beta);

}
model{
  int grainsize = 1;
  for(i in 1:d)phi[i] ~ normal(lengthscale_prior[i,1],lengthscale_prior[i,2]);
  sigma_e ~ normal(fscale_prior[1],fscale_prior[2]);
  intercept ~ normal(intercept_prior[1],intercept_prior[2]);
  for(i in 1:(d*D))to_vector(Amat)[i] ~ normal(a_mat_prior_mean[i],a_mat_prior_sd[i]);
  //to_vector(Amat) ~ std_normal();
  target += reduce_sum(partial_gauss_lpdf,to_array_1d(beta),grainsize,beta_prior[,1],beta_prior[,2]);
  target += reduce_sum(partial_sum_lpmf,y,grainsize,f);
}

