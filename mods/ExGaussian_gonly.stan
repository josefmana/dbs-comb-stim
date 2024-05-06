//
// This Stan program defines an Exponentially modified Gaussian
// model of Stop Signal Response Time data including
// global (i.e., 'fixed effect') Intercepts
// and participant-level deviations separately
// for experimental and control runs
//
// Learn more about the project at https://github.com/josefmana/dbs_combSTIM
//

data {

  int<lower=1> N_0;  // total number of observations in control condition
  vector[N_0] Y_0;  // response variable in control condition
  
  int<lower=1> N_1;  // total number of observations in experimental condition
  vector[N_1] Y_1;  // response variable in experimental condition
  
  // data for participant-level parameters (shared across conditions)
  int<lower=1> K;  // number of participants
  int<lower=1> M;  // number of coefficients per level
  array[N_0] int<lower=1> J_0;  // grouping indicator per observation in control condition
  array[N_1] int<lower=1> J_1;  // grouping indicator per observation in experimental condition
  
  // participant-level predictor values
  vector[N_0] Z_0_1; // control condition
  vector[N_1] Z_1_1; // experimental condition
  
  // prior specifications
  vector[2] mu_0p;
  vector[2] mu_1p;
  vector[2] sigma_0p;
  vector[2] sigma_1p;
  vector[2] beta_0p;
  vector[2] beta_1p;
  vector[2] tau_mu_p;
  //vector[2] tau_sigma_p;
  //vector[2] tau_beta_p;
  //vector[2] sd_1p;

}

parameters {
  
  // intercepts for control condition
  real<lower=0> InterceptMu_0;
  real InterceptSigma_0;
  real InterceptBeta_0;
  
  // intercepts for experimental condition
  real<lower=0> InterceptMu_1;
  real InterceptSigma_1;
  real InterceptBeta_1;
  
  // participant level standard deviations and standardized parameters
  vector<lower=0>[M] tau_mu;
  //vector<lower=0>[M] tau_sigma;
  //vector<lower=0>[M] tau_beta;
  array[M] vector[K] z_mu;
  //array[M] vector[K] z_sigma;
  //array[M] vector[K] z_beta;

}

transformed parameters {

  // actual group-level effects
  vector[K] r_mu;
  r_mu = ( tau_mu[1] * (z_mu[1]) );
  //vector[K] r_sigma;
  //r_sigma = ( tau_sigma[1] * (z_sigma[1]) );
  //vector[K] r_beta;
  //r_beta = ( tau_beta[1] * (z_beta[1]) );
  
  real lprior = 0;  // prior contributions to the log posterior
  
  // control condition
  lprior += normal_lpdf( InterceptMu_0 | mu_0p[1], mu_0p[2] );
  lprior += normal_lpdf( InterceptSigma_0 | sigma_0p[1], sigma_0p[2] );
  lprior += normal_lpdf( InterceptBeta_0 | beta_0p[1], beta_0p[2] );
  
   // experimental condition
  lprior += normal_lpdf( InterceptMu_1 | mu_1p[1], mu_1p[2] );
  lprior += normal_lpdf( InterceptSigma_1 | sigma_1p[1], sigma_1p[2] );
  lprior += normal_lpdf( InterceptBeta_1 | beta_1p[1], beta_1p[2] );
  
  // participant-level
  lprior += normal_lpdf( tau_mu | tau_mu_p[1], tau_mu_p[2] ) - 1 * normal_lccdf( 0 | tau_mu_p[1], tau_mu_p[2] );
  //lprior += normal_lpdf( tau_sigma | tau_sigma_p[1], tau_sigma_p[2] ) - 1 * normal_lccdf( 0 | tau_sigma_p[1], tau_sigma_p[2] );
  //lprior += normal_lpdf( tau_beta | tau_beta_p[1], tau_beta_p[2] ) - 1 * normal_lccdf( 0 | tau_beta_p[1], tau_beta_p[2] );

}

model {
  
  // likelihood including constants
  
  /// initialize ExGaussian means
  vector[N_0] mu_0 = rep_vector(0.0, N_0);
  vector[N_1] mu_1 = rep_vector(0.0, N_1);
  
  /// initialize ExGaussian standard deviations
  vector[N_0] sigma_0 = rep_vector(0.0, N_0);
  vector[N_1] sigma_1 = rep_vector(0.0, N_1);
  
  // initialize ExGaussian rates
  vector[N_0] beta_0 = rep_vector(0.0, N_0);
  vector[N_1] beta_1 = rep_vector(0.0, N_1);
  
  /// add measurement model for control condition
  mu_0 += InterceptMu_0;
  sigma_0 += InterceptSigma_0;
  beta_0 += InterceptBeta_0;
  // add more terms to the linear predictor
  for (n in 1:N_0) {
    mu_0[n] += r_mu[J_0[n]] * Z_0_1[n];
    //sigma_0[n] += r_sigma[J_0[n]] * Z_0_1[n];
    //beta_0[n] += r_beta[J_0[n]] * Z_0_1[n];
  }
  sigma_0 = exp(sigma_0);
  beta_0 = exp(beta_0);
  target += exp_mod_normal_lpdf(Y_0 | mu_0 - beta_0, sigma_0, inv(beta_0) );
  
  /// add measurement model for experimental condition
  mu_1 += InterceptMu_1;
  sigma_1 += InterceptSigma_1;
  beta_1 += InterceptBeta_1;
  // add more terms to the linear predictor
  for (n in 1:N_1) {
    mu_1[n] += r_mu[J_1[n]] * Z_1_1[n];
    //sigma_1[n] += r_sigma[J_1[n]] * Z_1_1[n];
    //beta_1[n] += r_beta[J_1[n]] * Z_1_1[n];
  }
  sigma_1 = exp(sigma_1);
  beta_1 = exp(beta_1);
  target += exp_mod_normal_lpdf(Y_1 | mu_1 - beta_1, sigma_1, inv(beta_1) );
  
  // add priors including constants
  target += lprior;
  target += std_normal_lpdf(z_mu[1]);
  //target += std_normal_lpdf(z_sigma[1]);
  //target += std_normal_lpdf(z_beta[1]);
  
}
