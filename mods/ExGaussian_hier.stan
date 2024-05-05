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
  
  // data for participant-level parameters in control condition
  //int<lower=1> K_0;  // number of participants
  //int<lower=1> M_0;  // number of coefficients per level
  //array[N_0] int<lower=1> J_0;  // grouping indicator per observation
  
  // data for participant-level parameters in experimental condition
  //int<lower=1> K_1;  // number of participants
  //int<lower=1> M_1;  // number of coefficients per level
  //array[N_1] int<lower=1> J_1;  // grouping indicator per observation
  
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
  vector[2] tau_p;
  //vector[2] sd_1p;

}

parameters {
  
  // intercepts for control condition
  real InterceptMu_0;
  real InterceptSigma_0;
  real InterceptBeta_0;
  
  // intercepts for experimental condition
  real InterceptMu_1;
  real InterceptSigma_1;
  real InterceptBeta_1;
  
  // participant level standard deviations and standardized parameters
  vector<lower=0>[M] tau;
  array[M] vector[K] z;
  
  //vector<lower=0>[M_0] sd_0;
  //array[M_0] vector[K_0] z_0;
  //vector<lower=0>[M_1] sd_1;
  //array[M_1] vector[K_1] z_1;
  
}

transformed parameters {

  
  vector[K] r;  // actual group-level effects
  r = ( tau[1] * (z[1]) );
  
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
  lprior += normal_lpdf( tau | tau_p[1], tau_p[2] ) - 1 * normal_lccdf( 0 | tau_p[1], tau_p[2] );

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
  for (n in 1:N_0) {
    // add more terms to the linear predictor
    mu_0[n] += r[J_0[n]] * Z_0_1[n];
  }
  sigma_0 = exp(sigma_0);
  beta_0 = exp(beta_0);
  target += exp_mod_normal_lpdf(Y_0 | mu_0 - beta_0, sigma_0, inv(beta_0) );
  
  /// add measurement model for experimental condition
  mu_1 += InterceptMu_1;
  sigma_1 += InterceptSigma_1;
  beta_1 += InterceptBeta_1;
  for (n in 1:N_1) {
    // add more terms to the linear predictor
    mu_1[n] += r[J_0[n]] * Z_1_1[n];
  }
  sigma_1 = exp(sigma_1);
  beta_1 = exp(beta_1);
  target += exp_mod_normal_lpdf(Y_1 | mu_1 - beta_1, sigma_1, inv(beta_1) );
  
  // add priors including constants
  target += lprior;
  target += std_normal_lpdf(z[1]);
  
}

generated quantities {
  
  // intercepts for mu
  //real b_InterceptMu_0 = InterceptMu_0;
  //real b_InterceptMu_1 = InterceptMu_1;
  
  // intercepts for sigma
  //real b_sigma_InterceptMu_0 = InterceptSigma_0;
  //real b_sigma_InterceptMu_1 = InterceptSigma_1;
  
  // intercepts for beta
  //real b_beta_InterceptMu_0 = InterceptBeta_0;
  //real b_beta_InterceptMu_1 = InterceptBeta_1;

}
