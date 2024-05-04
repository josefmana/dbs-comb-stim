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
  
  // data for participant-level parameters in control condition
  int<lower=1> K_0;  // number of participants
  int<lower=1> M_0;  // number of coefficients per level
  array[N_0] int<lower=1> J_0;  // grouping indicator per observation
  
  // data for participant-level parameters in experimental condition
  int<lower=1> K_1;  // number of participants
  int<lower=1> M_1;  // number of coefficients per level
  array[N_1] int<lower=1> J_1;  // grouping indicator per observation
  
  // participant-level predictor values
  vector[N0] Z_0_1; // control condition
  vector[N1] Z_1_1; // experimental condition

}

parameters {
  
  // temporary intercepts for control condition
  real Intercept_0;
  real InterceptSigma_0;
  real InterceptTau_0;
  
  // temporary intercepts for experimental condition
  real Intercept_1;
  real InterceptSigma_1;
  real InterceptTau_1;
  
  // participant level standard deviations and standardized parameters
  vector<lower=0>[M_0] sd_0;
  array[M_0] vector[N_0] z_0;
  vector<lower=0>[M_1] sd_1;
  array[M_1] vector[N_1] z_1;
  
}

transformed parameters {

  // control condition
  vector[N_0] r_0_1;  // actual group-level effects
  real lprior_0 = 0;  // prior contributions to the log posterior
  r_0_1 = ( sd_0[1] * (z_0[1]) );
  lprior_0 += normal_lpdf(Intercept_0 | 700, 300);
  lprior_0 += normal_lpdf(InterceptSigma_0 | 0, 2.5);
  lprior_0 += normal_lpdf(InterceptTau_0 | 1.7, 1.3);
  lprior_0 += normal_lpdf(sd_0 | 0, 500) - 1 * normal_lccdf(0 | 0, 500);
  
  // experimental condition
  vector[N_1] r_1_1;  // actual group-level effects
  real lprior_1 = 0;  // prior contributions to the log posterior
  r_1_1 = ( sd_1[1] * (z_1[1]) );
  lprior_1 += normal_lpdf(Intercept_1 | 700, 300);
  lprior_1 += normal_lpdf(InterceptSigma_1 | 0, 2.5);
  lprior_1 += normal_lpdf(InterceptTau_1 | 1.7, 1.3);
  lprior_1 += normal_lpdf(sd_0 | 0, 500) - 1 * normal_lccdf(0 | 0, 500);

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
  vector[N_0] tau_0 = rep_vector(0.0, N_0);
  vector[N_1] tau_1 = rep_vector(0.0, N_1);
  
  /// add measurement model for control condition
  mu_0 += Intercept_0;
  sigma_0 += InterceptSigma_0;
  tau_0 += InterceptTau_0;
  for (n in 1:N_0) {
    // add more terms to the linear predictor
    mu_0[n] += r_0_1[J_0[n]] * Z_0_1[n];
  }
  sigma_0 = exp(sigma_0);
  tau_0 = exp(tau_0);
  target += exp_mod_normal_lpdf(Y | mu_0 - tau_0, sigma_0, inv(tau_0) );
  
  /// add measurement model for experimental condition
  mu_1 += Intercept_1;
  sigma_1 += InterceptSigma_1;
  tau_1 += InterceptTau_1;
  for (n in 1:N_1) {
    // add more terms to the linear predictor
    mu_1[n] += r_1_1[J_0[n]] * Z_1_1[n];
  }
  sigma_1 = exp(sigma_1);
  tau_1 = exp(tau_1);
  target += exp_mod_normal_lpdf(Y | mu_1 - tau_1, sigma_1, inv(tau_1) );
  
  // add priors including constants
  target += lprior_0;
  target += std_normal_lpdf(z_0[1]);
  target += lprior_1;
  target += std_normal_lpdf(z_1[1]);
  
}

generated quantities {
  
  // intercepts for mu
  real b_Intercept_0 = Intercept_0;
  real b_Intercept_1 = Intercept_1;
  
  // intercepts for sigma
  real b_sigma_Intercept_0 = Intercept_sigma_0;
  real b_sigma_Intercept_1 = Intercept_sigma_1;
  
  // intercepts for tau
  real b_tau_Intercept_0 = Intercept_tau_0;
  real b_tau_Intercept_1 = Intercept_tau_1;

}