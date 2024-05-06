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

  int<lower=1> N_0_go;  // total number of observations in control GO condition
  vector[N_0_go] Y_0_go;  // response variable in control GO condition
  int<lower=1> N_0_sr;  // total number of observations in control condition for STOP-SIGNAL data
  vector[N_0_sr] Y_0_sr;  // response variable in control condition for STOP-SIGNAL data
  int<lower=1> N_0_na;  // total number of observations in control condition for SUCCESSFUL STOP data
  vector[N_0_na] Y_0_na;  // response variable in control condition for SUCCESSFUL STOP data
  
  int<lower=1> N_1_go;  // total number of observations in experimental GO condition
  vector[N_1_go] Y_1_go;  // response variable in experimental GO condition
  int<lower=1> N_1_sr;  // total number of observations in experimental condition for STOP-SIGNAL data
  vector[N_1_sr] Y_1_sr;  // response variable in experimental condition for STOP-SIGNAL data
  int<lower=1> N_1_na;  // total number of observations in experimental condition for SUCCESSFUL STOP data
  vector[N_1_na] Y_1_na;  // response variable in experimental condition for SUCCESSFUL STOP data
  
  // SSD
  vector[N_0_sr] SSD_0_sr;
  vector[N_1_sr] SSD_1_sr;
  vector[N_0_na] SSD_0_na;
  vector[N_1_na] SSD_1_na;
  
  // data for participant-level parameters (shared across conditions)
  int<lower=1> K;  // number of participants
  int<lower=1> M;  // number of coefficients per level
  array[N_0_go] int<lower=1> J_0_go;  // grouping indicator per observation in control condition
  array[N_1_go] int<lower=1> J_1_go;  // grouping indicator per observation in experimental condition
  array[N_0_sr] int<lower=1> J_0_sr;  // grouping indicator per observation in control condition
  array[N_1_sr] int<lower=1> J_1_sr;  // grouping indicator per observation in experimental condition
  array[N_0_na] int<lower=1> J_0_na;  // grouping indicator per observation in control condition
  array[N_1_na] int<lower=1> J_1_na;  // grouping indicator per observation in experimental condition
  
  // participant-level predictor values
  vector[N_0_go] Z_0_go; // control condition
  vector[N_1_go] Z_1_go; // experimental condition
  vector[N_0_sr] Z_0_sr; // control condition
  vector[N_1_sr] Z_1_sr; // experimental condition
  vector[N_0_na] Z_0_na; // control condition
  vector[N_1_na] Z_1_na; // experimental condition
  
  // prior specifications
  vector<lower=0>[2] mu_go_0_p;
  vector<lower=0>[2] mu_go_1_p;
  vector[2] sigma_go_0_p;
  vector[2] sigma_go_1_p;
  vector[2] beta_go_0_p;
  vector[2] beta_go_1_p;
  vector[2] tau_go_p;
  
  vector<lower=0>[2] mu_stop_0_p;
  vector<lower=0>[2] mu_stop_1_p;
  vector[2] sigma_stop_0_p;
  vector[2] sigma_stop_1_p;
  vector[2] beta_stop_0_p;
  vector[2] beta_stop_1_p;
  vector[2] tau_stop_p;

}

parameters {
  
  // intercepts for control condition
  real mu_go_0;
  real sigma_go_0;
  real beta_go_0;
  real mu_stop_0;
  real sigma_stop_0;
  real beta_stop_0;
  
  // intercepts for experimental condition
  real mu_go_1;
  real sigma_go_1;
  real beta_go_1;
  real mu_stop_1;
  real sigma_stop_1;
  real beta_stop_1;
  
  // participant level standard deviations and standardized parameters
  vector<lower=0>[M] tau_go;
  vector<lower=0>[M] tau_stop;
  array[M] vector[K] z_go;
  array[M] vector[K] z_stop;
  
}

transformed parameters {

  // actual group-level effects
  vector[K] r_go;
  r_go = ( tau_go[1] * (z_go[1]) );
  vector[K] r_stop;
  r_stop = ( tau_stop[1] * (z_stop[1]) );
  
  real lprior = 0;  // prior contributions to the log posterior
  
  // control condition
  lprior += normal_lpdf( mu_go_0 | mu_go_0_p[1], mu_go_0_p[2] );
  lprior += normal_lpdf( sigma_go_0 | sigma_go_0_p[1], sigma_go_0_p[2] );
  lprior += normal_lpdf( beta_go_0 | beta_go_0_p[1], beta_go_0_p[2] );
  lprior += normal_lpdf( mu_stop_0 | mu_stop_0_p[1], mu_stop_0_p[2] );
  lprior += normal_lpdf( sigma_stop_0 | sigma_stop_0_p[1], sigma_stop_0_p[2] );
  lprior += normal_lpdf( beta_stop_0 | beta_stop_0_p[1], beta_stop_0_p[2] );
  
   // experimental condition
  lprior += normal_lpdf( mu_go_1 | mu_go_1_p[1], mu_go_1_p[2] );
  lprior += normal_lpdf( sigma_go_1 | sigma_go_1_p[1], sigma_go_1_p[2] );
  lprior += normal_lpdf( beta_go_1 | beta_go_1_p[1], beta_go_1_p[2] );
  lprior += normal_lpdf( mu_stop_1 | mu_stop_1_p[1], mu_stop_1_p[2] );
  lprior += normal_lpdf( sigma_stop_1 | sigma_stop_1_p[1], sigma_stop_1_p[2] );
  lprior += normal_lpdf( beta_stop_1 | beta_stop_1_p[1], beta_stop_1_p[2] );
  
  // participant-level
  lprior += normal_lpdf( tau_go | tau_go_p[1], tau_go_p[2] ) - 1 * normal_lccdf( 0 | tau_go_p[1], tau_go_p[2] );
  lprior += normal_lpdf( tau_stop | tau_stop_p[1], tau_stop_p[2] ) - 1 * normal_lccdf( 0 | tau_stop_p[1], tau_stop_p[2] );

}

model {
  
  // measurement model for control condition GO trials
  for (n in 1:N_0_go) {
    mu_go_0[n] += r[J_0_go[n]] * Z_0_go[n];
  }
  sigma_go_0 = exp(sigma_go_0);
  beta_go_0 = exp(beta_go_0);
  target += exp_mod_normal_lpdf(Y_0_go | mu_go_0 - beta_go_0, sigma_go_0, inv(beta_go_0) );
  
  // measurement model for control condition STOP-RESPONSE trials
  for (n in 1:N_0_sr) {
    mu_stop_sr[n] += r[J_0_sr[n]] * Z_0_sr[n];
  }
  sigma_stop_0 = exp(sigma_stop_0);
  beta_stop_0 = exp(beta_stop_0);
  target += exp_mod_normal_lpdf( Y_0_sr | mu_go_0 - beta_go_0, sigma_go_0, inv(beta_go_0) );
  target += exp_mod_normal_lccdf( Y_0_sr - SSD_0_sr | mu_stop_0 - beta_stop_0, sigma_stop_0, inv(beta_stop_0) )
  
  // measurement model for control condition STOP-SUCCESS trials
  
  
  
  // add measurement model for experimental condition
  for (n in 1:N_1) {
    // add more terms to the linear predictor
    mu_1[n] += r[J_0[n]] * Z_1_1[n];
  }
  sigma_1 = exp(sigma_1);
  beta_1 = exp(beta_1);
  target += exp_mod_normal_lpdf(Y_1 | mu_1 - beta_1, sigma_1, inv(beta_1) );
  //target += exp_mod_normal_lpdf(Y_1 - SSD | mu_1 - beta_1, sigma_1, inv(beta_1) );
  
  // add priors including constants
  target += lprior;
  target += std_normal_lpdf(z[1]);
  
}

generated quantities {

}
