//
// This Stan program defines an Exponentially modified Gaussian
// model of Stop Signal Response Time data
// for a bunch of participants with hierarchical priors
//
// Learn more about the project at https://github.com/josefmana/dbs_combSTIM
//

functions {
  
  // define function computing integrand for numerical integration in succesful stop trials
  real integrand(real y, real yc, array[] real theta, array[] real x_r, array[] int x_i) {
    
    // definde parametes
    real mu_go = theta[1];
    real sigma_go = theta[2];
    real lambda_go = theta[3];
    real mu_stop = theta[4];
    real sigma_stop = theta[5];
    real lambda_stop = theta[6];
    real ssd = theta[7];
    
    // compute the integrand
    return exp( exp_mod_normal_lccdf( y | mu_go, sigma_go, lambda_go ) ) * exp( exp_mod_normal_lpdf( y - ssd | mu_stop, sigma_stop, lambda_stop ) );
    
  }

}

data {

  int<lower=1> N_0_go;  // total number of observations in GO condition
  vector[N_0_go] Y_0_go;  // response variable in GO condition
  int<lower=1> N_0_sr;  // total number of observations in STOP-RESPOND data
  vector[N_0_sr] Y_0_sr;  // response variable in STOP-RESPOND data
  int<lower=1> N_0_na;  // total number of observations in SUCCESSFUL STOP data
  
   // SSDs
  vector[N_0_sr] SSD_0_sr; // vector od SSDs in STOP-RESPOND data
  vector[N_0_na] SSD_0_na; // vector od SSDs in SUCCESSFUL STOP data
  
  // data for participant-level parameters
  int<lower=1> K;  // number of participants
  array[N_0_go] int<lower=1> J_0_go;  // grouping indicator per observation
  array[N_0_sr] int<lower=1> J_0_sr;  // grouping indicator per observation
  array[N_0_na] int<lower=1> J_0_na;  // grouping indicator per observation
  
  // prior specifications
  //
  // GO process
  //
  // global intercepts
  vector[2] alpha_go_0_p;
  vector[2] beta_go_0_p;
  vector[2] gamma_go_0_p;
  //
  // participant-level standard deviations
  vector[2] tau_go_0_p;
  vector[2] zeta_go_0_p;
  vector[2] epsilon_go_0_p;
  //
  // STOP process
  //
  // global intercepts
  vector[2] alpha_stop_0_p;
  vector[2] beta_stop_0_p;
  vector[2] gamma_stop_0_p;
  //
  // participant-level standard deviations
  vector[2] tau_stop_0_p;
  vector[2] zeta_stop_0_p;
  vector[2] epsilon_stop_0_p;
  
}

transformed data {
  
  // dummy data for the numerical integration
  array[0] real x_r;
  array[0] int x_i;

}

parameters {
  
  // global intercepts
  real alpha_go_0;
  real beta_go_0;
  real gamma_go_0;
  real alpha_stop_0;
  real beta_stop_0;
  real gamma_stop_0;
  //
  // participant-level standard deviations
  real tau_go_0;
  real zeta_go_0;
  real epsilon_go_0;
  real tau_stop_0;
  real zeta_stop_0;
  real epsilon_stop_0;
  //
  // standardized participant-level parameters
  vector[K] x_go_0;
  vector[K] y_go_0;
  vector[K] z_go_0;
  vector[K] x_stop_0;
  vector[K] y_stop_0;
  vector[K] z_stop_0;

}

transformed parameters {
  
  // actual group-level parameters
  //
  // GO process
  vector[K] r_mu_go_0 = rep_vector(0.0, K);
  vector[K] r_sigma_go_0 = rep_vector(0.0, K);
  vector[K] r_lambda_go_0 = rep_vector(0.0, K);
  r_mu_go_0 += ( exp(tau_go_0) * (x_go_0) );
  r_sigma_go_0 += ( exp(zeta_go_0) * (y_go_0) );
  r_lambda_go_0 += ( exp(epsilon_go_0) * (z_go_0) );
  //
  // STOP process
  vector[K] r_mu_stop_0 = rep_vector(0.0, K);
  vector[K] r_sigma_stop_0 = rep_vector(0.0, K);
  vector[K] r_lambda_stop_0 = rep_vector(0.0, K);
  r_mu_stop_0 += ( exp(tau_stop_0) * (x_stop_0) );
  r_sigma_stop_0 += ( exp(zeta_stop_0) * (y_stop_0) );
  r_lambda_stop_0 += ( exp(epsilon_stop_0) * (z_stop_0) );
  //
  // prior contributions to the log posterior
  real lprior = 0;
  //
  // sampling from prior
  //
  // global intercepts
  lprior += normal_lpdf( alpha_go_0 | alpha_go_0_p[1], alpha_go_0_p[2] );
  lprior += normal_lpdf( beta_go_0 | beta_go_0_p[1], beta_go_0_p[2] );
  lprior += normal_lpdf( gamma_go_0 | gamma_go_0_p[1], gamma_go_0_p[2] );
  lprior += normal_lpdf( alpha_stop_0 | alpha_stop_0_p[1], alpha_stop_0_p[2] );
  lprior += normal_lpdf( beta_stop_0 | beta_stop_0_p[1], beta_stop_0_p[2] );
  lprior += normal_lpdf( gamma_stop_0 | gamma_stop_0_p[1], gamma_stop_0_p[2] );
  //
  // participant-level standard deviations
  lprior += normal_lpdf( tau_go_0 | tau_go_0_p[1], tau_go_0_p[2] );
  lprior += normal_lpdf( zeta_go_0 | zeta_go_0_p[1], zeta_go_0_p[2] );
  lprior += normal_lpdf( epsilon_go_0 | epsilon_go_0_p[1], epsilon_go_0_p[2] );
  lprior += normal_lpdf( tau_stop_0 | tau_stop_0_p[1], tau_stop_0_p[2] );
  lprior += normal_lpdf( zeta_stop_0 | zeta_stop_0_p[1], zeta_stop_0_p[2] );
  lprior += normal_lpdf( epsilon_stop_0 | epsilon_stop_0_p[1], epsilon_stop_0_p[2] );

}

model {
  
  // likelihood
  //
  // measurement model for GO trials
  //
  // initialise parameters
  vector[N_0_go] mu_go_go_0 = rep_vector(0.0, N_0_go);
  vector[N_0_go] sigma_go_go_0 = rep_vector(0.0, N_0_go);
  vector[N_0_go] lambda_go_go_0 = rep_vector(0.0, N_0_go);
  //
  // add global intercepts
  mu_go_go_0 += alpha_go_0;
  sigma_go_go_0 += beta_go_0;
  lambda_go_go_0 += gamma_go_0;
  //
  // shift by participant-level parameters
  for (n in 1:N_0_go) {
    mu_go_go_0[n] += r_mu_go_0[J_0_go[n]];
    sigma_go_go_0[n] += r_sigma_go_0[J_0_go[n]];
    lambda_go_go_0[n] += r_lambda_go_0[J_0_go[n]];
  }
  //
  // exponentiate to get out of the log-scale
  mu_go_go_0 = exp(mu_go_go_0);
  sigma_go_go_0 = exp(sigma_go_go_0);
  lambda_go_go_0 = exp(lambda_go_go_0);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf(Y_0_go | mu_go_go_0-lambda_go_go_0, sigma_go_go_0, inv(lambda_go_go_0) );
  //
  // measurement model for STOP-RESPONSE trials
  //
  // initialise parameters
  vector[N_0_sr] mu_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] sigma_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] lambda_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] mu_sr_stop_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] sigma_sr_stop_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] lambda_sr_stop_0 = rep_vector(0.0, N_0_sr);
  //
  // add global Intercepts
  mu_sr_go_0 += alpha_go_0;
  sigma_sr_go_0 += beta_go_0;
  lambda_sr_go_0 += gamma_go_0;
  mu_sr_stop_0 += alpha_stop_0;
  sigma_sr_stop_0 += beta_stop_0;
  lambda_sr_stop_0 += gamma_stop_0;
  //
  // shift by participant-level parameters
  for (n in 1:N_0_sr) {
    mu_sr_go_0[n] += r_mu_go_0[J_0_sr[n]];
    sigma_sr_go_0[n] += r_sigma_go_0[J_0_sr[n]];
    lambda_sr_go_0[n] += r_lambda_go_0[J_0_sr[n]];
    mu_sr_stop_0[n] += r_mu_stop_0[J_0_sr[n]];
    sigma_sr_stop_0[n] += r_sigma_stop_0[J_0_sr[n]];
    lambda_sr_stop_0[n] += r_lambda_stop_0[J_0_sr[n]];
  }
  //
  // exponentiate to get out of the log-scale
  mu_sr_go_0 = exp(mu_sr_go_0);
  sigma_sr_go_0 = exp(sigma_sr_go_0);
  lambda_sr_go_0 = exp(lambda_sr_go_0);
  mu_sr_stop_0 = exp(mu_sr_stop_0);
  sigma_sr_stop_0 = exp(sigma_sr_stop_0);
  lambda_sr_stop_0 = exp(lambda_sr_stop_0);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf( Y_0_sr | mu_sr_go_0-lambda_sr_go_0, sigma_sr_go_0, inv(lambda_sr_go_0) );
  target += exp_mod_normal_lccdf( Y_0_sr-SSD_0_sr | mu_sr_stop_0-lambda_sr_stop_0, sigma_sr_stop_0, inv(lambda_sr_stop_0) );
  //
  // measurement model for SUCCESSFUL-STOP trials
  //
  // initialise parameters
  vector[N_0_na] mu_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] sigma_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] lambda_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] mu_na_stop_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] sigma_na_stop_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] lambda_na_stop_0 = rep_vector(0.0, N_0_na);
  //
  // add global Intercepts
  mu_na_go_0 += alpha_go_0;
  sigma_na_go_0 += beta_go_0;
  lambda_na_go_0 += gamma_go_0;
  mu_na_stop_0 += alpha_stop_0;
  sigma_na_stop_0 += beta_stop_0;
  lambda_na_stop_0 += gamma_stop_0;
  //
  // shift by participant-level parameters
  for (n in 1:N_0_na) {
    mu_na_go_0[n] += r_mu_go_0[J_0_na[n]];
    sigma_na_go_0[n] += r_sigma_go_0[J_0_na[n]];
    lambda_na_go_0[n] += r_lambda_go_0[J_0_na[n]];
    mu_na_stop_0[n] += r_mu_stop_0[J_0_na[n]];
    sigma_na_stop_0[n] += r_sigma_stop_0[J_0_na[n]];
    lambda_na_stop_0[n] += r_lambda_stop_0[J_0_na[n]];
  }
  // exponentiate to get out of the log-scale
  mu_na_go_0 = exp(mu_na_go_0);
  sigma_na_go_0 = exp(sigma_na_go_0);
  lambda_na_go_0 = exp(lambda_na_go_0);
  mu_na_stop_0 = exp(mu_na_stop_0);
  sigma_na_stop_0 = exp(sigma_na_stop_0);
  lambda_na_stop_0 = exp(lambda_na_stop_0);
  //
  // add to the log likelihood
  for (n in 1:N_0_na) {
    target += log( integrate_1d( integrand, 0, 10, { mu_na_go_0[n]-lambda_na_go_0[n], sigma_na_go_0[n], inv(lambda_na_go_0[n]), mu_na_stop_0[n]-lambda_na_stop_0[n], sigma_na_stop_0[n], inv(lambda_na_stop_0[n]), SSD_0_na[n] }, x_r, x_i, .001 ) );
  }
  //
  // add priors including constants
  target += lprior;
  target += std_normal_lpdf(x_go_0);
  target += std_normal_lpdf(y_go_0);
  target += std_normal_lpdf(z_go_0);
  target += std_normal_lpdf(x_stop_0);
  target += std_normal_lpdf(y_stop_0);
  target += std_normal_lpdf(z_stop_0);

}
