//
// This Stan program defines an Exponentially modified Gaussian
// model of Stop Signal Response Time data
// for a single participant,
// for experimental and control runs
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

  int<lower=1> N_0_go;  // total number of observations in control GO condition
  vector[N_0_go] Y_0_go;  // response variable in control GO condition
  int<lower=1> N_0_sr;  // total number of observations in control condition for STOP-SIGNAL data
  vector[N_0_sr] Y_0_sr;  // response variable in control condition for STOP-SIGNAL data
  int<lower=1> N_0_na;  // total number of observations in control condition for SUCCESSFUL STOP data
  //vector[N_0_na] Y_0_na;  // response variable in control condition for SUCCESSFUL STOP data
  
  int<lower=1> N_1_go;  // total number of observations in experimental GO condition
  vector[N_1_go] Y_1_go;  // response variable in experimental GO condition
  int<lower=1> N_1_sr;  // total number of observations in experimental condition for STOP-SIGNAL data
  vector[N_1_sr] Y_1_sr;  // response variable in experimental condition for STOP-SIGNAL data
  int<lower=1> N_1_na;  // total number of observations in experimental condition for SUCCESSFUL STOP data
  //vector[N_1_na] Y_1_na;  // response variable in experimental condition for SUCCESSFUL STOP data
  
   // SSDs
  vector[N_0_sr] SSD_0_sr;
  vector[N_1_sr] SSD_1_sr;
  vector[N_0_na] SSD_0_na;
  vector[N_1_na] SSD_1_na;
  
  // prior specifications
  // GO process
  vector[2] mu_go_0_p;
  vector[2] mu_go_1_p;
  vector[2] sigma_go_0_p;
  vector[2] sigma_go_1_p;
  vector[2] beta_go_0_p;
  vector[2] beta_go_1_p;
  
  // prior specifications
  // STOP process
  vector[2] mu_stop_0_p;
  vector[2] mu_stop_1_p;
  vector[2] sigma_stop_0_p;
  vector[2] sigma_stop_1_p;
  vector[2] beta_stop_0_p;
  vector[2] beta_stop_1_p;
  
}

transformed data {
  
  // dummy data for the numerical integration
  array[0] real x_r;
  array[0] int x_i;

}

parameters {
  
  // intercepts for control condition
  real Int_mu_go_0;
  real Int_sigma_go_0;
  real Int_beta_go_0;
  real Int_mu_stop_0;
  real Int_sigma_stop_0;
  real Int_beta_stop_0;
  
  // intercepts for experimental condition
  real Int_mu_go_1;
  real Int_sigma_go_1;
  real Int_beta_go_1;
  real Int_mu_stop_1;
  real Int_sigma_stop_1;
  real Int_beta_stop_1;

}

transformed parameters {
  
  // prior contributions to the log posterior
  real lprior = 0;
  
  // control condition
  lprior += normal_lpdf( Int_mu_go_0 | mu_go_0_p[1], mu_go_0_p[2] );
  lprior += normal_lpdf( Int_sigma_go_0 | sigma_go_0_p[1], sigma_go_0_p[2] );
  lprior += normal_lpdf( Int_beta_go_0 | beta_go_0_p[1], beta_go_0_p[2] );
  lprior += normal_lpdf( Int_mu_stop_0 | mu_stop_0_p[1], mu_stop_0_p[2] );
  lprior += normal_lpdf( Int_sigma_stop_0 | sigma_stop_0_p[1], sigma_stop_0_p[2] );
  lprior += normal_lpdf( Int_beta_stop_0 | beta_stop_0_p[1], beta_stop_0_p[2] );
  
   // experimental condition
  lprior += normal_lpdf( Int_mu_go_1 | mu_go_1_p[1], mu_go_1_p[2] );
  lprior += normal_lpdf( Int_sigma_go_1 | sigma_go_1_p[1], sigma_go_1_p[2] );
  lprior += normal_lpdf( Int_beta_go_1 | beta_go_1_p[1], beta_go_1_p[2] );
  lprior += normal_lpdf( Int_mu_stop_1 | mu_stop_1_p[1], mu_stop_1_p[2] );
  lprior += normal_lpdf( Int_sigma_stop_1 | sigma_stop_1_p[1], sigma_stop_1_p[2] );
  lprior += normal_lpdf( Int_beta_stop_1 | beta_stop_1_p[1], beta_stop_1_p[2] );

}

model {
  
  // likelihood for the CONTROL condition
  //
  // measurement model for GO trials
  //
  // initialise parameters
  vector[N_0_go] mu_go_go_0 = rep_vector(0.0, N_0_go);
  vector[N_0_go] sigma_go_go_0 = rep_vector(0.0, N_0_go);
  vector[N_0_go] beta_go_go_0 = rep_vector(0.0, N_0_go);
  //
  // add global Intercepts
  mu_go_go_0 += Int_mu_go_0;
  sigma_go_go_0 += Int_sigma_go_0;
  beta_go_go_0 += Int_beta_go_0;
  //
  // exponentiate to get out of the log-scale
  mu_go_go_0 = exp(mu_go_go_0);
  sigma_go_go_0 = exp(sigma_go_go_0);
  beta_go_go_0 = exp(beta_go_go_0);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf(Y_0_go | mu_go_go_0-beta_go_go_0, sigma_go_go_0, inv(beta_go_go_0) );
  //
  // measurement model for control condition STOP-RESPONSE trials
  //
  // initialise parameters
  vector[N_0_sr] mu_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] sigma_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] beta_sr_go_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] mu_sr_stop_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] sigma_sr_stop_0 = rep_vector(0.0, N_0_sr);
  vector[N_0_sr] beta_sr_stop_0 = rep_vector(0.0, N_0_sr);
  //
  // add global Intercepts
  mu_sr_go_0 += Int_mu_go_0;
  sigma_sr_go_0 += Int_sigma_go_0;
  beta_sr_go_0 += Int_beta_go_0;
  mu_sr_stop_0 += Int_mu_stop_0;
  sigma_sr_stop_0 += Int_sigma_stop_0;
  beta_sr_stop_0 += Int_beta_stop_0;
  //
  // exponentiate to get out of the log-scale
  mu_sr_go_0 = exp(mu_sr_go_0);
  sigma_sr_go_0 = exp(sigma_sr_go_0);
  beta_sr_go_0 = exp(beta_sr_go_0);
  mu_sr_stop_0 = exp(mu_sr_stop_0);
  sigma_sr_stop_0 = exp(sigma_sr_stop_0);
  beta_sr_stop_0 = exp(beta_sr_stop_0);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf( Y_0_sr | mu_sr_go_0-beta_sr_go_0, sigma_sr_go_0, inv(beta_sr_go_0) );
  target += exp_mod_normal_lccdf( Y_0_sr-SSD_0_sr | mu_sr_stop_0-beta_sr_stop_0, sigma_sr_stop_0, inv(beta_sr_stop_0) );
  //
  // measurement model for control condition SUCCESSFUL-STOP trials
  //
  // initialise parameters
  vector[N_0_na] mu_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] sigma_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] beta_na_go_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] mu_na_stop_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] sigma_na_stop_0 = rep_vector(0.0, N_0_na);
  vector[N_0_na] beta_na_stop_0 = rep_vector(0.0, N_0_na);
  //
  // add global Intercepts
  mu_na_go_0 += Int_mu_go_0;
  sigma_na_go_0 += Int_sigma_go_0;
  beta_na_go_0 += Int_beta_go_0;
  mu_na_stop_0 += Int_mu_stop_0;
  sigma_na_stop_0 += Int_sigma_stop_0;
  beta_na_stop_0 += Int_beta_stop_0;
  //
  // exponentiate to get out of the log-scale
  mu_na_go_0 = exp(mu_na_go_0);
  sigma_na_go_0 = exp(sigma_na_go_0);
  beta_na_go_0 = exp(beta_na_go_0);
  mu_na_stop_0 = exp(mu_na_stop_0);
  sigma_na_stop_0 = exp(sigma_na_stop_0);
  beta_na_stop_0 = exp(beta_na_stop_0);
  //
  // add to the log likelihood
  for (n in 1:N_0_na) {
    target += log( integrate_1d( integrand, 0, 10, { mu_na_go_0[n]-beta_na_go_0[n], sigma_na_go_0[n], inv(beta_na_go_0[n]), mu_na_stop_0[n]-beta_na_stop_0[n], sigma_na_stop_0[n], inv(beta_na_stop_0[n]), SSD_0_na[n] }, x_r, x_i, .001 ) );
  }
  //
  // likelihood for the EXPERIMENTAL condition
  //
  // measurement model for GO trials
  //
  // initialise parameters
  vector[N_1_go] mu_go_go_1 = rep_vector(0.0, N_1_go);
  vector[N_1_go] sigma_go_go_1 = rep_vector(0.0, N_1_go);
  vector[N_1_go] beta_go_go_1 = rep_vector(0.0, N_1_go);
  //
  // add global Intercepts
  mu_go_go_1 += Int_mu_go_1;
  sigma_go_go_1 += Int_sigma_go_1;
  beta_go_go_1 += Int_beta_go_1;
  //
  // exponentiate to get out of the log-scale
  mu_go_go_1 = exp(mu_go_go_1);
  sigma_go_go_1 = exp(sigma_go_go_1);
  beta_go_go_1 = exp(beta_go_go_1);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf(Y_1_go | mu_go_go_1-beta_go_go_1, sigma_go_go_1, inv(beta_go_go_1) );
  //
  // measurement model for control condition STOP-RESPONSE trials
  //
  // initialise parameters
  vector[N_1_sr] mu_sr_go_1 = rep_vector(0.0, N_1_sr);
  vector[N_1_sr] sigma_sr_go_1 = rep_vector(0.0, N_1_sr);
  vector[N_1_sr] beta_sr_go_1 = rep_vector(0.0, N_1_sr);
  vector[N_1_sr] mu_sr_stop_1 = rep_vector(0.0, N_1_sr);
  vector[N_1_sr] sigma_sr_stop_1 = rep_vector(0.0, N_1_sr);
  vector[N_1_sr] beta_sr_stop_1 = rep_vector(0.0, N_1_sr);
  //
  // add global Intercepts
  mu_sr_go_1 += Int_mu_go_1;
  sigma_sr_go_1 += Int_sigma_go_1;
  beta_sr_go_1 += Int_beta_go_1;
  mu_sr_stop_1 += Int_mu_stop_1;
  sigma_sr_stop_1 += Int_sigma_stop_1;
  beta_sr_stop_1 += Int_beta_stop_1;
  //
  // exponentiate to get out of the log-scale
  mu_sr_go_1 = exp(mu_sr_go_1);
  sigma_sr_go_1 = exp(sigma_sr_go_1);
  beta_sr_go_1 = exp(beta_sr_go_1);
  mu_sr_stop_1 = exp(mu_sr_stop_1);
  sigma_sr_stop_1 = exp(sigma_sr_stop_1);
  beta_sr_stop_1 = exp(beta_sr_stop_1);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf( Y_1_sr | mu_sr_go_1-beta_sr_go_1, sigma_sr_go_1, inv(beta_sr_go_1) );
  target += exp_mod_normal_lccdf( Y_1_sr-SSD_1_sr | mu_sr_stop_1-beta_sr_stop_1, sigma_sr_stop_1, inv(beta_sr_stop_1) );
  //
  // measurement model for control condition SUCCESSFUL-STOP trials
  //
  // initialise parameters
  vector[N_1_na] mu_na_go_1 = rep_vector(0.0, N_1_na);
  vector[N_1_na] sigma_na_go_1 = rep_vector(0.0, N_1_na);
  vector[N_1_na] beta_na_go_1 = rep_vector(0.0, N_1_na);
  vector[N_1_na] mu_na_stop_1 = rep_vector(0.0, N_1_na);
  vector[N_1_na] sigma_na_stop_1 = rep_vector(0.0, N_1_na);
  vector[N_1_na] beta_na_stop_1 = rep_vector(0.0, N_1_na);
  //
  // add global Intercepts
  mu_na_go_1 += Int_mu_go_1;
  sigma_na_go_1 += Int_sigma_go_1;
  beta_na_go_1 += Int_beta_go_1;
  mu_na_stop_1 += Int_mu_stop_1;
  sigma_na_stop_1 += Int_sigma_stop_1;
  beta_na_stop_1 += Int_beta_stop_1;
  //
  // exponentiate to get out of the log-scale
  mu_na_go_1 = exp(mu_na_go_1);
  sigma_na_go_1 = exp(sigma_na_go_1);
  beta_na_go_1 = exp(beta_na_go_1);
  mu_na_stop_1 = exp(mu_na_stop_1);
  sigma_na_stop_1 = exp(sigma_na_stop_1);
  beta_na_stop_1 = exp(beta_na_stop_1);
  //
  // add to the log likelihood
  for (n in 1:N_1_na) {
    target += log( integrate_1d( integrand, 0, 10, { mu_na_go_1[n]-beta_na_go_1[n], sigma_na_go_1[n], inv(beta_na_go_1[n]), mu_na_stop_1[n]-beta_na_stop_1[n], sigma_na_stop_1[n], inv(beta_na_stop_1[n]), SSD_1_na[n] }, x_r, x_i, .001 ) );
  }

  // add priors including constants
  target += lprior;
  
}
