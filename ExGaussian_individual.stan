//
// This Stan program defines an Exponentially modified Gaussian
// model of Stop Signal Response Time data
// for a single participant, for a single run
//
// Learn more about the project at https://github.com/josefmana/dbs_combSTIM
//

functions {
/*  my_exp_mod_normal_lccdf(real y, real mu, real sigma, real lambda) {
    real term1 = normal_lccdf(y | mu, sigma);
    real term2a = ((lambda / 2) * (2*mu-lambda * sigma^2 - 2 * y))
    //erfc(x) = 2(1 - phi(x*sqrt(2)))
    //log(erfc(x)) = log(2) + std_normal_lccdf(x * sqrt(2))
    real term2b = std_normal_lccdf((mu + lambda * sigma^2 -  y)/ sigma)
    return(log_sum_exp()
  }
  */
  // define function computing integrand for numerical integration in succesful stop trials
  real integrand(real y, real yc, array[] real theta, array[] real x_r, array[] int x_i) {
    
    // define parametes
    real stan_muGO = theta[1];
    real stan_sigmaGO = theta[2];
    real stan_lambdaGO = theta[3];
    real stan_muSTOP = theta[4];
    real stan_sigmaSTOP = theta[5];
    real stan_lambdaSTOP = theta[6];
    real ssd = x_r[1];
    int integrand_variant = x_i[1];
    int debug_integrand = x_i[2];
    
    // compute the integrand
    real value;
    
    if (integrand_variant == 0) {
      
      real log_term1 = exp_mod_normal_lccdf(y | stan_muGO, stan_sigmaGO, stan_lambdaGO);
      real log_term2 = exp_mod_normal_lpdf(y - ssd | stan_muSTOP, stan_sigmaSTOP, stan_lambdaSTOP);
      value = exp(log_term1 + log_term2);
      
      if( debug_integrand && ( is_nan(value) || is_inf(value) ) ) {
        print( "T1: ", log_term1,", T2: ", log_term2, ", y: ", y, ", muGO: ", stan_muGO, ", sigmaGO: ", stan_sigmaGO,", lambdaGO: ", stan_lambdaGO );
        print(
          "phi: ", normal_cdf(y | stan_muGO, stan_sigmaGO),
          ", exp: ", exp( (stan_lambdaGO / 2) * (2 * stan_muGO - stan_lambdaGO * stan_sigmaGO^2 - 2 * y) ),
          ", erfc: ", erfc( (stan_muGO + stan_lambdaGO * stan_sigmaGO^2 -  y) / (sqrt(2)*stan_sigmaGO) ) );
        print( "muSTOP: ", stan_muSTOP, ", sigmaSTOP", stan_sigmaSTOP, ", lambdaSTOP: ", stan_lambdaSTOP );

      }

    } else if(integrand_variant == 1) {
      
      real term1 = (1 - exp_mod_normal_cdf( y | stan_muGO, stan_sigmaGO, stan_lambdaGO ) );
      real log_term2 = exp_mod_normal_lpdf( y - ssd | stan_muSTOP, stan_sigmaSTOP, stan_lambdaSTOP );
      value = term1 * exp(log_term2);    
      
      if( debug_integrand && ( is_nan(value) || is_inf(value) ) ) {
        print( "T1:", term1, ", T2:", log_term2, ", y:", y, ", muGO: ", stan_muGO, ", sigmaGO: ", stan_sigmaGO,", lambdaGO: ", stan_lambdaGO );
        print(
          "phi: ", normal_cdf(y | stan_muGO, stan_sigmaGO),
          "exp: ", exp( (stan_lambdaGO / 2) * (2 * stan_muGO - stan_lambdaGO * stan_sigmaGO^2 - 2 * y) ),
          ", erfc: ", erfc( (stan_muGO + stan_lambdaGO * stan_sigmaGO^2 -  y)/ ( sqrt(2) * stan_sigmaGO) )
        );
        print( "muSTOP: ", stan_muSTOP, ", sigmaSTOP: ", stan_sigmaSTOP, ", lambdaSTOP: ", stan_lambdaSTOP );
      }
    }
    
    return value;

  }

}

data {

  int<lower=1> N_go;  // total number of observations in GO condition
  vector[N_go] Y_go;  // response variable in GO condition
  int<lower=1> N_sr;  // total number of observations in STOP-RESPOND data
  vector[N_sr] Y_sr;  // response variable in STOP-RESPOND data
  int<lower=1> N_na;  // total number of observations in SUCCESSFUL STOP data
  
   // SSDs
  vector[N_sr] SSD_sr; // vector od SSDs in STOP-RESPOND data
  vector[N_na] SSD_na; // vector od SSDs in SUCCESSFUL STOP data
  
  // prior specifications
  // GO process
  vector[2] prior_muGO;
  vector[2] prior_sigmaGO;
  vector[2] prior_lambdaGO;
  //
  // prior specifications
  // STOP process
  vector[2] prior_muSTOP;
  vector[2] prior_sigmaSTOP;
  vector[2] prior_lambdaSTOP;
  
  // integrand variant & debugging option
  int<lower=0, upper=1> integration_mode;
  int<lower=0, upper=1> integrand_variant;
  int<lower=0, upper=1> debug_integrand;
}

transformed data {
  
  // dummy data for the numerical integration
  array[0] real x_r;
  array[0] int x_i;
  real integration_upper = integration_mode ? 10 : positive_infinity();
  
}

parameters {
  
  // intercepts
  real Int_mu_go;
  real Int_sigma_go;
  real Int_lambda_go;
  real Int_mu_stop;
  real Int_sigma_stop;
  real Int_lambda_stop;

}

transformed parameters {
  
  // prior contributions to the log posterior
  real lprior = 0;
  lprior += normal_lpdf( Int_mu_go | prior_muGO[1], prior_muGO[2] );
  lprior += normal_lpdf( Int_sigma_go | prior_sigmaGO[1], prior_sigmaGO[2] );
  lprior += normal_lpdf( Int_lambda_go | prior_lambdaGO[1], prior_lambdaGO[2] );
  lprior += normal_lpdf( Int_mu_stop | prior_muSTOP[1], prior_muSTOP[2] );
  lprior += normal_lpdf( Int_sigma_stop | prior_sigmaSTOP[1], prior_sigmaSTOP[2] );
  lprior += normal_lpdf( Int_lambda_stop | prior_lambdaSTOP[1], prior_lambdaSTOP[2] );

}

model {
  
  // likelihood
  //
  // measurement model for GO trials
  //
  // initialise parameters
  real muGO_go = Int_mu_go;
  real sigmaGO_go = Int_sigma_go;
  real lambdaGO_go = Int_lambda_go;
  //
  // exponentiate to get out of the log-scale
  muGO_go = exp(muGO_go);
  sigmaGO_go = exp(sigmaGO_go);
  lambdaGO_go = exp(lambdaGO_go);
  //
  // compute Stan ExGaussian parameters
  real stan_muGO_go = muGO_go - lambdaGO_go;
  real stan_sigmaGO_go = sigmaGO_go;
  real stan_lambdaGO_go = inv(lambdaGO_go);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf(Y_go | stan_muGO_go, stan_sigmaGO_go, stan_lambdaGO_go);
  //
  // measurement model for STOP-RESPONSE trials
  //
  // initialise parameters
  real muGO_sr = Int_mu_go;
  real sigmaGO_sr = Int_sigma_go;
  real lambdaGO_sr = Int_lambda_go;
  real muSTOP_sr = Int_mu_stop;
  real sigmaSTOP_sr = Int_sigma_stop;
  real lambdaSTOP_sr = Int_lambda_stop;
  //
  // exponentiate to get out of the log-scale
  muGO_sr = exp(muGO_sr);
  sigmaGO_sr = exp(sigmaGO_sr);
  lambdaGO_sr = exp(lambdaGO_sr);
  muSTOP_sr = exp(muSTOP_sr);
  sigmaSTOP_sr = exp(sigmaSTOP_sr);
  lambdaSTOP_sr = exp(lambdaSTOP_sr);
  //
  // compute Stan ExGaussian parameters
  real stan_muGO_sr = muGO_sr - lambdaGO_sr;
  real stan_sigmaGO_sr = sigmaGO_sr;
  real stan_lambdaGO_sr = inv(lambdaGO_sr);
  real stan_muSTOP_sr = muSTOP_sr - lambdaSTOP_sr;
  real stan_sigmaSTOP_sr = sigmaSTOP_sr;
  real stan_lambdaSTOP_sr = inv(lambdaSTOP_sr);
  //
  // add to the log likelihood
  target += exp_mod_normal_lpdf(Y_sr | stan_muGO_sr, stan_sigmaGO_sr, stan_lambdaGO_sr);
  target += exp_mod_normal_lccdf(Y_sr - SSD_sr | stan_muSTOP_sr, stan_sigmaSTOP_sr, stan_lambdaSTOP_sr);
  //
  // measurement model for SUCCESSFUL-STOP trials
  real muGO_na = Int_mu_go;
  real sigmaGO_na = Int_sigma_go;
  real lambdaGO_na = Int_lambda_go;
  real muSTOP_na = Int_mu_stop;
  real sigmaSTOP_na = Int_sigma_stop;
  real lambdaSTOP_na = Int_lambda_stop;
  //
  // exponentiate to get out of the log-scale
  muGO_na = exp(muGO_na);
  sigmaGO_na = exp(sigmaGO_na);
  lambdaGO_na = exp(lambdaGO_na);
  muSTOP_na = exp(muSTOP_na);
  sigmaSTOP_na = exp(sigmaSTOP_na);
  lambdaSTOP_na = exp(lambdaSTOP_na);
  //
  // compute Stan ExGaussian parameters
  real stan_muGO_na = muGO_na - lambdaGO_na;
  real stan_sigmaGO_na = sigmaGO_na;
  real stan_lambdaGO_na = inv(lambdaGO_na);
  real stan_muSTOP_na = muSTOP_na - lambdaSTOP_na;
  real stan_sigmaSTOP_na = sigmaSTOP_na;
  real stan_lambdaSTOP_na = inv(lambdaSTOP_na);
  //
  // add to the log likelihood
  for (n in 1:N_na) {
    target += log( integrate_1d( integrand, 0, integration_upper, {stan_muGO_na, stan_sigmaGO_na, stan_lambdaGO_na, stan_muSTOP_na, stan_sigmaSTOP_na, stan_lambdaSTOP_na}, {SSD_na[n]}, {integrand_variant, debug_integrand}, .001 ) );
  }
  // add priors including constants
  target += lprior;
  
}
