// matrix normal model to estimate covariances only, with negative binomial
//   likelihood
data {
  
  // main indices
  int<lower=1> N;                     // num observations
  int<lower=1> K;                     // num covariates
  int<lower=1> nsp;                   // number of species
  int<lower=1> nclass;                // number of size classes
  int<lower=1> nq;                    // nsp * nclass
  
  // response, covariate, and effort data  
  int<lower=0> nflat;
  int<lower=0> yflat[nflat];
  row_vector[N] log_effort;               // effort data
  
  // random effects
  int<lower=1> nriver;                // number of rivers
  int<lower=1> nreach;                // number of reaches
  int<lower=1> nsite;                 // number of sites
  int<lower=1> nyear;                 // number of years
  int<lower=1> ngear;
  int<lower=1,upper=nriver> river[N]; // river identifier for random effects
  int<lower=1,upper=nreach> reach[N]; // reach identifier for random effects
  int<lower=1,upper=nsite> site[N];   // site identifier for random effects
  int<lower=1,upper=nyear> year[N];   // year identifier for random effects
  int<lower=1,upper=ngear> gear[N];
  
  // scale of main effects and overall variance
  real<lower=0> main_scale;           // scale of main effects (alpha)
  real<lower=0> sigma_scale;          // scale of residuals
  real<lower=0> sigma_scale2;         // scale of random effects

}

transformed data {
  real log_half = -0.693147180559945286;
}

parameters {
  
  // regression coefficients with hierarchical structure
  vector<lower=0>[nsp] sigma_alpha;
  vector[nsp] z_alpha_mean;
  matrix[nsp, nclass] z_alpha;

  // random effects
  matrix[nq, nriver] gamma_river;
  matrix[nq, nreach] gamma_reach;
  matrix[nq, nsite] gamma_site;
  matrix[nq, nyear] gamma_year;
  matrix[nq, ngear] gamma_gear;
  
  // prior for correlation matrix (Cholesky decomposition)
  cholesky_factor_corr[nq] L;

  // random variates used to define matrix/multinormal  
  matrix[nq, N] eps;
  
  // priors for variances (added to correlation matrices)
  vector<lower=0>[nq] z_sigma;

  // variances for random terms
  real sigma_main_river;
  real sigma_main_reach;
  real sigma_main_site;
  real sigma_main_year;
  real sigma_main_gear;
  vector<lower=0>[nq] sigma_river;
  vector<lower=0>[nq] sigma_reach;
  vector<lower=0>[nq] sigma_site;
  vector<lower=0>[nq] sigma_year;
  vector<lower=0>[nq] sigma_gear;

  // over-dispersion parameter
  vector<lower=0>[nq] phi;

}

transformed parameters {
  vector[nq] alpha;
  matrix[nq, N] mu;
  vector<lower=0>[nq] sigma;
  matrix[nq, nq] W;
  vector[nflat] mu_flat;
  vector[nflat] phi_flat;

  // define covariance matrix from correlations and standard deviations
  sigma = sigma_scale * z_sigma;
  W = diag_pre_multiply(sigma, L);

  // rescale z-transformed covariate effects
  // with exchangeable priors within species (shared mean across all classes)
  alpha = to_vector(
    rep_matrix(sigma_alpha, nclass) .* z_alpha +
    rep_matrix(main_scale * z_alpha_mean, nclass)
  );

  // combine linear predictor
  mu = rep_matrix(alpha, N) +
    sigma_scale2 *
    (sigma_main_river * rep_matrix(sigma_river, N) .* gamma_river[, river] +
     sigma_main_reach * rep_matrix(sigma_reach, N) .* gamma_reach[, reach] +
     sigma_main_site * rep_matrix(sigma_site, N) .* gamma_site[, site] +
     sigma_main_year * rep_matrix(sigma_year, N) .* gamma_year[, year] +
     sigma_main_gear * rep_matrix(sigma_gear, N) .* gamma_gear[, gear]) +
    rep_matrix(log_effort, nq) +
    W * eps;

  // flatten mu and theta_zero for likelihood calcs
  mu_flat = to_vector(mu);
  phi_flat = to_vector(rep_matrix(phi, N));

}

model {
  
  // priors for fixed part of linear predictor
  z_alpha_mean ~ std_normal();
  to_vector(z_alpha) ~ std_normal();
  sigma_alpha ~ std_normal();

  // random effects
  to_vector(gamma_river) ~ std_normal();
  to_vector(gamma_reach) ~ std_normal();
  to_vector(gamma_site) ~ std_normal();
  to_vector(gamma_year) ~ std_normal();
  to_vector(gamma_gear) ~ std_normal();
  
  // and their variances
  sigma_main_river ~ std_normal();
  sigma_main_reach ~ std_normal();
  sigma_main_site ~ std_normal();
  sigma_main_year ~ std_normal();
  sigma_main_gear ~ std_normal();
  sigma_river ~ std_normal();
  sigma_reach ~ std_normal();
  sigma_site ~ std_normal();
  sigma_year ~ std_normal();
  sigma_gear ~ std_normal();
  
  // prior for variance
  z_sigma ~ std_normal();

  // prior for correlation matrix
  L ~ lkj_corr_cholesky(1.);

  // residuals
  to_vector(eps) ~ std_normal();
  
  // prior for phi
  target += student_t_lpdf(phi | 3.0, 0, 1) - log_half;

  // calculate likelihood of response (negative binomial)
  yflat ~ neg_binomial_2_log(mu_flat, phi_flat);

}

generated quantities {
  matrix[nq, nq] Sigma;

  // calculate correlation matrix from Cholesky decomposition
  Sigma = W * W';

}
