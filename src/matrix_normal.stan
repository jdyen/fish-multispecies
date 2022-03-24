// matrix normal model with simplified parameterisation of covariances
data {
  int<lower=1> N;                     // num observations
  int<lower=1> K;                     // num covariates
  int<lower=1> nsp;                   // number of species
  int<lower=1> nclass;                // number of size classes
  int<lower=1> nq;                    // nsp * nclass
  matrix[N, K] X;                     // design matrix of covariates (fixed)
  int<lower=0> y[N, nq];              // response
  int<lower=1> yp1[N, nq];            // because matrices and arrays don't like each other
  vector[N] log_effort;               // effort data
  int<lower=1> nriver;                // number of rivers
  int<lower=1> nreach;                // number of reaches
  int<lower=1> nsite;                 // number of sites
  int<lower=1> nyear;                 // number of years
  int<lower=1,upper=nriver> river[N]; // river identifier for random effects
  int<lower=1,upper=nreach> reach[N]; // reach identifier for random effects
  int<lower=1,upper=nsite> site[N];   // site identifier for random effects
  int<lower=1,upper=nyear> year[N];   // year identifier for random effects
  int<lower=1> prev_idx[N];           // index of previous year's survey
  int<lower=0> visited_prev[N];       // was site surveyed in previous year?
  int<lower=1> nmissing;              // number of obs with no survey in previous year
  real<lower=0> phi_scale;            // scale of dispersion parameter
  real<lower=0> main_scale;           // scale of main effects (alpha/beta)
  real<lower=0> sigma_scale;          // scale of residuals
  real<lower=0> sigma_scale2;         // scale of random effects
  real<lower=0> ar_scale;             // scale of AR processes
}

transformed data {
  int ystart;
  int yend;
  int yflat[N * nq];
  real log_y[N, nq];
  real log_y_shift[N, nq];

  // flatten response
  for (i in 1:N) {
    ystart = (i - 1) * nq + 1;
    yend = (i - 1) * nq + nq;
    yflat[ystart:yend] = y[i];
  }
  
  // define log-transformed version of (y + 1) and a size-shifted version of this
  //   for use in the AR model
  log_y = log(yp1);
  for (i in 1:nsp) {
    log_y_shift[, (i - 1) * nclass + 1] = log_y[, i * nclass];
    log_y_shift[, ((i - 1) * nclass + 2):(i * nclass)] = log_y[, ((i - 1) * nclass + 1):(i * nclass - 1)];
  }

}

parameters {
  
  // overall mean
  row_vector[nq] z_alpha;
  
  // z-scaled regression coefficients in Q space
  matrix[K, nq] z_beta;

  // AR parameter
  row_vector[nq] z_rho;
  row_vector[nq] z_tau;
  
  // random effects
  matrix[nriver, nq] gamma_river;
  matrix[nreach, nq] gamma_reach;
  matrix[nsite, nq] gamma_site;
  matrix[nyear, nq] gamma_year;
  
  // priors for correlation matrices (Cholesky decomposition)
  cholesky_factor_corr[nsp] L_sp;
  cholesky_factor_corr[nclass] L_class;
  
  // random variates used to define matrix/multinormal  
  matrix[nsp, nclass] eps[N];

  // priors for variances (added to correlation matrices)
  vector<lower=0>[nsp] z_sigma_sp;
  vector<lower=0>[nclass] z_sigma_class;

  // variances for random terms
  row_vector<lower=0>[nq] sigma_river;
  row_vector<lower=0>[nq] sigma_reach;
  row_vector<lower=0>[nq] sigma_site;
  row_vector<lower=0>[nq] sigma_year;

  // dispersion parameter
  real<lower=0> phi;

  // define some matrices for unobserved responses (for AR term)
  matrix[nmissing, nq] log_y_missing;
  matrix[nmissing, nq] log_y_shift_missing;
  
}

transformed parameters {
  row_vector[nq] alpha;
  matrix[K, nq] beta;
  matrix[N, nq] mu;
  matrix[N, nq] mu_shift;
  matrix[N, nq] beta_term;
  row_vector[nq] rho;
  row_vector[nq] tau;
  row_vector[nq] ar_term;
  vector<lower=0>[nsp] sigma_sp;
  vector<lower=0>[nclass] sigma_class;
  matrix[nsp, nsp] W_sp;
  matrix[nclass, nclass] W_class;

  // define covariance matrix from correlations and standard deviations
  sigma_sp = sigma_scale * z_sigma_sp;
  sigma_class = sigma_scale * z_sigma_class;
  W_sp = diag_pre_multiply(sigma_sp, L_sp);
  W_class = diag_pre_multiply(sigma_class, L_class);
  
  // rescale z-transformed covariate effects
  alpha = main_scale * z_alpha;
  beta = main_scale * z_beta;

  // rescale z-transformed AR parameter
  rho = ar_scale * z_rho;
  tau = ar_scale * z_tau;
  
  // define mean value mu for each observation
  beta_term = X * beta;
  for (i in 1:N) {
    ar_term = rep_row_vector(0.0, nq);
    if (visited_prev[i])
      ar_term = rho .* to_row_vector(log_y[prev_idx[i]]) + tau .* to_row_vector(log_y_shift[prev_idx[i]]);
    else
      ar_term = rho .* to_row_vector(log_y_missing[prev_idx[i]]) + tau .* to_row_vector(log_y_shift_missing[prev_idx[i]]);
    mu[i] = alpha + beta_term[i, ] + ar_term +
      sigma_scale2 * sigma_river .* gamma_river[river[i]] +
      sigma_scale2 * sigma_reach .* gamma_reach[reach[i]] + 
      sigma_scale2 * sigma_site .* gamma_site[site[i]] + 
      sigma_scale2 * sigma_year .* gamma_year[year[i]] +
      rep_row_vector(log_effort[i], nq) + 
      to_row_vector((W_sp * (eps[i] * W_class)));
  }

}

model {
  
  // priors for fixed part of linear predictor
  z_alpha ~ std_normal();
  to_vector(z_beta) ~ std_normal();
  
  // prior for AR term
  z_rho ~ std_normal();
  z_tau ~ std_normal();
  
  // random effects
  to_vector(gamma_river) ~ std_normal();
  to_vector(gamma_reach) ~ std_normal();
  to_vector(gamma_site) ~ std_normal();
  to_vector(gamma_year) ~ std_normal();
  
  // and their variances
  sigma_river ~ std_normal();
  sigma_reach ~ std_normal();
  sigma_site ~ std_normal();
  sigma_year ~ std_normal();

  // priors to species- and class-level variances
  z_sigma_sp ~ std_normal();
  z_sigma_class ~ std_normal();
  
  // priors for correlation matrices
  L_sp ~ lkj_corr_cholesky(2.);
  L_class ~ lkj_corr_cholesky(2.);
  
  // priors for missing observations
  to_vector(log_y_missing) ~ std_normal();
  to_vector(log_y_shift_missing) ~ std_normal();
  
  // dispersion prior
  phi ~ std_normal();
  
  // residuals
  for (i in 1:N)
    to_vector(eps[i]) ~ std_normal();
    
  // main likelihood term
  yflat ~ neg_binomial_2_log(to_vector(mu), phi_scale * phi);

}

generated quantities {
  matrix[nsp, nsp] Sigma_sp;
  matrix[nclass, nclass] Sigma_class;

  // calculate correlation matrices from Cholesky decompositions  
  Sigma_sp = W_sp * W_sp';
  Sigma_class = W_class * W_class';

}
