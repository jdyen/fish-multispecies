// matrix normal model with simplified parameterisation of covariances
data {
  int<lower=1> N;                     // num observations
  int<lower=1> K;                     // num covariates
  int<lower=1> ntarget;               // number of target classes
  int<lower=1> nother;                // number of off-target classes
  int<lower=1> nq;                    // ntarget * nother
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
}

transformed data {
  matrix[N, K] Q_ast;
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  int ystart;
  int yend;
  int yflat[N * nq];
  real log_y[N, nq];
  real log_y_shift[N, nq];

  // thin and scale the QR decomposition
  Q_ast = qr_Q(X)[, 1:K] * sqrt(N - 1);
  R_ast = qr_R(X)[1:K, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
  
  // flatten response
  for (i in 1:N) {
    ystart = (i - 1) * nq + 1;
    yend = (i - 1) * nq + nq;
    yflat[ystart:yend] = y[i];
  }
  
  // define log-transformed version of (y + 1) and a size-shifted version of this
  //   for use in the AR model
  log_y = log(yp1);
  for (i in 1:ntarget) {
    log_y_shift[, (i - 1) * nother + 1] = log_y[, i * nother];
    log_y_shift[, ((i - 1) * nother + 2):(i * nother)] = log_y[, ((i - 1) * nother + 1):(i * nother - 1)];
  }

}

parameters {
  
  // overall mean
  row_vector[nq] z_alpha;
  
  // z-scaled regression coefficients in Q space
  matrix[K, nq] z_theta;

  // AR parameter
  row_vector[nq] z_rho;
  row_vector[nq] z_tau;
  
  // random effects
  matrix[nriver, nq] gamma_river;
  matrix[nreach, nq] gamma_reach;
  matrix[nsite, nq] gamma_site;
  matrix[nyear, nq] gamma_year;
  
  // priors for correlation matrices (Cholesky decomposition)
  cholesky_factor_corr[ntarget] L_target;

  // random variates used to define matrix/multinormal  
  matrix[ntarget, nother] eps[N];

  // priors for variances (added to correlation matrices)
  vector<lower=0>[ntarget] sigma_target;
  vector<lower=0>[nother] sigma_other;

  // variances for random terms
  row_vector<lower=0>[nq] sigma_river;
  row_vector<lower=0>[nq] sigma_reach;
  row_vector<lower=0>[nq] sigma_site;
  row_vector<lower=0>[nq] sigma_year;

  // define some matrices for unobserved responses (for AR term)
  matrix[nmissing, nq] log_y_missing;
  matrix[nmissing, nq] log_y_shift_missing;
  
}

transformed parameters {
  row_vector[nq] alpha;
  matrix[K, nq] theta;
  matrix[N, nq] mu;
  matrix[N, nq] mu_shift;
  matrix[N, nq] theta_term;
  row_vector[nq] rho;
  row_vector[nq] tau;
  row_vector[nq] ar_term;

  // define covariance matrix from correlations and standard deviations
  matrix[ntarget, ntarget] W_target = diag_pre_multiply(sigma_target, L_target);

  // rescale z-transformed covariate effects
  alpha = 2. * z_alpha;
  theta = 2. * z_theta;

  // rescale z-transformed AR parameter
  // check these priors with Stan AR1 models, maybe too narrow?
  rho = 0.5 * z_rho;
  tau = 0.5 * z_tau;
  
  // define mean value mu for each observation
  theta_term = Q_ast * theta;
  for (i in 1:N) {
    ar_term = rep_row_vector(0.0, nq);
    if (visited_prev[i])
      ar_term = rho .* to_row_vector(log_y[prev_idx[i]]) + tau .* to_row_vector(log_y_shift[prev_idx[i]]);
    else
      ar_term = rho .* to_row_vector(log_y_missing[prev_idx[i]]) + tau .* to_row_vector(log_y_shift_missing[prev_idx[i]]);
    mu[i] = alpha + theta_term[i, ] + ar_term +
      sigma_river .* gamma_river[river[i]] +
      sigma_reach .* gamma_reach[reach[i]] + 
      sigma_site .* gamma_site[site[i]] + 
      sigma_year .* gamma_year[year[i]] +
      rep_row_vector(log_effort[i], nq) + 
      to_row_vector((W_target * diag_post_multiply(eps[i], sigma_other)));
  }

}

model {
  
  // priors for fixed part of linear predictor
  z_alpha ~ std_normal();
  to_vector(z_theta) ~ std_normal();
  
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
  sigma_target ~ std_normal();
  sigma_other ~ std_normal();
  
  // priors for correlation matrices
  L_target ~ lkj_corr_cholesky(2.);

  // priors for missing observations
  to_vector(log_y_missing) ~ std_normal();
  to_vector(log_y_shift_missing) ~ std_normal();
  
  // increment log prob
  for (i in 1:N)
    to_vector(eps[i]) ~ std_normal();
  yflat ~ poisson_log(to_vector(mu));

}

generated quantities {
  matrix[ntarget, ntarget] Sigma_target;
  matrix[K, nq] beta;

  // calculate correlation matrices from Cholesky decompositions  
  Sigma_target = W_target * W_target';

  // coefficients converted back to X space
  beta = R_ast_inverse * theta;

}
