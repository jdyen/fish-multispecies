functions {
  // function to calculate log pdf of matrix normal
  //   from Cholesky decomposition of covariance
  //   matrices
  //     - y: nsp x nclass matrix (observations)
  //     - mu: nsp x nclass matrix (linear predictor)
  //     - U_chol: Cholesky decomposition of covariance for
  //          nsp dimension
  //     - V_chol: Cholesky decomposition of covariance for
  //          nclass dimension
  real matrix_normal_cholesky_lpdf(matrix y, matrix mu, matrix U_chol, matrix V_chol) {
  
    // work out dimensions
    int n = rows(U_chol);
    int m = rows(V_chol);

    // calculate full covariance matrices (inverses used below)
    matrix[n, n] U = multiply_lower_tri_self_transpose(U_chol);
    matrix[m, m] V = multiply_lower_tri_self_transpose(V_chol);

    // calculate residual
    matrix[n, m] y_resid = y - mu;
  
    // calculate numerator of log pdf (inverse of covariances multiplied
    //   by residuals)
    matrix[n, m] P = U \ y_resid; // inverse(U) * y_resid
    matrix[m, n] Q = V \ y_resid'; // inverse(V) * y_resid'
    // possible alternative but seems to be slower:
    //matrix[n, m] P = mdivide_left_tri_low(U_chol, mdivide_right_tri_low(y_resid', U_chol)'); // inverse(U) * y_resid 

    real numerator = trace(Q * P);
  
    // calculate normalisation term (denominator of log pdf)
    real log2pi = 1.8378770664093453;
    real denom = n * m * log2pi + 2 * n * sum(diagonal(U_chol)) + 2 * m * sum(diagonal(V_chol));
  
    // and return calculated log pdf
    return -0.5 * (numerator + denom);
  
  }
  
}

// matrix normal model with AR1 abundances
data {
  int<lower=1> N;               // num observations
  int<lower=1> K;               // num covariates
  int<lower=1> nsp;             // number of species
  int<lower=1> nclass;          // number of size classes
  int<lower=2> ntime;           // number of observed time steps
  int<lower=1> ar_order;        // order of AR process
  int<lower=1> nriver;          // number of rivers
  int<lower=1> nreach;          // number of reaches
  matrix[N * ntime, K] X;       // design matrix of covariates (fixed)
  int<lower=1> river[N];        // river ID for each observation
  int<lower=1> reach[N];        // reach ID for each observation
  int y[N, ntime, nsp * nclass]; // response
}

transformed data {
  matrix[N * ntime, K] Q_ast;
  matrix[N, K] Q_ast_time[ntime];
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;

  // thin and scale the QR decomposition
  Q_ast = qr_Q(X)[, 1:K] * sqrt(N - 1);
  R_ast = qr_R(X)[1:K, ] / sqrt(N - 1);
  R_ast_inverse = inverse(R_ast);
  
  // reformat predictor matrix to have time slices
  for (t in 1:ntime) {
    Q_ast_time[t] = X[((t - 1) * N + 1):(t * N)];
  }
  
}

parameters {
  
  // linear predictor with gaussian noise
  matrix[nsp, nclass] log_lambda[N, ntime];

  // overall mean
  matrix[nsp, nclass] z_alpha;
  
  // autoregressive parameters
  matrix[nsp, nclass] rho[ar_order];
  matrix[nsp, nclass] tau[ar_order];
  
  // priors for correlation matrices (Cholesky decomposition)
  cholesky_factor_corr[nsp] L_sp;
  cholesky_factor_corr[nclass] L_class;
  
  // priors for variances (added to correlation matrices)
  vector<lower=0>[nsp] sigma_sp;
  vector<lower=0>[nclass] sigma_class;
  
  // z-scaled regression coefficients in Q space
  matrix[nsp, nclass] z_theta[K];
  
  // z-scaled gamma terms with MVN prior on species term
  matrix[nriver, nsp] z_gamma_river_sp;
  matrix[nriver, nclass] z_gamma_river_class;
  matrix[nreach, nsp] z_gamma_reach_sp;
  matrix[nreach, nclass] z_gamma_reach_class;
  matrix[ntime, nsp] z_gamma_year_sp;
  matrix[ntime, nclass] z_gamma_year_class;

  // and their standard deviations
  vector<lower=0>[nsp] sigma_gamma_river_sp;
  vector<lower=0>[nclass] sigma_gamma_river_class;
  vector<lower=0>[nsp] sigma_gamma_reach_sp;
  vector<lower=0>[nclass] sigma_gamma_reach_class;
  vector<lower=0>[nsp] sigma_gamma_year_sp;
  vector<lower=0>[nclass] sigma_gamma_year_class;

}

transformed parameters {
  matrix[nsp, nclass] alpha;
  matrix[nsp, nclass] theta[K];
  matrix[nsp, nclass] gamma_river[nriver];
  matrix[nsp, nclass] gamma_reach[nreach];
  matrix[nsp, nclass] gamma_year[ntime];
  matrix[nsp, nclass] mu[N, ntime];
  matrix[nsp, nclass] theta_term;
  matrix[nsp, nclass] ar_term;
  
  // define covariance matrix from correlations and standard deviations
  matrix[nsp, nsp] W_sp = diag_pre_multiply(sigma_sp, L_sp);
  matrix[nclass, nclass] W_class = diag_pre_multiply(sigma_class, L_class);

  // rescale z-transformed covariate effects
  alpha = 5.0 * z_alpha;
  for (k in 1:K) {
    theta[k] = 5.0 * z_theta[k];
  }
  
  // rescale gamma (random) terms
  for (i in 1:nriver) {
    gamma_river[i] = rep_matrix(sigma_gamma_river_sp .* z_gamma_river_sp[i]', nclass) +
      rep_matrix(sigma_gamma_river_class .* z_gamma_river_class[i]', nsp)';
  }
  for (i in 1:nreach) {
    gamma_reach[i] = rep_matrix(sigma_gamma_reach_sp .* z_gamma_reach_sp[i]', nclass) +
      rep_matrix(sigma_gamma_reach_class .* z_gamma_reach_class[i]', nsp)';
  }
  for (i in 1:ntime) {
    gamma_year[i] = rep_matrix(sigma_gamma_year_sp .* z_gamma_year_sp[i]', nclass) +
      rep_matrix(sigma_gamma_year_class .* z_gamma_year_class[i]', nsp)';
  }

  // define mean value mu for each observation
  for (i in 1:N) {

    // deal with terms too early to specify an AR process
    for (t in 1:ar_order) {
      
      // need to pull out each element of theta and add
      //   them manually because it is a 3D array
      theta_term = rep_matrix(0.0, nsp, nclass);
      for (k in 1:K) {
        theta_term += Q_ast_time[t, i, k] * theta[k];
      }
      
      // add mean, theta, and random river effects
      mu[i, t] = alpha + theta_term + gamma_river[river[i]] + gamma_reach[reach[i]] + gamma_year[t];
      
    }
    
    // then add AR process for later time steps
    for (t in (ar_order + 1):ntime) {
      
      // as above
      theta_term = rep_matrix(0.0, nsp, nclass);
      for (k in 1:K) {
        theta_term += Q_ast_time[t, i, k] * theta[k];
      }
      
      // define AR term that depends on same and preceding size class
      ar_term = rep_matrix(0.0, nsp, nclass);
      for (w in 1:ar_order) {
        ar_term += rho[w] .* mu[i, t-w];
        ar_term[1:nsp, 1] += tau[w][1:nsp, 1] .* mu[i, t-w][1:nsp, nclass];
        ar_term[1:nsp, 2:nclass] += tau[w][1:nsp, 2:nclass] .* mu[i, t-w][1:nsp, 1:(nclass-1)];
      }
      
      // as above
      mu[i, t] = ar_term + alpha + theta_term + gamma_river[river[i]] + gamma_reach[reach[i]] + gamma_year[t];
      
    }
    
  }
  
}

model {
  
  // priors for linear predictor
  sigma_gamma_river_sp ~ std_normal();
  sigma_gamma_river_class ~ std_normal();
  to_vector(z_gamma_river_sp) ~ std_normal();
  to_vector(z_gamma_river_class) ~ std_normal();
  sigma_gamma_reach_sp ~ std_normal();
  sigma_gamma_reach_class ~ std_normal();
  to_vector(z_gamma_reach_sp) ~ std_normal();
  to_vector(z_gamma_reach_class) ~ std_normal();
  sigma_gamma_year_sp ~ std_normal();
  sigma_gamma_year_class ~ std_normal();
  to_vector(z_gamma_year_sp) ~ std_normal();
  to_vector(z_gamma_year_class) ~ std_normal();
  
  // priors for fixed part of linear predictor
  for (k in 1:K) {
    to_vector(z_theta[k]) ~ std_normal();
  }
  to_vector(z_alpha) ~ std_normal();
  
  // prior for AR1 term
  for (w in 1:ar_order) {
    to_vector(rho[w]) ~ std_normal();
    to_vector(tau[w]) ~ std_normal();
  }
  
  // priors to species- and class-level variances
  sigma_sp ~ std_normal();
  sigma_class ~ std_normal();
  
  // priors for correlation matrices
  L_sp ~ lkj_corr_cholesky(1);
  L_class ~ lkj_corr_cholesky(1);
  
  // increment log prob
  for (i in 1:N) {
    for (t in 1:ntime) {
      log_lambda[i, t] ~ matrix_normal_cholesky(mu[i, t], W_sp, W_class);
      target += poisson_log_lpmf(y[i, t] | to_vector(log_lambda[i, t]));
    }
  }
  
}

generated quantities {
  matrix[nsp, nsp] Omega_sp;
  matrix[nclass, nclass] Omega_class;
  matrix[nsp, nclass] beta[K];

  // calculate correlation matrices from Cholesky decompositions  
  Omega_sp = L_sp * L_sp';
  Omega_class = L_class * L_class';

  // coefficients converted back to X space
  for (k in 1:K) {
    beta[k] = rep_matrix(0.0, nsp, nclass);
    for (w in 1:K) {
      beta[k] += R_ast_inverse[k, w] * theta[w];
    }
  }
  
}
