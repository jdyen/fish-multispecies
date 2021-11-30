// matrix response with AR1 abundances and unstructured MVN
//   across both dimensions
data {
  int<lower=1> N;               // num observations
  int<lower=1> K;               // num covariates
  int<lower=1> nn;              // number of rows in response
  int<lower=1> nm;              // number of columns in response
  int<lower=1> nq;              // nn * nm
  int<lower=2> ntime;           // number of observed time steps
  int<lower=1> ar_order;        // order of AR process
  int<lower=1> nriver;          // number of rivers
  int<lower=1> nreach;          // number of reaches
  matrix[N * ntime, K] X;       // design matrix of covariates (fixed)
  int<lower=1> river[N];        // river ID for each observation
  int<lower=1> reach[N];        // reach ID for each observation
  int y[N, ntime, nq]; // response
  real<lower=0> eta;
}

transformed data {
  matrix[N * ntime, K] Q_ast;
  matrix[N, K] Q_ast_time[ntime];
  matrix[K, K] R_ast;
  matrix[K, K] R_ast_inverse;
  real<lower = 0> eta2 = eta + (nq - 2) / 2.0;
  vector<lower = 0>[nq-1] shape1;
  vector<lower = 0>[nq-1] shape2;

  shape1[1] = eta2;
  shape2[1] = eta2;
  for(k in 2:(nq-1)) {
    eta2 -= 0.5;
    shape1[k] = k / 2.0;
    shape2[k] = eta2;
  }
  
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
  matrix[ntime, nq] log_lambda[N];
  
  // overall mean
  matrix[nn, nm] z_alpha;
  
  // autoregressive parameters
  matrix[nn, nm] rho[ar_order];
  matrix[nn, nm] tau[ar_order];
  
  // priors for correlation matrices (Cholesky decomposition)
  row_vector[choose(nq, 2) - 1]  l;         // do NOT init with 0 for all elements
  vector<lower = 0,upper = 1>[nq-1] R2; // first element is not really a R^2 but is on (0,1)  

  // priors for variances (added to correlation matrices)
  real<lower=0> sigma_main;
  vector<lower=0>[nq] sigma;
  vector[nq] eps[N, ntime];

  // z-scaled regression coefficients in Q space
  matrix[nn, nm] z_theta[K];
  
  // z-scaled gamma terms with MVN prior on species term
  matrix[nriver, nn] z_gamma_river_n;
  matrix[nriver, nm] z_gamma_river_m;
  matrix[nreach, nn] z_gamma_reach_n;
  matrix[nreach, nm] z_gamma_reach_m;
  matrix[ntime, nn] z_gamma_year_n;
  matrix[ntime, nm] z_gamma_year_m;
  
  // and their standard deviations
  vector<lower=0>[nn] sigma_gamma_river_n;
  vector<lower=0>[nm] sigma_gamma_river_m;
  vector<lower=0>[nn] sigma_gamma_reach_n;
  vector<lower=0>[nm] sigma_gamma_reach_m;
  vector<lower=0>[nn] sigma_gamma_year_n;
  vector<lower=0>[nm] sigma_gamma_year_m;
  
}

transformed parameters {
  matrix[nn, nm] alpha;
  matrix[nn, nm] theta[K];
  matrix[nn, nm] gamma_river[nriver];
  matrix[nn, nm] gamma_reach[nreach];
  matrix[nn, nm] gamma_year[ntime];
  matrix[nn, nm] mu[N, ntime];
  matrix[ntime, nq] mu_flat[N];
  matrix[nn, nm] theta_term;
  matrix[nn, nm] ar_term;

  // calculate full covariance
  // matrix[nq,nq] W;
  matrix[nq,nq] L = rep_matrix(0, nq, nq);
  {
    int start = 1;
    int end = 2;

    L[1,1] = 1.0;
    L[2,1] = 2.0 * R2[1] - 1.0;
    L[2,2] = sqrt(1.0 - square(L[2,1]));
    for(k in 2:(nq-1)) {
      int kp1 = k + 1;
      row_vector[k] l_row = segment(l, start, k);
      real scale = sqrt(R2[k] / dot_self(l_row));
      L[kp1, 1:k] = l_row[1:k] * scale;
      L[kp1, kp1] = sqrt(1.0 - R2[k]);
      start = end + 1;
      end = start + k - 1;

    }
  }

  // rescale z-transformed covariate effects
  alpha = 5.0 * z_alpha;
  for (k in 1:K) {
    theta[k] = 5.0 * z_theta[k];
  }
  
  // rescale gamma (random) terms
  for (i in 1:nriver) {
    gamma_river[i] = rep_matrix(sigma_gamma_river_n .* z_gamma_river_n[i]', nm) +
      rep_matrix(sigma_gamma_river_m .* z_gamma_river_m[i]', nn)';
  }
  for (i in 1:nreach) {
    gamma_reach[i] = rep_matrix(sigma_gamma_reach_n .* z_gamma_reach_n[i]', nm) +
  rep_matrix(sigma_gamma_reach_m .* z_gamma_reach_m[i]', nn)';
  }
  for (i in 1:ntime) {
    gamma_year[i] = rep_matrix(sigma_gamma_year_n .* z_gamma_year_n[i]', nm) +
      rep_matrix(sigma_gamma_year_m .* z_gamma_year_m[i]', nn)';
  }

  // define mean value mu for each observation
  for (i in 1:N) {

    // deal with terms too early to specify an AR process
    for (t in 1:ar_order) {
      
      // need to pull out each element of theta and add
      //   them manually because it is a 3D array
      theta_term = rep_matrix(0.0, nn, nm);
      for (k in 1:K) {
        theta_term += Q_ast_time[t, i, k] * theta[k];
      }
      
      // add mean, theta, and random river effects
      mu[i, t] = alpha + theta_term + gamma_river[river[i]] + gamma_reach[reach[i]] + gamma_year[t];
      
      // define flattened version of mu
      mu_flat[i][t] = (to_vector(mu[i, t]) + sigma .* (L * eps[i, t]))';

    }
  
    // then add AR process for later time steps
    for (t in (ar_order + 1):ntime) {
      
      // as above
      theta_term = rep_matrix(0.0, nn, nm);
      for (k in 1:K) {
        theta_term += Q_ast_time[t, i, k] * theta[k];
      }
      
      // define AR term that depends on same and preceding size class
      ar_term = rep_matrix(0.0, nn, nm);
      for (w in 1:ar_order) {
        ar_term += rho[w] .* mu[i, t-w];
        ar_term[1:nn, 1] += tau[w][1:nn, 1] .* mu[i, t-w][1:nn, nm];
        ar_term[1:nn, 2:nm] += tau[w][1:nn, 2:nm] .* mu[i, t-w][1:nn, 1:(nm-1)];
      }
      
      // as above
      mu[i, t] = ar_term + alpha + theta_term + gamma_river[river[i]] + gamma_reach[reach[i]] + gamma_year[t];

      // define flattened version of mu
      mu_flat[i][t] = (to_vector(mu[i, t]) + sigma .* (L * eps[i, t]))';

    }
    
    
  }
  
}

model {
  
  // priors for linear predictor
  sigma_gamma_river_n ~ std_normal();
  sigma_gamma_river_m ~ std_normal();
  to_vector(z_gamma_river_n) ~ std_normal();
  to_vector(z_gamma_river_m) ~ std_normal();
  sigma_gamma_reach_n ~ std_normal();
  sigma_gamma_reach_m ~ std_normal();
  to_vector(z_gamma_reach_n) ~ std_normal();
  to_vector(z_gamma_reach_m) ~ std_normal();
  sigma_gamma_year_n ~ std_normal();
  sigma_gamma_year_m ~ std_normal();
  to_vector(z_gamma_year_n) ~ std_normal();
  to_vector(z_gamma_year_m) ~ std_normal();
  
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
  sigma ~ std_normal();
  sigma_main ~ std_normal();

  // priors that imply L ~ lkj_corr_cholesky(eta)
  // L ~ lkj_corr_cholesky(2.0);
  l ~ std_normal();
  R2 ~ beta(shape1, shape2);

  // increment log prob
  for (i in 1:N) {
    to_vector(log_lambda[i]) ~ normal(to_vector(mu_flat[i]), sigma_main);
    for (t in 1:ntime) {
      eps[i, t] ~ std_normal();
      target += poisson_log_lpmf(y[i, t] | log_lambda[i, t]);
    }
  }
  
}

generated quantities {
  matrix[nq, nq] Omega;
  matrix[nn, nm] beta[K];
  
  // calculate correlation matrices from Cholesky decompositions  
  Omega = L * L';
  
  // coefficients converted back to X space
  for (k in 1:K) {
    beta[k] = rep_matrix(0.0, nn, nm);
    for (w in 1:K) {
      beta[k] += R_ast_inverse[k, w] * theta[w];
    }
  }
  
}
