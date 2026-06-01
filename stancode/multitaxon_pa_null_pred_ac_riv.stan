data {
  int<lower=1> n_site;                   // Number of sites
  int<lower=1> n_taxa;                   // Number of species
  int<lower=1> n_pred;                   // Number of fixed predictors
  int<lower=1> n_obs;                    // Total number of samples 
  
  array[n_obs] int<lower=1, upper=n_site> site; // Site index for each sample

  matrix[n_obs, n_pred] u;               // Matrix of scaled predictors for each observation
  matrix[n_taxa, n_taxa] phylo_cor;      // Phylogenetic correlation matrix
  matrix[n_site, n_site] riv_dist_mat;   // River distance matrix (either up or down)
  array[n_obs, n_taxa] int<lower=0, upper=1> y; // Presence/absence matrix
}

transformed data {
  matrix[n_taxa, n_taxa] L_phylo;
  // Pre-decomposing static data here saves massive computational overhead
  L_phylo = cholesky_decompose(phylo_cor);
}

parameters {
  // Fixed effects
  vector[n_pred] mu_beta;                
  matrix[n_pred, n_taxa] beta_raw;       
  vector<lower=0>[n_pred] scale_beta;    

  // Uncorrelated error (Overdispersion)
  vector[n_taxa] mu_epsilon;             
  matrix[n_obs, n_taxa] epsilon_raw;     
  real<lower=0> sigma_epsilon;           // Fixed type & bounds for non-centering
  
  // Between-predictor correlation matrix prior
  cholesky_factor_corr[n_pred] L_Omega;           

  // Site random effects
  matrix[n_site, n_taxa] a_site_raw;           
  vector[n_taxa] mu_site;                     
  vector<lower=0>[n_taxa] sigma_site;    // Added variance tracking for site effects

  // Taxon intercept
  vector[n_taxa] a_taxon_raw;                 
  real mu_taxon;                              
  real<lower=0> sigma_taxon;             // Added variance tracking for taxon intercepts

  // Spatial autocorrelation:
  matrix[n_site, n_taxa] delta_raw;    
  real<lower=0, upper=0.99> rho_space;         
  real<lower=0> sigma_space;                   
}

transformed parameters {
  matrix[n_pred, n_taxa] beta;                
  matrix[n_obs, n_taxa] mu;                   

  matrix[n_site, n_taxa] a_site;              
  matrix[n_obs, n_taxa] epsilon;              
  matrix[n_site, n_taxa] delta;               
  vector[n_taxa] a_taxon;

  // 1. Correctly Non-Centered Phylogenetic Regression Slopes
  {
    matrix[n_pred, n_pred] L_Sigma_beta = diag_pre_multiply(scale_beta, L_Omega);
    for (k in 1:n_pred) {
      for (s in 1:n_taxa) {
        beta[k, s] = mu_beta[k] + dot_product(L_Sigma_beta[k, ], beta_raw[, s]);
      }
    }
    beta = beta * L_phylo'; 
  }

  // 2. Non-centered asymmetric spatial drift (Matrix Solver)
  {
    matrix[n_site, n_site] A;
    A = diag_matrix(rep_vector(1.0, n_site)) - rho_space * riv_dist_mat;
    delta = A \ (delta_raw * sigma_space); 
  }

  // 3. Genuine Non-Centered Transformations for Hierarchical Variances
  a_taxon = mu_taxon + a_taxon_raw * sigma_taxon;

  for (j in 1:n_taxa) {
    // Multiplied raw standard normals by standard deviations for mathematical validity
    a_site[, j] = mu_site[j] + a_site_raw[, j] * sigma_site[j];                    
    epsilon[, j] = mu_epsilon[j] + epsilon_raw[, j] * sigma_epsilon;               
  }

  // 4. Corrected Linear Predictor Mapping Loop
  for (j in 1:n_taxa) {
    for (i in 1:n_obs) { // Fixed: Changed from n_pred to n_obs
      // Solved: Vector multiplication u[i] * col(beta, j) evaluates cleanly to a scalar real.
      // Fixed: Indexed a_site and delta using row AND column coordinates to pull a scalar.
      mu[i, j] = a_taxon[j] + a_site[site[i], j] + (u[i] * col(beta, j)) + 
                 delta[site[i], j] + epsilon[i, j];
    }
  }
}

model {
  // --- Priors for the betas ---
  to_vector(beta_raw) ~ std_normal();
  mu_beta ~ normal(0, 5);
  scale_beta ~ normal(0, 3.5);
  L_Omega ~ lkj_corr_cholesky(2.0); 

  // --- Spatial Priors ---
  to_vector(delta_raw) ~ std_normal();
  sigma_space ~ normal(0, 3.5);
  rho_space ~ beta(2, 2); 

  // --- Genuine Random Effects Priors ---
  mu_site ~ normal(0, 5);                 
  to_vector(a_site_raw) ~ std_normal();   
  sigma_site ~ normal(0, 2.5);

  mu_taxon ~ normal(0, 5);                
  a_taxon_raw ~ std_normal();             
  sigma_taxon ~ normal(0, 2.5);
  
  // --- Overdispersion Priors (Fixed Scalar Syntax) ---
  mu_epsilon ~ normal(0, 1);              
  to_vector(epsilon_raw) ~ std_normal();  
  sigma_epsilon ~ normal(0, 0.5);  

  // --- Likelihood Estimation ---
  for (j in 1:n_taxa) {
    y[, j] ~ bernoulli_logit(mu[, j]); 
  }
}

generated quantities {
  array[n_obs, n_taxa] int<lower=0, upper=1> y_rep; // Simulated replica data
  matrix[n_obs, n_taxa] log_lik;                    // Log-likelihood for WAIC/LOO-CV
  
  // Custom metrics for checking performance
  vector[n_taxa] tjurs_r2;                         // Explanatory power per taxon
  
  {
    // Temporary tracking vectors to calculate Tjur's R2 per taxon
    // Tjur's R2 = (Mean predicted probability when present) - (Mean predicted probability when absent)
    for (j in 1:n_taxa) {
      real sum_prob_pres = 0.0;
      real sum_prob_abs = 0.0;
      real n_pres = 0.0;
      real n_abs = 0.0;
      
      for (i in 1:n_obs) {
        real prob = inv_logit(mu[i, j]);
        
        // 1. Generate posterior predictive data
        y_rep[i, j] = bernoulli_rng(prob);
        
        // 2. Calculate point-wise log-likelihood for LOO-CV
        log_lik[i, j] = bernoulli_logit_lpmf(y[i, j] | mu[i, j]);
        
        // 3. Track values for Tjur's R2 calculation
        if (y[i, j] == 1) {
          sum_prob_pres += prob;
          n_pres += 1.0;
        } else {
          sum_prob_abs += prob;
          n_abs += 1.0;
        }
      }
      
      // Calculate discrimination power safely to avoid dividing by zero if a taxon is always present/absent
      if (n_pres > 0 && n_abs > 0) {
        tjurs_r2[j] = (sum_prob_pres / n_pres) - (sum_prob_abs / n_abs);
      } else {
        tjurs_r2[j] = 0.0; 
      }
    }
  }
}
