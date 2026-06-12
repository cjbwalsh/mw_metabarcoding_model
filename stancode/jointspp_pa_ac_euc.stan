data {
  int<lower=1> n_site;                   // Number of sites
  int<lower=1> n_taxa;                   // Number of species
  int<lower=1> n_pred;                   // Number of fixed predictors
  int<lower=1> n_obs;                    // Total number of samples 
  
  array[n_obs] int<lower=1, upper=n_site> site; // Site index for each sample

  matrix[n_obs, n_pred] u;               // Matrix of scaled predictors for each observation
  matrix[n_taxa, n_taxa] phylo_cor;      // Phylogenetic correlation matrix
  matrix[n_site, n_site] euc_dist_mat;    // Euclidean distance between sites standardized so maax = 1
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

  // Uncorrelated error (Overdispersion) - removed mu_epsilon to eliminate identical collinearity between the species intercept and the overdispersion mean
  // vector[n_taxa] mu_epsilon;             
  matrix[n_obs, n_taxa] epsilon_raw;     
  real<lower=0> sigma_epsilon;           // Fixed type & bounds for non-centering
  
  // Between-predictor correlation matrix prior
  cholesky_factor_corr[n_pred] L_Omega;           

  // Site random effects - removed mu_site to eliminate a second source of identical collinearity with a_taxon
  matrix[n_site, n_taxa] a_site_raw;           
  // vector[n_taxa] mu_site;                     
  vector<lower=0>[n_taxa] sigma_site;    // Added variance tracking for site effects

  // Taxon intercept
  vector[n_taxa] a_taxon_raw;                 
  real mu_taxon;                              
  real<lower=0> sigma_taxon;             // Added variance tracking for taxon intercepts

  // Spatial autocorrelation:
  matrix<lower=0>[n_site,n_taxa] delta_raw;    // rawspatial autocorrelation adjustment
  real<lower=0> sigma_space;                   // scale parameter
  real<lower=0> rho_space;                     // length scale parameter, controlling decay with distance
}

transformed parameters {
  matrix[n_obs, n_taxa] mu;                   

  matrix[n_pred, n_taxa] beta;                
  matrix[n_site, n_taxa] a_site;              
  vector[n_taxa] a_taxon;
  matrix[n_obs, n_taxa] epsilon;              

  // Spatial autocorrelation parameters
  real d_adj = 1e-5;          // Jitter added to distance matrix diagonal for stability
  matrix[n_site, n_site] Sigma_space;
  matrix[n_site, n_site] L_Sigma_space;
  matrix[n_site, n_taxa] delta;             

  // 1. Non-Centered Phylogenetic Regression Slopes
  {
  matrix[n_pred, n_taxa] beta_phylo_raw = beta_raw * L_phylo'; 
  matrix[n_pred, n_pred] L_Sigma_beta = diag_pre_multiply(scale_beta, L_Omega);
  
  for (k in 1:n_pred) {
    for (s in 1:n_taxa) {
      beta[k, s] = mu_beta[k] + dot_product(L_Sigma_beta[k, ], beta_phylo_raw[, s]);
    }
  }
}

  // Spatial autocorrelation effect (squared exponential decay)
  Sigma_space = square(sigma_space) * exp(-0.5 * square(euc_dist_mat / rho_space));
  // Make the diagonal elements are strictly positive
  for (i in 1:n_site) {
    Sigma_space[i, i] = Sigma_space[i, i] + d_adj;
  }
  // Compute the lower-triangular Cholesky factor
  L_Sigma_space = cholesky_decompose(Sigma_space);
  delta = L_Sigma_space * delta_raw;
  
  // Non-Centered Transformations for Hierarchical Variances
  a_taxon = mu_taxon + a_taxon_raw * sigma_taxon;

  // Vectorized site and overdispersion transformations (Removes loop overhead)
  a_site = a_site_raw .* rep_matrix(sigma_site', n_site);
  epsilon = epsilon_raw * sigma_epsilon; // removing mu_epsilon centers at 0.

  for (j in 1:n_taxa) {
    // Vectorized row-wise assignment significantly helps energy/BFMI geometry
    mu[, j] = a_taxon[j] + a_site[site, j] + (u * col(beta, j)) + delta[site, j] + epsilon[, j];
  }
}

model {
  // --- Priors for the betas ---
  to_vector(beta_raw) ~ std_normal();
  mu_beta ~ normal(0, 2);
  scale_beta ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(2.0); 

  // --- Spatial Priors ---
  to_vector(delta_raw) ~ std_normal();
  sigma_space ~ normal(0, 0.5);
  rho_space ~ inv_gamma(2, 0.5); // Prior concentrated on the 0-1 scaled distance range

  // --- Genuine Random Effects Priors ---
  // mu_site ~ normal(0, 5);    //removed above             
  to_vector(a_site_raw) ~ std_normal();   
  sigma_site ~ normal(0, 1.5);

  mu_taxon ~ normal(0, 3);                
  a_taxon_raw ~ std_normal();             
  sigma_taxon ~ normal(0, 1.5);
  
  // --- Overdispersion Priors (Fixed Scalar Syntax) ---
  // mu_epsilon ~ normal(0, 1);  //removed above             
  to_vector(epsilon_raw) ~ std_normal();  
  sigma_epsilon ~ normal(0, 0.5);  

  // --- Likelihood Estimation ---
  for (j in 1:n_taxa) {
    y[, j] ~ bernoulli_logit(mu[, j]); 
  }
}

generated quantities {
  array[n_obs, n_taxa] int<lower=0, upper=1> y_rep; // Simulated replica data
  vector[n_obs * n_taxa] log_lik;                  // Flattened observation-by-taxon log-likelihood
  vector[n_taxa] tjurs_r2;                         // Explanatory power per taxon
 
  {
    int idx = 1; // Counter to flatten the log_lik vector
    
    for (j in 1:n_taxa) {
      real sum_prob_pres = 0.0;
      real sum_prob_abs = 0.0;
      real n_pres = 0.0;
      real n_abs = 0.0;
      
      for (i in 1:n_obs) {
        real prob = inv_logit(mu[i, j]);
        
        // 1. Generate posterior predictive data
        y_rep[i, j] = bernoulli_rng(prob); 
        
        // 2. Calculate point-level log-likelihood (Observation x Taxon)
        log_lik[idx] = bernoulli_logit_lpmf(y[i, j] | mu[i, j]);
        idx += 1;
        
        // 3. Track values for Tjur's R2 calculation
        if (y[i, j] == 1) {
          sum_prob_pres += prob;
          n_pres += 1.0;
        } else {
          sum_prob_abs += prob;
          n_abs += 1.0;
        }
      }
      
      // Calculate final Tjur's R2 for the taxon
      if (n_pres > 0 && n_abs > 0) {
        tjurs_r2[j] = (sum_prob_pres / n_pres) - (sum_prob_abs / n_abs);
      } else {
        tjurs_r2[j] = 0.0; 
      }
    }
  }
}
