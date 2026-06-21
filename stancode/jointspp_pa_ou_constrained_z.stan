data {
  int<lower=1> n_obs;                     // Total number of samples 
  int<lower=1> n_site;                    // Number of sites
  int<lower=1> n_taxa;                    // Number of species
  int<lower=1> n_site_pred;               // Number of site-level predictors
  int<lower=1> n_obs_pred;                // Number of observation-level predictors

  matrix[n_site, n_site_pred] u_site;     // Site-level predictors
  matrix[n_obs, n_obs_pred] u_obs;        // Sample-level predictors
  int<lower=1> n_latent;                  // Number of latent factors (set at 2)
  matrix[n_site, n_site] M_site_res;      // Residual projection matrix for latent-factor estimation
  matrix[n_taxa, n_taxa] phylo_cor;       // Phylogenetic correlation matrix

  array[n_obs] int<lower=1, upper=n_site> site; // Site index for each sample

  array[n_obs, n_taxa] int<lower=0, upper=1> y; // Presence/absence matrix
}

parameters {
  // Slopes for static site-specific predictors
  vector[n_site_pred] mu_beta_site;                
  matrix[n_site_pred, n_taxa] beta_site_raw;       
  vector<lower=0>[n_site_pred] scale_beta_site;    

  // Slopes for sample-specific predictors (season, riff)
  vector[n_obs_pred] mu_beta_obs;                
  matrix[n_obs_pred, n_taxa] beta_obs_raw;       
  vector<lower=0>[n_obs_pred] scale_beta_obs;    
  
  // strength of phylogenetic influence on betas (exponential decay parameter for Ornstein-Uhlenbeck divergence)
  real<lower=0, upper=20> alpha_phylo; 

  // Latent variables
  matrix[n_site, n_latent] z_raw;
  real<lower=0> sigma_z;
  real<lower=0> sigma_lambda;
  
  vector<lower=0>[n_latent] lambda_diag; // Diagonal elements forced positive
  real lambda_lower;                     // Row 2, Column 1 elements
  matrix[n_taxa - 2, n_latent] lambda_remaining; // Rows 3 to n_taxa are unconstrained

  // Between-predictor correlation matrices prior
  cholesky_factor_corr[n_site_pred] L_Omega_site;
  cholesky_factor_corr[n_obs_pred]  L_Omega_obs;

  // Taxon intercept
  vector[n_taxa] a_taxon_raw;                 
  real mu_taxon;                              
  real<lower=0> sigma_taxon;             // Added variance tracking for taxon intercepts

}

transformed parameters {
  // 1. Non-Centered Phylogenetic Regression Slopes
  matrix[n_site_pred, n_taxa] beta_site;
  matrix[n_obs_pred, n_taxa] beta_obs;
  vector[n_taxa] a_taxon;
  
  matrix[n_taxa, n_latent] lambda;
  matrix[n_obs, n_latent] z_expanded;
   // 1. NATIVELY STITCH THE IDENTIFIABLE LOWER-TRIANGULAR MATRIX
    // Row 1: Species 1 loads ONLY on Axis 1 (Axis 2 is forced to exactly 0.0)
    lambda[1, 1] = lambda_diag[1] * sigma_lambda;
    lambda[1, 2] = 0.0; 
    
    // Row 2: Species 2 loads on Axis 1 (free) and Axis 2 (forced positive)
    lambda[2, 1] = lambda_lower * sigma_lambda;
    lambda[2, 2] = lambda_diag[2] * sigma_lambda;
    
    // Rows 3 through 875: Drop the rest of the community into place unconstrained
    lambda[3:n_taxa, ] = lambda_remaining * sigma_lambda;
    
    {
      matrix[n_site, n_latent] z_clean = M_site_res * z_raw;
      z_expanded = z_clean[site] * sigma_z;
    }

  // Non-Centered Transformations for Hierarchical Variances
  a_taxon = mu_taxon + a_taxon_raw * sigma_taxon;

    {
    matrix[n_taxa, n_taxa] sigma_ou;
    matrix[n_taxa, n_taxa] L_phylo_mixed;
    sigma_ou= exp(-2.0 * alpha_phylo * (1.0 - phylo_cor));
    for (i in 1:n_taxa) {
        sigma_ou[i, i] = 1.0;
    }
   L_phylo_mixed = cholesky_decompose(sigma_ou);

    // Match site-level slopes to the mixed tree
    beta_site = rep_matrix(mu_beta_site, n_taxa) + 
                (diag_pre_multiply(scale_beta_site, L_Omega_site) * beta_site_raw * L_phylo_mixed);
    // Match sample-level slopes to the exact same mixed tree
    beta_obs = rep_matrix(mu_beta_obs, n_taxa) + 
               (diag_pre_multiply(scale_beta_obs, L_Omega_obs) * beta_obs_raw * L_phylo_mixed);
  }
}

model {
  // --- Priors for the betas ---
  to_vector(beta_site_raw) ~ std_normal();
  mu_beta_site ~ normal(0, 2);
  scale_beta_site ~ normal(0, 2);
  L_Omega_site ~ lkj_corr_cholesky(2.0); 

  to_vector(beta_obs_raw) ~ std_normal();
  mu_beta_obs ~ normal(0, 2);
  scale_beta_obs ~ normal(0, 2);
  L_Omega_obs ~ lkj_corr_cholesky(2.0); 

  mu_taxon ~ normal(0, 3);                
  a_taxon_raw ~ std_normal();             
  sigma_taxon ~ normal(0, 1.5);
  
  // --- Latent Factor Priors ---
  to_vector(z_raw) ~ std_normal();      // Restricts site coordinates to standard normal space
  // Regularising priors on the scales to keep the logit linear predictors stable
  // Tight half-normal priors act as an L2 penalty to keep the geometry clean
  sigma_z ~ normal(0, 0.5);      
  sigma_lambda ~ normal(0, 0.5); 
  
  // Standard normal priors apply directly across the split components
  lambda_diag ~ std_normal();
  lambda_lower ~ std_normal();
  to_vector(lambda_remaining) ~ std_normal();

  alpha_phylo ~ normal(0, 5);

  // by declaring U_site_expanded and mu inside these braces, they are not saved in the model fit.
  {
    matrix[n_obs, n_taxa] mu;      
    matrix[n_obs, n_site_pred] U_site_expanded = u_site[site]; // Optimization block

    // Compile mu column-by-column using your corrected mathematical orientation
    for (j in 1:n_taxa) {
      mu[, j] = a_taxon[j] +                          // Intercept for species j
                (U_site_expanded * col(beta_site, j)) +  // Site predictor effects
                (u_obs * col(beta_obs, j)) +          // Sample predictor effects
                (z_expanded * lambda[j, ]');          // Latent site effects
    }
    // --- Likelihood Estimation ---
   for (j in 1:n_taxa) {
      y[, j] ~ bernoulli_logit(mu[, j]); 
      }
  }
}

generated quantities {
  array[n_obs, n_taxa] int<lower=0, upper=1> y_rep; // Simulated replica data
  vector[n_obs * n_taxa] log_lik;                   // Flattened log-likelihood
  vector[n_taxa] tjurs_r2_total;                    // Explanatory strength of the full model
  vector[n_taxa] tjurs_r2_env;                      // Explanatory strength using only fixed predictors 
  matrix[n_obs, n_taxa] p_total_out;                // For calculation of full-model AUC in R
  matrix[n_obs, n_taxa] p_env_out;                  // For calculation of AUC using only fixed predictors in R
 
  {
    matrix[n_obs, n_taxa] mu;
    matrix[n_obs, n_taxa] mu_env; 
    matrix[n_obs, n_site_pred] U_site_expanded = u_site[site];

   for (j in 1:n_taxa) {
      // mu_env isolates u_site and u_obs effects
      mu_env[, j] = a_taxon[j] + (U_site_expanded * col(beta_site, j)) + (u_obs * col(beta_obs, j));
      // mu adds the unmeasured residual latent site matrix 
      mu[, j]     = mu_env[, j] + (z_expanded * lambda[j, ]'); 
    }

    // Validation metrics
    int idx = 1; 
    for (j in 1:n_taxa) {
      real sum_p_pres_total = 0.0;  real sum_p_abs_total = 0.0;
      real sum_p_pres_env = 0.0;    real sum_p_abs_env = 0.0;
      real n_pres = 0.0;            real n_abs = 0.0;
      
      for (i in 1:n_obs) {
        real p_total = inv_logit(mu[i, j]);
        real p_env   = inv_logit(mu_env[i, j]);
        
        y_rep[i, j] = bernoulli_rng(p_total); 
        log_lik[idx] = bernoulli_logit_lpmf(y[i, j] | mu[i, j]);
        idx += 1;
        
        if (y[i, j] == 1) {
          n_pres += 1.0;
          sum_p_pres_total += p_total;
          sum_p_pres_env   += p_env;
        } else {
          n_abs += 1.0;
          sum_p_abs_total += p_total;
          sum_p_abs_env   += p_env;
        }
      }
      
      if (n_pres > 0 && n_abs > 0) {
        tjurs_r2_total[j] = (sum_p_pres_total / n_pres) - (sum_p_abs_total / n_abs);
        tjurs_r2_env[j]   = (sum_p_pres_env / n_pres) - (sum_p_abs_env / n_abs);
      } else {
        tjurs_r2_total[j] = 0.0; 
        tjurs_r2_env[j]   = 0.0; 
      }
    }
  }
}
