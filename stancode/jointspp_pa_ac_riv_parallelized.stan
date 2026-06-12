functions {
  /**
   * Partial sum function for taxon-wise parallel likelihood
   * @param taxa_slice Subset array of taxon indices processed by this thread
   * @param start Starting index of the current slice
   * @param end Ending index of the current slice
   * @param y Master presence/absence matrix [n_obs, n_taxa]
   * @param mu Master linear predictor matrix [n_obs, n_taxa]
   * @return Log-likelihood contribution of this slice
   */
  real partial_taxa_lpmf(array[] int taxa_slice, int start, int end,
                         matrix mu, array[,] int y) {
    real lp = 0;
    for (i in start:end) {
      int j = taxa_slice[i - start + 1]; // Translate to actual taxon index
      lp += bernoulli_logit_lpmf(y[, j] | mu[, j]);
    }
    return lp;
  }
}

data {
  int<lower=1> n_site;                   
  int<lower=1> n_taxa;                   
  int<lower=1> n_pred;                   
  int<lower=1> n_obs;                    
  
  array[n_obs] int<lower=1, upper=n_site> site; 

  matrix[n_obs, n_pred] u;               
  matrix[n_taxa, n_taxa] phylo_cor;      
  matrix[n_site, n_site] riv_dist_mat;   
  matrix[n_site, n_site] is_downstream;   

  array[n_obs, n_taxa] int<lower=0, upper=1> y; 
  
  // Controls memory output in Generated Quantities
  int<lower=0, upper=1> save_y_rep;     // 1 = save massive y_rep matrix, 0 = skip
}

transformed data {
  matrix[n_taxa, n_taxa] L_phylo;
  array[n_taxa] int taxa_indices;
  
  L_phylo = cholesky_decompose(phylo_cor);
  
  // Construct the array of targets for parallel slicing
  for(j in 1:n_taxa) {
    taxa_indices[j] = j;
  }
}

parameters {
  vector[n_pred] mu_beta;                
  matrix[n_pred, n_taxa] beta_raw;       
  vector<lower=0>[n_pred] scale_beta;    

  matrix[n_obs, n_taxa] epsilon_raw;     
  real<lower=0> sigma_epsilon;           
  
  cholesky_factor_corr[n_pred] L_Omega;           

  matrix[n_site, n_taxa] a_site_raw;           
  vector<lower=0>[n_taxa] sigma_site;    

  vector[n_taxa] a_taxon_raw;                 
  real mu_taxon;                              
  real<lower=0> sigma_taxon;             

  real<lower=1e-4> rho_space;         
  real<lower=0> sigma_space;   
  matrix[n_site, n_taxa] delta_raw;       
}

transformed parameters {
  matrix[n_obs, n_taxa] mu;                   
  matrix[n_pred, n_taxa] beta;                
  matrix[n_site, n_taxa] a_site;              
  vector[n_taxa] a_taxon;
  matrix[n_obs, n_taxa] epsilon;              
  matrix[n_site, n_site] Sigma_space;
  matrix[n_site, n_taxa] delta;               

  {
    matrix[n_pred, n_taxa] beta_phylo_raw = beta_raw * L_phylo'; 
    matrix[n_pred, n_pred] L_Sigma_beta = diag_pre_multiply(scale_beta, L_Omega);
    
    for (k in 1:n_pred) {
      for (s in 1:n_taxa) {
        beta[k, s] = mu_beta[k] + dot_product(L_Sigma_beta[k, ], beta_phylo_raw[, s]);
      }
    }
  }

  for (i in 1:n_site) {
    for (j in 1:n_site) {
      if (i == j) {
        Sigma_space[i, j] = square(sigma_space) + (sigma_space * 1e-4) + 1e-5;  
      } else if (is_downstream[i, j] == 1) {
        Sigma_space[i, j] = square(sigma_space) * exp(-square(riv_dist_mat[i, j]) / (2 * square(rho_space)));
      } else {
        Sigma_space[i, j] = 1e-6; 
      }
    }
  }
  
  delta = Sigma_space * delta_raw;
  a_taxon = mu_taxon + a_taxon_raw * sigma_taxon;
  a_site = a_site_raw .* rep_matrix(sigma_site', n_site);
  epsilon = epsilon_raw * sigma_epsilon; 

  for (j in 1:n_taxa) {
    mu[, j] = a_taxon[j] + a_site[site, j] + (u * col(beta, j)) + delta[site, j] + epsilon[, j];
  }
}

model {
  // --- Priors ---
  to_vector(beta_raw) ~ std_normal();
  mu_beta ~ normal(0, 2);
  scale_beta ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(2.0); 

  to_vector(delta_raw) ~ std_normal();
  sigma_space ~ normal(0, 0.5);
  rho_space ~ inv_gamma(2, 0.5); 

  to_vector(a_site_raw) ~ std_normal();   
  sigma_site ~ normal(0, 1.5);

  mu_taxon ~ normal(0, 3);                
  a_taxon_raw ~ std_normal();             
  sigma_taxon ~ normal(0, 1.5);
  
  to_vector(epsilon_raw) ~ std_normal();  
  sigma_epsilon ~ normal(0, 0.5);  

  // --- PARALLELISED LIKELIHOOD ---
  // Grainsize = 1 instructs Stan to stop trying to automate the chunk sizes 
  // dynamically and simply hand out individual species columns to available 
  // processor threads as soon as they become free
  target += reduce_sum(partial_taxa_lpmf, taxa_indices, 1, mu, y);
}

generated quantities {
  // RAM Protection: allocates 0 blocks if save_y_rep is toggled off
  array[save_y_rep ? n_obs : 0, save_y_rep ? n_taxa : 0] int<lower=0, upper=1> y_rep; 
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
        
        // 1. SAFE Posterior predictive data generation
        if (save_y_rep) {
          y_rep[i, j] = bernoulli_rng(prob); 
        }
        
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
