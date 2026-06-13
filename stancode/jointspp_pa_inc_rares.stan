data {
  int<lower=1> n_site;                   // Number of sites
  int<lower=1> n_taxa;                   // Number of species
  int<lower=1> n_pred;                   // Number of fixed predictors
  int<lower=1> n_obs;                    // Total number of samples 
  
  array[n_obs] int<lower=1, upper=n_site> site; // Site index for each sample

  matrix[n_obs, n_pred] u;               // Matrix of scaled predictors for each observation
  matrix[n_taxa, n_taxa] phylo_cor;      // Phylogenetic correlation matrix

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

  // Between-predictor correlation matrix prior
  cholesky_factor_corr[n_pred] L_Omega;           

  // Random site effect shared across all taxa
  vector[n_site] a_site_raw;
  real<lower=0> sigma_site;

  // Taxon intercept
  vector[n_taxa] a_taxon_raw;                 
  real mu_taxon;                              
  real<lower=0> sigma_taxon;

}

transformed parameters {
  matrix[n_obs, n_taxa] mu;                   
  matrix[n_pred, n_taxa] beta;                
  vector[n_site] a_site = a_site_raw * sigma_site;
  // sites share a common site-level intercept (some sites richer than others)
  // avoids overparameterising for rare species
  // could explore a latent factor analysis (ordination) approach to compress 
  // the site matrix if we really want to have taxon-specific site-level intercepts
  vector[n_taxa] a_taxon;

  // 1. Non-Centered Phylogenetic Regression Slopes
  {
  matrix[n_pred, n_taxa] beta_phylo_raw = beta_raw * L_phylo'; 
  matrix[n_pred, n_pred] L_Sigma_beta = diag_pre_multiply(scale_beta, L_Omega);
  
  beta = rep_matrix(mu_beta, n_taxa) + 
                (diag_pre_multiply(scale_beta, L_Omega) * beta_raw * L_phylo');

  // Vectorized calculation of linear predictor matrix
  // u * beta calculates site predictions for all species simultaneously
  // rep_matrix adds the overall taxon intercept across all sites
  // rep_matrix(a_site, n_taxa) adds the shared site quality across all species
  // epsilon has been removed - not needed by Bernoulli to handle overdispersion
  mu = rep_matrix(a_taxon, n_site) + rep_matrix(a_site, n_taxa) + (u * beta);
  
}

model {
  // beta priors
  to_vector(beta_raw) ~ std_normal();
  mu_beta ~ normal(0, 2);
  scale_beta ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2.0); 

  // Random effect Priors
  a_site_raw ~ std_normal();   
  sigma_site ~ normal(0, 1); 

  mu_taxon ~ normal(0, 3);                
  a_taxon_raw ~ std_normal();             
  sigma_taxon ~ normal(0, 1); // tightened to avoid funnel that wasn't a problem with AC
  
  // Vectorized likelihood
  to_array_1d(y) ~ bernoulli_logit(to_vector(mu));

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
