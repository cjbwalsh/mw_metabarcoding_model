data {
  int<lower=1> n_site;                   // Number of sites
  int<lower=1> n_taxa;                   // Number of species
  int<lower=1> n_pred;                   // Number of fixed predictors
  int<lower=1> n_obs;                    // Total number of samples 
  
  array[n_obs] int site;                 // Site index for each sample

  matrix[n_obs,n_pred] u;                // Matrix of scaled, centred predictors for each observation
  matrix[n_taxa, n_taxa] phylo_cor;      // Phylogenetic correlation matrix
  matrix[n_site, n_site] riv_dist_mat;   // River distance matrix (either up or down)
  array[n_obs,n_taxa] int y;             // Presence/absence matrix for each taxon in each sample

}

parameters {
  // Fixed effects
  vector[n_pred] mu_beta;                // Means of fixed-effect beta parameter hyperdistributions
  matrix[n_pred, n_taxa] beta_raw;       // Uncentered regression coefficients

  // Uncorrelated error:
  vector[n_taxa] mu_epsilon;             // Means of fixed-effect beta parameter hyperdistributions
  matrix[n_obs, n_taxa] epsilon_raw;     // uncentered uncorrelated errors
  
  // Between-predictor correlation matrix prior
  cholesky_factor_corr[n_pred] L_Omega;           // Cholesky factor of the correlation matrix
  vector<lower=0>[n_pred] tau;                    // Scale parameter of the correlation matrix

  // site (random) effects (shared across species)
  matrix[n_site,n_taxa] a_site_raw;           // Uncentered random site effect
  vector[n_taxa] mu_site;                     // Mean of hyperdistribution of a_site among taxa

  // taxon intercept
  vector[n_taxa] a_taxon_raw;                 // Uncentered intercept for each taxon
  real mu_taxon;                              // Mean of taxon intercepts hyperdistribution

  // Spatial autocorrelation:
  vector[n_taxa] mu_delta;                     // Means of spatial autocorrelations among taxa
  matrix<lower=0>[n_site,n_taxa] delta_raw;    // rawspatial autocorrelation adjustment
  real<lower=0> alpha;                         // marginal variability of the spatial function
  real<lower=0> rho;                           // length scale parameter, controlling decay with distance

  // Spatial autocorrelation:
  vector[n_taxa] mu_delta;                     // Means of spatial autocorrelations among taxa
  matrix<lower=0>[n_site,n_taxa] delta_raw;    // uncentered spatial autocorrelation adjustment
  // real<lower=0> alpha;                      // marginal variability of the spatial function
  // Leaving alpha out in the first instance, thinking that delta does this job by taxon
  real<lower=0> rho_space;                     // length scale parameter, controlling decay with distance
  real<lower=0> sigma_spatial;                 // Spatial variance

  // Phylogenetic effects
  matrix[n_site, n_taxa] phylo_raw;           // uncentred phylogenetic effects
  real<lower=0> sigma_phylo;                  // Phylogenetic variance
  
}

transformed parameters {
  // Transformed fixed effects
  matrix[n_pred,n_taxa] beta;                // Individual beta parameters of fixed effects in u
  matrix[n_obs, n_taxa] mu;                   // Linear predictor (on logit scale) for each observation and taxon

  // Transformed random effects
  vector[n_taxa] a_taxon;                     // Coefficient of taxon intercept
  matrix[n_site, n_taxa] a_site;              // Coefficient of random site effect
  matrix[n_obs, n_taxa] epsilon;              // Coefficient of random site effect
  matrix[n_site, n_taxa] delta;               // Spatially autocorrelated error term
  
  // --- Non-Centered Phylogenetic Matrix-Variate Regression Slopes ---
  matrix[n_taxa, n_taxa] L_phylo;
  // Compute Cholesky factors for non-centered multivariate effects
  L_phylo = cholesky_decompose(phylo_cor);
  // This propagates both the LKJ cross-predictor correlation AND the phylogenetic tree structure.
  // beta = mu_beta + diag_pre_multiply(scale_beta, L_Omega) * beta_raw * L_phylo'
  {
  matrix[n_pred, n_pred] L_Sigma_beta = diag_pre_multiply(scale_beta, L_Omega);
    
    for (k in 1:n_pred) {
      for (s in 1:n_taxa) {
        beta[k, s] = mu_beta[k] + dot_product(L_Sigma_beta[k, ], beta_raw[, s]);
      }
    }
    // Post-multiply by the transposed Phylogenetic Cholesky factor to induce tree-correlation
    beta = beta * L_phylo';
  }

  // --- Non-centered Spatial Covariance ---
  {
    matrix[n_site,n_site] K_spatial;
    for (i in 1:n_site) {
      for (j in 1:n_site) {
        K_spatial[i, j] = exp(-0.5 * square(riv_dist_mat[i, j] / rho_spatial));
      }
      K_spatial[i, i] = 1.0 + 1e-9; 
    }
    delta = cholesky_decompose(K_spatial) * delta_raw * sigma_spatial;
  }
 
    // Random effects specified with hyperpriors governing the location of the estimate (fixed scale)
  for (j in 1:n_taxa){
    // See https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    // And https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/
    a_taxon[j] = mu_taxon + a_taxon_raw[j];                      // Global taxon intercept
    a_site[,j] = mu_site[j] + a_site_raw[,j];                    // Site effect
    epsilon[,j] = mu_epsilon[j] + epsilon_raw[,j];               // Uncorrelated error term
}

  // Compute the linear predictor for each site and species
  for (j in 1:n_taxa) {
    for (i in 1:n_pred) {
      mu_logit[i,j] = a_taxon[j] + a_site[site[i],j] + u[i,] * beta[,j] + delta[site[i],j] + epsilon[i,j];
    }
  }
}

model {
    // --- Priors for the betas ---
  to_vector(beta_raw) ~ std_normal();
  mu_beta ~ normal(0, 5);
  L_Omega ~ lkj_corr_cholesky(2.0); // Prior on correlation between environmental responses

  to_vector(delta_raw) ~ std_normal();
  sigma_space ~ normal(0, 3);
  
  // Inverse-gamma or Gamma priors work well for length-scale parameters
  // Ensure the prior matches the scale/units of your distance matrix
  rho_space ~ inv_gamma(2.5, 3); 

  // --- Random effects priors
  // Site effect priors
  mu_site ~ normal(0, 5);                 // mean of site effects
  to_vector(a_site_raw) ~ std_normal();   // raw site effect
  
  // Taxon intercept priors
  a_taxon_raw ~ std_normal();             // raw taxon intercept
  mu_taxon ~ std_normal();                // mean of taxon intercepts
  
  // --- Residual effect priors
  // Uncorrelated error term
  to_vector(epsilon_raw) ~ std_normal();  // raw uncorrelated error term
  mu_epsilon ~ std_normal();              // mean of uncorrelated error term
  
   // Estimated taxon occurrences from likelihood estimator
  for (j in 1:n_taxa) {
     // Pull the entire column for taxon j, and match it against the calculated vector mu for taxon j
     // faster than a double loop
      y[, j] ~ bernoulli_logit(mu[, j]); 
          }
 }
