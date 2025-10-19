// --- Multitaxon Joint-Species Distribution Model with spatial autocorrelation
// Correlation specified by species' responses to environmental predictors
// See Ovaskainen & Soininen 2011, Ovaskainen et al. 2016, Wilkinson et al. 2019

data {
  int<lower=0> n_obs;                             // Number of observations (total samples across all sites):
  int<lower=1> n_pred;                            // Number of predictor variables in matrix u
  int<lower=1> n_taxa;                            // Number of taxa in matrix c
  int<lower=1> n_site;                            // Number of sites
  int<lower=1> n_seasons;                         // Number of seasons included in the dataset
  matrix[n_obs,n_pred] u;                         // Matrix of predictor variables for each sample (2 per site)
  array[n_obs,n_taxa] int y;                      // Presence/absence matrix for each taxon in each sample
  array[n_obs,n_taxa] real s;                     // Subsample proportion for each sample (all coarsepicks = 1)
  array[n_obs] int site;                          // Site number for each sample
  array[n_obs] int season;                        // Season for each sample
  matrix[n_site, n_site] M;                       // Haversine distance matrix for sites (for spatial autocorrelation)
}
parameters {
  // fixed-effect slope parameters
  vector[n_pred] mu_gamma;                        // Means of fixed-effect beta parameter hyperdistributions
  matrix[n_pred, n_taxa] gamma_raw;               // Raw individual beta parameters
  
  // Uncorrelated error:
  vector[n_taxa] mu_epsilon;                      // Means of fixed-effect beta parameter hyperdistributions
  matrix[n_obs, n_taxa] epsilon_raw;              // Raw individual beta parameters
  
  // Correlation matrix prior
  cholesky_factor_corr[n_pred] L_Omega;           // Cholesky factor of the correlation matrix
  vector<lower=0>[n_pred] tau;                    // Scale parameter of the correlation matrix

  // Adaptive priors for correlation matrix (alternative to fixed priors)
  // Try these if the model has convergence issues surrounding L_Omega, tau and/or gamma:
//  real<lower=0> sigma_Omega;                      // Hyperprior for scale parameter of Omega
//  real<lower=0> sigma_tau;                        // Hyperprior for scale parameter of tau

  // site effects
  matrix[n_site,n_taxa] a_site_raw;           // Raw coefficient of random site effect
  vector[n_taxa] mu_site;                     // Mean of hyperdistribution of a_site among taxa

  // season effects
  matrix[n_seasons,n_taxa] a_season_raw;      // Raw coefficient of random seasons effect
  vector[n_taxa] mu_season;                   // Mean of hyperdistribution of a_site among taxa
  
  // taxon intercept
  vector[n_taxa] a_taxon_raw;                 // Raw taxon intercept
  real mu_taxon;                              // Mean of taxon intercepts hyperdistribution
  
  // Spatial autocorrelation:
  vector[n_taxa] sp_delta;                    // Raw spatially autocorrelated error term
  vector[n_site] mu_delta;                    // Spatially autocorrelated mean
  real<lower=0> lambda;                       // Spatial exponential decay rate parameter
}
transformed parameters {
  // Transformed fixed effects
  matrix[n_pred,n_taxa] gamma;                // Individual beta parameters of fixed effects in u
  matrix[n_obs, n_taxa] mu;                   // Linear predictor (on logit scale) for each observation and taxon

  // Transformed random effects
  vector[n_taxa] a_taxon;                     // Coefficient of taxon intercept
  matrix[n_site, n_taxa] a_site;              // Coefficient of random site effect
  matrix[n_seasons, n_taxa] a_season;         // Coefficient of random season effect
  matrix[n_obs, n_taxa] epsilon;              // Coefficient of random site effect
  matrix[n_site, n_taxa] delta;               // Spatially autocorrelated error term
  
  // LKJ correlation matrix density
  matrix[n_pred, n_pred] L_Sigma = diag_pre_multiply(tau, L_Omega);
  
  // Inverse distance-weighted matrix
  matrix[n_site, n_site] W = exp(-lambda * M);
  
  // Setting the diagonal to zero to avoid a site being its own neighbor
  for (k in 1:n_site){
    W[k,k] = 0; // Zero diagonal
  }
  
  for (j in 1:n_taxa){
    // Random effects specified with hyperpriors governing the location of the estimate (fixed scale)
    // See https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    // And https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/
    a_taxon[j] = mu_taxon + a_taxon_raw[j];                      // Global taxon intercept
    a_site[,j] = mu_site[j] + a_site_raw[,j];                    // Site effect
    a_season[,j] = mu_season[j] + a_season_raw[,j];              // Season effect
    epsilon[,j] = mu_epsilon[j] + epsilon_raw[,j];               // Uncorrelated error term
    
    // Noncentered parameterisation for gamma based on Christopher-Peterson
    // https://discourse.mc-stan.org/t/joint-species-distribution-model-performance/19799/3
    gamma[,j] = mu_gamma + L_Sigma * gamma_raw[,j]; // Fixed effect slopes for predictors in u
    
    // Spatially autocorrelated error term
    delta[,j] = mu_delta + sp_delta[j]; // Spatially autocorrelated error term
  }
// Non-centered parameterisation for hierarchical models fit with Hamiltonian Monte Carlo
// (e.g.) a_site_raw ~ std_normal() implies a_site ~ normal(0, sigma_site)
// See Neal's Funnel: https://mc-stan.org/docs/stan-users-guide/reparameterization.html

  for (i in 1:n_obs){
    for (j in 1:n_taxa){
      // Model for linear predictor (mu):
      mu[i,j] = a_taxon[j] + u[i,] * gamma[,j] + a_site[site[i],j] + a_season[season[i],j] + epsilon[i,j] + delta[site[i],j];
    }
  }
}
model {
  // --- Priors
  
  // All fixed and random effects are specified with a non-centered parameterization on the individual effects
  // and a standard normal prior on the hyperdistribution mean.
  // Because predictors were scaled prior to modelling, this is equivalent to a weakly informative prior
  // with mean std_normal() and standard deviation 1 on the logit scale for each taxon.
  // This is a more HMC-friendly approach to modelling the fixed effects than the common practice of using
  // flat or very wide priors with an SD hyperdistribution, which can lead to autocorrelation issues.
  
  // --- Fixed effects priors
  
  // Non-centered hierarchical gamma:
  mu_gamma ~ std_normal();                // mean of gamma hyperdistribution
  to_vector(gamma_raw) ~ std_normal();    // raw gamma parameters

  // --- Correlation matrix priors
  // Student's t distribution prior for tau
  // see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // "Aki prefers student_t(3,0,1), something about some shape of some curve, 
  // he put it on the blackboard and I can't remember"
  
  // LKJ prior with fixed scale:
  tau ~ student_t(3,0,1);                 // Scaling parameter of the correlation matrix
  L_Omega ~ lkj_corr_cholesky(2);         // Cholesky factor of the correlation matrix
  // Prior of 2 suggested in McElreath 2022, Statistical Rethinking 2nd ed.
  
  // Alternative scaling priors that might have better convergence properties:
  // Try these if the model has convergence issues surrounding L_Omega, tau and/or gamma:
//  sigma_tau ~ cauchy(1.5, 0.25); // Check posterior - should be greater than ~1
//  tau ~ student_t(3,0,sigma_tau);
//  sigma_Omega ~ cauchy(2.5, 0.25); // Check posterior - should be greater than ~2
//  L_Omega ~ lkj_corr_cholesky(sigma_Omega);

  // --- Random effects priors
  // Site effect priors
  mu_site ~ std_normal();                 // mean of site effects
  to_vector(a_site_raw) ~ std_normal();   // raw site effect
  
  // Seasonal effect priors
  to_vector(a_season_raw) ~ std_normal(); // raw season effect
  mu_season ~ std_normal();               // mean of season effects
  
  // Taxon intercept priors
  a_taxon_raw ~ std_normal();             // raw taxon intercept
  mu_taxon ~ std_normal();                // mean of taxon intercepts
  
  // --- Residual effect priors
  // Uncorrelated error term
  to_vector(epsilon_raw) ~ std_normal();  // raw uncorrelated error term
  mu_epsilon ~ std_normal();              // mean of uncorrelated error term
  
  // Multivariate Conditional Autoregressive model (MCAR) for spatial autocorrelation
  // In Gelfand & Vounatsou 2003 (eq. 4), the autoregressive matrix is defined as
  // tau * (D - lambda * W) where D is a diagonal matrix with the row sums of W
  lambda ~ std_normal();  // spatial exponential decay rate parameter
  for (k in 1:n_site){
    // Spatially autocorrelated mean
      mu_delta[k] ~ normal(W[k,] * mu_delta, 1); // Recursive priors can be solved by HMC
  }
  to_vector(sp_delta) ~ std_normal();  // raw spatially autocorrelated error term
  
  // see https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  // approximation in https://users.aalto.fi/~ave/casestudies/Priors/negbinomial_shape_prior.html

  // Estimation of correlated beta-coefficient parameters in gamma:
  // blocked out - equivalent to non-centered form in transformed parameters block
//   for (j in 1:n_taxa){
//     target += multi_normal_prec_lpdf(gamma[,j] | mu_gamma, quad_form_diag(Omega, tau));
//   }

  // Estimated taxon occurrences from likelihood estimator
  for (i in 1:n_obs){
    for (j in 1:n_taxa){
      target += bernoulli_logit_lpmf(y[i,j] | mu[i,j] + log(s[i,j]));

  // This parameterization adds the marginal log-binomial-probability
  // resulting from subsampling error to the marginal logit-probability
  // of the linear model. It is equivalent to a (50 times) slower
  // parameterization modelling the two marginal binomial probabilities 
  // separately, by looping through all feasible occurrence probabilities
  // given each presence/absence and subsample proportion.
     }
   }
}
generated quantities {
  // log-likelihood for model comparisons (unblock during model development).
  // matrix[n_obs,n_taxa] log_lik; // Log-likelihood for each observation and taxon
  // for (i in 1:n_obs){
  //   for (j in 1:n_taxa){
  //     log_lik[i,j] = bernoulli_logit_lpmf(y[i,j] | mu[i,j] + log(s[i,j]));
  //   }
  // }
  
  // Posterior-predictive distribution
  matrix[n_obs,n_taxa] y_pred; // Predicted occurrence for each observation and taxon
  for (i in 1:n_obs){
    for (j in 1:n_taxa){
      // Posterior-predictive distribution including subsampling error
      y_pred[i,j] = bernoulli_rng(inv_logit(mu[i,j] + log(s[i,j])));
    }
  }
}
