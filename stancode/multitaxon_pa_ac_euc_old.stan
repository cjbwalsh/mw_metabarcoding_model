// --- Multitaxon Joint-Species Distribution Model with spatial autocorrelation
// Correlation specified by species' responses to environmental predictors
// See Ovaskainen & Soininen 2011, Ovaskainen et al. 2016, Wilkinson et al. 2019

data {
  int<lower=0> n_obs;                             // Number of observations (total samples across all sites):
  int<lower=1> n_pred;                            // Number of predictor variables in matrix u
  int<lower=1> n_taxa;                            // Number of taxa in matrix c
  int<lower=1> n_site;                            // Number of sites
  matrix[n_obs,n_pred] u;                         // Matrix of predictor variables for each sample (2 per site)
  array[n_obs,n_taxa] int y;                      // Presence/absence matrix for each taxon in each sample
  array[n_obs] int site;                          // Site number for each sample
//  matrix M[n_site, n_site];                       // matrix of distances between sites in km (capped at 100)
  array[n_site] vector[2] coords;                 // coordinates of the sites (this won't work for the network distances)
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

  // taxon intercept
  vector[n_taxa] a_taxon_raw;                 // Raw taxon intercept
  real mu_taxon;                              // Mean of taxon intercepts hyperdistribution
  
  // Spatial autocorrelation:
  vector[n_taxa] mu_delta;                     // Means of spatial autocorrelations among taxa
  matrix<lower=0>[n_site,n_taxa] delta_raw;    // rawspatial autocorrelation adjustment
//  real<lower=0> alpha;                         // marginal variability of the spatial function
  real<lower=0> rho;                           // length scale parameter, controlling decay with distance
}

transformed parameters {
  // Transformed fixed effect
  matrix[n_pred,n_taxa] gamma;                // Individual beta parameters of fixed effects in u
  matrix[n_obs, n_taxa] mu;                   // Linear predictor (on logit scale) for each observation and taxon

  // Transformed random effects
  vector[n_taxa] a_taxon;                     // Coefficient of taxon intercept
  matrix[n_site, n_taxa] a_site;              // Coefficient of random site effect
  matrix[n_obs, n_taxa] epsilon;              // Coefficient of random site effect
  matrix[n_site, n_taxa] delta;               // Spatially autocorrelated error term
  
  // LKJ correlation matrix density
  matrix[n_pred, n_pred] L_Sigma = diag_pre_multiply(tau, L_Omega);
  
  // Spatial autocorrelation effect
  // squared exponential decay (sharper shift at zero than exponential decay) after 
  // https://www.r-bloggers.com/2023/11/using-stan-to-model-geostatistical-count-data-with-distance-matrices/
  real d_adj = 1e-5;          // Jitter added to distance matrix diagonal for stability
  matrix[n_site, n_site] S_Sigma;
  matrix[n_site, n_site] L_S_Sigma;
  // Transform distances into the spatial covariance matrix
  # In the following, I'm trialling fixed 1 instead of parameter alpha
   S_Sigma = gp_exp_quad_cov(coords, 1, rho);//squared exponential process (failed to converge with alpha fixed or not)
  //  S_Sigma = gp_exponential_cov(coords, 1, rho);  // 1st-order exponential process (aka Ornstein-Uhlenbeck)
  // 2. Ensure the diagonal elements are strictly positive
    for (i in 1:n_site) {
      S_Sigma[i, i] = S_Sigma[i, i] + d_adj;
    }
  // 3. Compute the lower-triangular Cholesky factor
    L_S_Sigma = cholesky_decompose(S_Sigma);
   // 
   // matrix[n_site, n_site] L_space = cholesky_decompose( //square(alpha) * 
   //                                       exp(-0.5 * (square(M / rho))) +
   //                                             diag_matrix(rep_vector(d_adj, n_site)));
// the following did not work well
 // I used pow(M, lambda) because M is supplied as exp(-distance) and 
  // and exponential weighting is usually exp(-lambda * distance)
//matrix[n_site, n_site] L_space = cholesky_decompose(pow(M, 0.25))); 
#trial setting lambda = 0.25 to get model to work HDD ~ 2.5 km
  //                                +  diag_matrix(rep_vector(1 - 0.25, n_site)));

 for (j in 1:n_taxa){
    // Random effects specified with hyperpriors governing the location of the estimate (fixed scale)
    // See https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
    // And https://statmodeling.stat.columbia.edu/2018/04/03/justify-my-love/
    a_taxon[j] = mu_taxon + a_taxon_raw[j];                      // Global taxon intercept
    a_site[,j] = mu_site[j] + a_site_raw[,j];                    // Site effect
    epsilon[,j] = mu_epsilon[j] + epsilon_raw[,j];               // Uncorrelated error term
    
    // Noncentered parameterisation for gamma based on Christopher-Peterson
    // https://discourse.mc-stan.org/t/joint-species-distribution-model-performance/19799/3
    gamma[,j] = mu_gamma + L_Sigma * gamma_raw[,j]; // Fixed effect slopes for predictors in u
    delta[,j] = mu_delta[j] + L_S_Sigma * delta_raw[,j];
        }
// Non-centered parameterisation for hierarchical models fit with Hamiltonian Monte Carlo
// (e.g.) a_site_raw ~ std_normal() implies a_site ~ normal(0, sigma_site)
// See Neal's Funnel: https://mc-stan.org/docs/stan-users-guide/reparameterization.html

  for (i in 1:n_obs){
    for (j in 1:n_taxa){
      // Model for linear predictor (mu):
      mu[i,j] = a_taxon[j] + u[i,] * gamma[,j] + a_site[site[i],j] + delta[site[i],j] + epsilon[i,j];
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
  
  // Taxon intercept priors
  a_taxon_raw ~ std_normal();             // raw taxon intercept
  mu_taxon ~ std_normal();                // mean of taxon intercepts
  
  // --- Residual effect priors
  // Uncorrelated error term
  to_vector(epsilon_raw) ~ std_normal();  // raw uncorrelated error term
  mu_epsilon ~ std_normal();              // mean of uncorrelated error term
  
  // spatial autocorrelation
  // prior to block values <1 km while leaving a long right tail up to ~50 km
  // for squared exponential decay
  rho ~ inv_gamma(2.5, 3); //
  // For 1st-order exponential decay, the following better covers the range of distances in the dataset
  // rho ~ inv_gamma(6, 10);  // bounds to range of neighbour distances, with bulk of sampling mass b/n 2 and 25 
  // Prior for the GP covariance magnitude
#  alpha ~ std_normal();
  to_vector(mu_delta) ~ std_normal();
  
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
      target += bernoulli_logit_lpmf(y[i,j] | mu[i,j]);

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
      y_pred[i,j] = bernoulli_rng(inv_logit(mu[i,j]));
    }
  }
}
