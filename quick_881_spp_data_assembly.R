suppressPackageStartupMessages(library(cmdstanr))
# Appends the threading variable to your Stan makefile configuration
cpp_options <- list(stan_threads = TRUE)
options(cmdstanr_max_num_threads = 4)
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(loo))

# Random seed to fix Stan results
proj_seed <- 12345 
set.seed(proj_seed)

# Network functions used herein: 
source("https://tools.thewerg.unimelb.edu.au/mwstr/data/mwstr_network_functions.R")
# Bug database functions used herein: ct(), 
source("https://tools.thewerg.unimelb.edu.au/mwbugs/data/bug_database_functions.R")
right <- function(x, n) {substr(x, nchar(x) - n + 1, nchar(x))}

####

# Connect to MW databases:
mwbugs <- DBI::dbConnect(RPostgres::Postgres(), dbname = "mwbugs", host = "localhost")
mwstr <- DBI::dbConnect(RPostgres::Postgres(), dbname = "mwstr", host = "localhost")
# Project 3:  Metabarcoded combined-RBA-sample-pair survey of 223 sites across
#   Melbourne 2021-2024 (Nickname MW metabarcoding survey)
samples <- DBI::dbGetQuery(mwbugs, 
                           "SELECT * FROM samples 
                           JOIN sample_project_groups
                           ON samples.smpcode = sample_project_groups.smpcode
                           WHERE sample_project_groups.project_code = 3;")
samples$sitecode[grep("TOO",samples$sitecode)] <- "TOO-4334-2"
samples$sitecode_v12[grep("TOO",samples$sitecode_v12)] <- "TOO_4332c"
sites <- sf::st_read(mwbugs, query = 
                       paste0("SELECT * FROM sites WHERE sitecode_v12 IN ('",
                              paste(samples$sitecode_v12, collapse = "','"),"');"))
# Create an integer site number for use in the Stan model
sites <- dplyr::mutate(sites, site_no = 1:nrow(sites), .before = sitecode)

####

biota_all <- DBI::dbGetQuery(mwbugs, 
                             "SELECT * FROM biota 
                           JOIN sample_project_groups
                           ON biota.smpcode = sample_project_groups.smpcode
                           WHERE sample_project_groups.project_code = 3;")
biota_all <- biota_all[!grepl("EXCLUDE this record from 3-rep", 
                              biota_all$notes),]  #20584
species <- unique(biota_all[c("taxoncode","taxon")])

####

samples$reach_v12 <- sites$reach_v12[match(samples$sitecode, sites$sitecode)]
biota_all$reach <- samples$reach_v12[match(biota_all$smpcode, samples$smpcode)]
biot_by_site <- unique(biota_all[c("taxoncode","taxon","reach")])
biot_prev <- aggregate(biot_by_site$taxoncode, 
                       by = list(taxoncode = biot_by_site$taxoncode, taxon = biot_by_site$taxon), 
                       FUN = length)
biot_prev$prevalence <- biot_prev$x/length(unique(samples$reach_v12))

biota <- biota_all # 20566 records

### Biota presence/absence matrix
biota_ct <- with(biota, ct(smpcode, shortcode, count))

#### standardized phylogenetic matrix phylo_cor

spp_class_all <- readRDS(
  url("https://tools.thewerg.unimelb.edu.au/mwbugs/data/spp_classes_itis.rds"))
species <- species[match(colnames(biota_ct),species$taxoncode),]
spp_class <- spp_class_all[species$taxon]
# Select the subset of taxa relevant for the current analysis (or use all 982 spp)
spp_tree <- taxize::class2tree(spp_class)
# Scale the tree so total root-to-tip height is exactly 1.0
spp_tree_raw <- spp_tree$phylo
max_height <- max(ape::node.depth.edgelength(spp_tree_raw))
scaled_tree <- spp_tree_raw
scaled_tree$edge.length <- scaled_tree$edge.length / max_height
spp_vcv_ou <- ape::vcv(scaled_tree)

####

# Add fields for matching with other tables to the samples table
samples$reach_v12 <- substr(samples$sitecode_v12,1,nchar(samples$sitecode_v12)-1)
# numeric site identifier for stan (site_no)
samples$site_no <- sites$site_no[match(samples$sitecode_v12, sites$sitecode_v12)]

# Derive season and riffle fields from sample specifications
samples$season <- 0; samples$season[lubridate::month(samples$date) > 6] <- 1 
# Season - 0 = Autumn; 1 = Spring
# samples with collection_method mentionnig riffle, include 1 riffle RBA sample
samples$riffle <- 0; 
samples$riffle[grep("riffle", samples$collection_method)] <- 1 

# Reduce the samples table to relevant fields and add sample-level environmental 
# data ("season" and "riffle")
samples <- samples[,c("smpcode","date","sitecode_v12","reach_v12",
                      "site_no", "season", "riffle")]

# Attenuated Forest cover: AF - optimal parameterization from wilkins_etal_2026 
# (See also Walsh (2023) Chapter 7)-  drawn from the mwbugs database,
# which accounts for location of the sample site within the reach
sites_af <- DBI::dbGetQuery(mwbugs, 
                            paste0("SELECT sitecode_v12, af_l5_w633_1750, ",
                                   "af_l5_w633_2018, af_l5_w633_2022 ",
                                   "FROM sitecode_env WHERE sitecode_v12 IN ('",
                                   paste(sites$sitecode_v12, collapse = "', '"),"');"))
# Check all sites are in sites_af: 
# sum(!sites$sitecode_v12 %in% sites_af$sitecode_v12)  # zero missing. all good.
sites$af <- sites_af$af_l5_w633_2022[match(sites$sitecode_v12, sites_af$sitecode_v12)]
# A major storm on 9 June 2021 caused much forest loss https://www.theguardian.com/australia-news/2021/jun/24/dandenong-ranges-power-outage-storm-debris-victoria-emergency-services-warning
# I originally used af_2018 before this date and af_2022 after, but on checking noted that it only changed af values for 
# urban sites for which the forest had no effect (and DNG_7957, where af was slightly higher in 2022 than 2018)
# Therefore leave af as invariate across samples from each site (Which also makes modelling easier)
# samples$af[samples$date < "2021-06-10"] <- 
# sites_af$af_l5_w633_2018[match(samples$sitecode_v12[samples$date < "2021-06-10"],
#                          sites_af$sitecode_v12)]
sites$af_ref <- sites_af$af_l5_w633_1750[match(sites$sitecode_v12, sites_af$sitecode_v12)] 
# reference state for AF to be used for predictions from the ultimate model

# Impervious cover and other environmental variables from the mwstr database
cat_env <- DBI::dbGetQuery(mwstr,
                           paste0("SELECT subcs.reach, subcs.carea_km2, ei_l44_sw1942_2022, c_basalt, ", 
                                  "ti_2022, tf_2018, tf_2022, meanq_mm, dor_2010 FROM cat_env ",
                                  "JOIN subcs ON cat_env.site = subcs.site ","WHERE cat_env.reach IN ('",
                                  paste(sites$reach_v12, collapse = "', '"), "');"))

subc_env <- DBI::dbGetQuery(mwstr,
                            paste0("SELECT reach, meant_y30_2022, elev_asl_m FROM subc_env ",
                                   "WHERE reach IN ('",
                                   paste(sites$reach_v12, collapse = "', '"), "');"))
# Catchment area in sq km. 
sites$carea_km2 <- cat_env$carea_km2[match(sites$reach_v12, cat_env$reach)] 
# Proportion of catchment with basalt geology See Walsh (2023) Chapter 5.7
sites$basalt <- cat_env$c_basalt[match(sites$reach_v12, cat_env$reach)]
# Degree of regulation: see Walsh (2023) Chapter 9
sites$dor <- cat_env$dor_2010[match(sites$reach_v12,cat_env$reach)] 
# Effective imperviousness: see Walsh (2023) Chapter 6
sites$ei <- cat_env$ei_l44_sw1942_2022[match(sites$reach_v12,cat_env$reach)] 
# Total imperviousness:  see Walsh (2023) Chapter 6
sites$ti <- cat_env$ti_2022[match(sites$reach_v12,cat_env$reach)] 
# Total forest cover:  see Walsh (2023) Chapter 7
sites$tf <- cat_env$tf_2022[match(sites$reach_v12,cat_env$reach)] 
# Mean annual runoff depth (mm). See Walsh (2023) Chapter 8
sites$meanq <- cat_env$meanq_mm[match(sites$reach_v12,cat_env$reach)] 
# Mean annual air temperature (degrees C). See Walsh (2023) Chapter 5.8
sites$meant <- subc_env$meant_y30_2022[match(sites$reach_v12,subc_env$reach)] 
# Elevation above sea-level (m). See Walsh (2023) Chapter 5.9
sites$elev <- subc_env$elev_asl_m[match(sites$reach_v12,subc_env$reach)]
# A lowland/estuarine influenc indicator
sites$lowness <- -1
sites$lowness[sites$elev < 30] <- -0.5
sites$lowness[sites$elev < 15] <- 0
sites$lowness[sites$elev < 7.5] <- 0.5
sites$lowness[sites$elev < 2.5] <- 1

# Transformations and interactions prior to scaling
sites$lcarea <- log(sites$carea_km2)
sites$smeanq <- sites$meanq^0.5
sites$dor1 <- pmin(sites$dor,1)
sites$lei <- log(sites$ei + 0.001)
sites$q_af <- sites$smeanq * sites$af
sites$q_ei <- sites$smeanq * sites$lei
sites$b_ei <- sites$basalt * sites$lei
sites$q_d <- log(sites$smeanq * sites$dor1 + 0.1)
sites$af_ei <- sites$af * sites$lei

####

biota_ct <- biota_ct[match(samples$smpcode,row.names(biota_ct)),] # ensure c and u rows are ordered the same way
# make first two well-predicted prevalent species to stabilize the lambda_diag constraint

# Ensure samples and biota_1_ct rows are in the same order
samples <- samples[match(row.names(biota_ct), samples$smpcode),]
sum(!samples$smpcode == row.names(biota_ct))  # should be zero

# Index for sites in samples table
site <- samples$site_no 

### Two predictor matrices: one for site-specific variables (u_site), 
###                         one for sample(obs)-specific variables (u_obs)

# scale fixed predictor variables in 
# keeping them as scaled objects for ease of back-transforming later
I <- scale(sites$lei)
F <- scale(sites$af)
Q <- scale(sites$smeanq)
Q2 <- scale(sites$smeanq^2)
C <- scale(sites$lcarea)
T <- scale(sites$meant)
D <- scale(sites$dor1)
B <- scale(sites$basalt)
Q_F <- scale(sites$q_af)
Q_I <- scale(sites$q_ei)
Q_D <- scale(sites$q_d)
F_I <- scale(sites$af_ei)
B_I <- scale(sites$b_ei)
L <- sites$lowness

# Save scaling attributes for all scaled parameters
pred_specs <- data.frame(db_name = c("ei_l44_sw1942_2022","af_l5_w633_2022","meanq_mm","meanq_mm","carea_km2","dor_2010","meant_y30_2022",
                                     "elev_asl_m","basalt","ei*af","meanq*af","meanq*ei","basalt*ei","meanq*dor"),
                         short_name = c("ei","af","meanq","meanq2","carea_km2","dor","meant","lowness","c_basalt",
                                        "if","qf","qi","bi","qd"),
                         transf_name = c("lei","af","smeanq","meanq","lcarea","dor1","meant","lowness","basalt",
                                         "af_ei","q_af","q_ei","b_ei","q_d"),
                         transformation = c("log(x+0.001)",NA,"x^0.5",NA,"log(x)","pmin(x,1)",NA,"bespoke scale",NA,"lei*af","smeanq*af","smeanq*lei","basalt*lei","log(smeanq*dor1+0.1"),
                         scaled_name = c("I","F","Q","Q2","C","D","T","L","B","F_I","Q_F","Q_I","B_I","Q_D"),
                         scale = NA, center = NA)

for(i in 1:nrow(pred_specs)){
  pred_specs$scale[i] <- ifelse(pred_specs$scaled_name[i] %in% c("B","L"), NA,
                                attr(get(pred_specs$scaled_name[i]), "scaled:scale"))
  pred_specs$center[i] <- ifelse(pred_specs$scaled_name[i] %in% c("B","L"), NA,
                                 attr(get(pred_specs$scaled_name[i]), "scaled:center"))
}

# Save the scaled parameters for use in predictions
write.csv(pred_specs, "~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/parameter_scaling_attributes.csv")

u_site <- unique(as.matrix(data.frame(I = I,      # Effective imperviousness
                                      F = F,                   # Attenuated Forest cover
                                      Q = Q,                   # mean annual runoff
                                      Q2 = Q2,                 # square of Q for modal response
                                      C = C,                   # catchment area
                                      D = D,                   # Degree of regulation
                                      T = T,                   # Temperature 
                                      L = L,                   # lowness
                                      B = B,                   # Basalt
                                      F_I = F_I,                # Interaction betwen I and F
                                      Q_F = Q_F,                # Interaction betwen Q and F
                                      Q_I = Q_I)))             # Interaction between Q and D
# diff(site) # all zeros and 1s meaning that unique(u_site) is ordered by site
u_obs <- as.matrix(data.frame(riff = samples$riffle, season = samples$season))

# Matrix for modelling latent variation unrelated to predictors
# Compute the Hat projection matrix onto the site covariate space
H_site <- u_site %*% solve(t(u_site) %*% u_site) %*% t(u_site)
# Compute the Residual Projection Matrix (The Null-Space Operator)
M_site_res <- diag(nrow(u_site)) - H_site

# Assemble biota and predictor data for post-processing purposes
species_set <- colnames(biota_ct)
species_set <- biot_prev[match(species_set, biot_prev$taxoncode),]
sample_set <- row.names(biota_ct)
sample_set <- samples[match(sample_set, samples$smpcode),c("smpcode","reach_v12","site_no")]

# Assemble stan data for second trial set
master_data <- list(n_obs = nrow(biota_ct),          # no. samples
                    n_site_pred = ncol(u_site),           # no. site-specific predictors
                    n_obs_pred = ncol(u_obs),             # no. sample-specific predictors
                    n_taxa = ncol(biota_ct),            # no. taxa
                    n_site = max(site),                   # no. sites with a site upstream
                    n_latent = 2,                         # no. latent factors for site effect
                    u_site = u_site,                      # site-specific predictor matrix
                    u_obs = u_obs,                        # sample-specific predictor matrix
                    y = as.matrix(biota_ct),
                    site  = site,
                    M_site_res = M_site_res,
                    phylo_cor = spp_vcv_ou)

params_to_summarise <- c(
  "mu_beta_site", "scale_beta_site", "beta_site",  # Fixed site-level effects
  "mu_beta_obs", "scale_beta_obs", "beta_obs",  # Fixed sample-level effects
  "mu_taxon", "sigma_taxon",        # Species hyper-parameters
  "z_expanded", "lambda",          # latent factors - these do not converge as they are unconstrained
  "a_taxon", "alpha_phylo",           # Random species effect, strength of phylogeny effect
  "tjurs_r2_total",                 # Explanatory power using fixed predictors and latent random site variation
  "tjurs_r2_env",                   # Explanatory power solely using fixed predictors
  "y_rep"                           # For observed vs predicted plots
) # to reduce the RAM required to build the summary table

pred_site <- c("I","F","Q","Q2","C","D","T","B","L","F_I","Q_F","Q_I")
pred_obs <- c("riff","season") # included in all candidate models

run_881spp_model <- function(pred_site, mod_code, mod_path = "stancode/jointspp_pa_no_eps_ou.stan"){
  pred_set <- u_site[,match(pred_site, colnames(u_site))]
  mod_data <- master_data
  mod_data$n_site_pred <- length(pred_site)
  mod_data$u_site <- pred_set
  mod_bundle <- list(species_set = species_set, sample_set = sample_set, 
                     pred_set = pred_set, data = mod_data)
  mod <- cmdstan_model(stan_file = mod_path,
                       compile_model_methods = TRUE) #, 
  #  force_recompile = TRUE )
  mod_dir <-paste0("~/uomShare/wergStaff/ChrisW/git-data/",
                   "mw_metabarcoding_model/full_881spp_models/", mod_code)
  if (!dir.exists(mod_dir)) dir.create(mod_dir)
  model_fit <- mod$sample(data = mod_data, seed = proj_seed,
                          output_dir = mod_dir,
                          chains = 4, 
                          parallel_chains = 4,
                          iter_warmup = 500,  
                          iter_sampling = 500)
  mod_bundle$fit <- model_fit
  mod_bundle$diagnostics <- model_fit$diagnostic_summary()
  saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  mod_summary <- model_fit$summary(variables = params_to_summarise); gc()
  mod_bundle$summary <- mod_summary
  saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  # The following is approach is taken instead of using cmdstanr's wrapper fit$loo()
  # to avoid spikes in RAM usage.
  # (ess_bulk is used because elpd is an expected mean. bulk ess is a measure of 
  # sampling efficiency of central tendenciees)
  log_lik_matrix <- model_fit$draws(variables = "log_lik", format = "draws_matrix")
  n_chains <- posterior::nchains(log_lik_matrix)
  n_draws <- posterior::ndraws(log_lik_matrix)
  r_eff <- posterior::ess_bulk(log_lik_matrix) / (n_draws / n_chains)
  mod_loo <- loo::loo(log_lik_matrix, r_eff = r_eff, cores = 1) # 4 cores sails too close to the wind o 64 Gb RAM machine
  mod_bundle$loo <- mod_loo
  saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  rm(mod_bundle,mod_summary, model_fit, mod_loo, log_lik_matrix, r_eff); gc()
}
