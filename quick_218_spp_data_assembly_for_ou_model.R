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
# Remove problem UYT, LOL and SHE samples (problem source data)
samples <- samples[!samples$smpcode %in% c("384-UYT-199-4-DP","384-LOL-11922-8-EP","384-SHE-1665-6-DP"),]
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
sum(grepl("; likely false positive", biota_all$notes))  #Remove the 2182 likely false positives
biota_all <- biota_all[-grep("; likely false positive", biota_all$notes),] # reduced to 18186
# Exclude all 
biota_all <- biota_all[biota_all$smpcode %in% samples$smpcode,]
species <- unique(biota_all[c("taxoncode","taxon")])  #878 species

####

samples$reach_v12 <- sites$reach_v12[match(samples$sitecode, sites$sitecode)]
biota_all$reach <- samples$reach_v12[match(biota_all$smpcode, samples$smpcode)]
biot_by_site <- unique(biota_all[c("taxoncode","taxon","reach")])
biot_prev <- aggregate(biot_by_site$taxoncode, 
                       by = list(taxoncode = biot_by_site$taxoncode, taxon = biot_by_site$taxon), 
                       FUN = length)
biot_prev$prevalence <- biot_prev$x/length(unique(samples$reach_v12))
max(biot_prev$prevalence) # 0.84 is the highest prevalence
sum(biot_prev$prevalence >= 0.2)  #74 spp with prevalence > 0.2
sum(biot_prev$prevalence >= 0.10 & biot_prev$prevalence < 0.2) # 134 sp with prevalence 0.1-0.2
sum(biot_prev$prevalence >= 0.05 & biot_prev$prevalence < 0.1) # 149 sp with prevalence 0.05-0.1
sum(biot_prev$prevalence >= 0.01 & biot_prev$prevalence < 0.05) # 300 sp with prevalence 0.01-0.05 (0.01 = 3 or more records)
# Given more information in more prevalent spp, let's take all of group 1, half of group 2, a third of group 3 and a twelfth of group 4
# Build a ~200-spp dataset for mode development
# 87 + 134/2 + 145/3 + 294/12 = 226 spp # 5 specials added below...
group_1 <- biot_prev[biot_prev$prevalence >= 0.2,]  # 74 species
group_2 <- biot_prev[biot_prev$prevalence >= 0.10 & biot_prev$prevalence < 0.2,]  # 134 species
group_3 <- biot_prev[biot_prev$prevalence >= 0.05 & biot_prev$prevalence < 0.1,]  # 149 species
group_4 <- biot_prev[biot_prev$prevalence >= 0.01 & biot_prev$prevalence < 0.05,]  # 300 species

spp_trial <- rbind(group_1, group_2[seq(2,134,2),], group_3[seq(3,145,3),], group_4[seq(12,294,12),])  # 213
# sum(duplicated($taxoncode=) # 0 # all good
# I also want to include some lowland species to give the "L" predictor a fair go
lowland_spp <- c("OP0599A1",  # 1 record - Cannings St ford
                 "IB020101","KG021202") # both with lots of lowland records
# sum(lowland_spp %in%spp_trial$taxoncode )  # None included so far, so add them
# Also the two Pontogen_genus_A species from the Dandenongs to bulk up congenerics
pont_spp <- c("OP03A1A1","OP03A1A3")
# sum(pont_spp %in%spp_trial$taxoncode )  # None included so far, so add them
spp_trial<- unique(rbind(spp_trial, biot_prev[biot_prev$taxoncode %in% c(lowland_spp,pont_spp),]))  #218 spp

biota_1 <- biota_all[biota_all$taxoncode %in% spp_trial$taxoncode,]  #10986 records

### Biota presence/absence matrix
biota_1_ct <- with(biota_1, ct(smpcode, shortcode, count))

#### standardized phylogenetic matrix phylo_cor
spp_class_all <- readRDS("~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/spp_class_all_aug2026.rds")
# Updated with revised Chironomidae and Pontogeneiidae.  need to update the https version as well....
#   url("https://tools.thewerg.unimelb.edu.au/mwbugs/data/spp_classes_itis.rds"))
# the non-https version removes subgenus, subtribe and section, which are non-informative
spp_trial <- spp_trial[match(colnames(biota_1_ct),spp_trial$taxoncode),]
spp_class_200 <- spp_class_all[spp_trial$taxon]
# Select the subset of taxa relevant for the current analysis (or use all 982 spp)
spp_tree_1 <- taxize::class2tree(spp_class_200)
# Scale the tree so total root-to-tip height is exactly 1.0
spp_tree_1_raw <- spp_tree_1$phylo
max_height <- max(ape::node.depth.edgelength(spp_tree_1_raw))
scaled_tree <- spp_tree_1_raw
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

miseq_samples <-DBI::dbGetQuery(mwbugs, "SELECT DISTINCT smpcode, miseq_run FROM miseq_samples;")
samples$run <- miseq_samples$miseq_run[match(samples$smpcode, miseq_samples$smpcode)]

# Reduce the samples table to relevant fields and add sample-level environmental 
# data ("season" and "riffle")
samples <- samples[,c("smpcode","date","sitecode_v12","reach_v12",
                      "site_no", "season", "riffle","run")]

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
# A lowland/estuarine influence indicator
sites$lowness <- -1
sites$lowness[sites$elev < 30] <- -0.5
sites$lowness[sites$elev < 15] <- 0
sites$lowness[sites$elev < 7.5] <- 0.5
sites$lowness[sites$elev < 2.5] <- 1
# An upstream habitat indicator from 'calculate length sampleable upstream.R'
sampleable_us <- readRDS("~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/us_sampleable_l_km.rds")
sites$us_sampleable_l_km <- sampleable_us$us_sampleable_l_km[match(sites$reach_v12, sampleable_us$reach)]
# # Transformations prior to scaling
sites$lcarea <- log(sites$carea_km2)
sites$smeanq <- sites$meanq^0.5
sites$dor1 <- pmin(sites$dor,1)
sites$lei <- log(sites$ei + 0.001)
# I had a play around with this, but it was too correlated with C.
# sites$lus <- 1
# sites$lus[sites$us_sampleable_l_km < 10] <- 0.5
# sites$lus[sites$us_sampleable_l_km < 5] <- 0
# sites$lus[sites$us_sampleable_l_km < 2] <- -0.5
# sites$lus[sites$us_sampleable_l_km < 0.5] <- -1

####

biota_1_ct <- biota_1_ct[match(samples$smpcode,row.names(biota_1_ct)),] # ensure c and u rows are ordered the same way
# make first two well-predicted prevalent species to stabilize the lambda_diag constraint

# Ensure samples and biota_1_ct rows are in the same order
samples <- samples[match(row.names(biota_1_ct), samples$smpcode),]
sum(!samples$smpcode == row.names(biota_1_ct))  # should be zero

# Index for sites in samples table
site <- samples$site_no 

### Two predictor matrices: one for site-specific variables (u_site), 
###                         one for sample(obs)-specific variables (u_obs)

# scale fixed predictor variables in 
# keeping them as scaled objects for ease of back-transforming later
I <- scale(sites$lei)
F <- scale(sites$af)
Q <- scale(sites$smeanq)
Q2 <- scale(sites$smeanq)^2
C <- scale(sites$lcarea)
T <- scale(sites$meant)
T2 <- scale(sites$meant)^2
D <- scale(sites$dor1)
B <- scale(sites$basalt)
L <- sites$lowness
Q_F <- Q * F  
Q_I <- Q * I  
Q_D <- Q * D
F_I <- F * I
B_I <- B * I
F_I_Q <- F * I * Q

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}
pairs(data.frame(I = I, F = F, Q = Q, Q2 = Q2, C = C, T = T, T2 = T2, D = D,
                 B = B, L = L, Q_F = Q_F, Q_I = Q_I, Q_D = Q_D, 
                 F_I = F_I, B_I = B_I, F_I_Q = F_I_Q),
      lower.panel = panel.smooth, upper.panel = panel.cor, gap=0, row1attop=FALSE)
# Multicollinearity not too bad. Most correlated is Q2 and Q_F at 0.73

# Save scaling attributes for all scaled parameters
pred_specs <- data.frame(db_name = c("ei_l44_sw1942_2022","af_l5_w633_2022","meanq_mm","carea_km2","dor_2010","meant_y30_2022",
                                     "elev_asl_m","basalt"),
                         short_name = c("ei","af","meanq","carea_km2","dor","meant","lowness","c_basalt"),
                         transf_name = c("lei","af","smeanq","lcarea","dor1","meant","lowness","basalt"),
                         transformation = c("log(x+0.001)","x^0.5",NA,"log(x)","pmin(x,1)",NA,"bespoke scale",NA),
                         scaled_name = c("I","F","Q","C","D","T","L","B"),
                         scale = NA, center = NA)

for(i in 1:nrow(pred_specs)){
  pred_specs$scale[i] <- ifelse(pred_specs$scaled_name[i] %in% c("B","L"), NA,
                                attr(get(pred_specs$scaled_name[i]), "scaled:scale"))
  pred_specs$center[i] <- ifelse(pred_specs$scaled_name[i] %in% c("B","L"), NA,
                                 attr(get(pred_specs$scaled_name[i]), "scaled:center"))
}

# Save the scaled parameters for use in predictions
write.csv(pred_specs, "~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/parameter_scaling_attributes_v2.csv")

u_site <- unique(as.matrix(data.frame(I = I,      # Effective imperviousness
                                      F = F,                   # Attenuated Forest cover
                                      Q = Q,                   # mean annual runoff
                                      Q2 = Q2,                 # square of Q for modal response
                                      C = C,                   # catchment area
                                      D = D,                   # Degree of regulation
                                      T = T,                   # Temperature 
                                      T2 = T2,                 # square of T for modal response
                                      L = L,                   # lowness
                                      B = B,                   # Basalt
                                      Q_D = Q_D,                # Interaction betwen I and F
                                      F_I = F_I,                # Interaction betwen I and F
                                      Q_F = Q_F,                # Interaction betwen Q and F
                                      Q_I = Q_I,
                                      F_I_Q = F_I_Q)))             # Interaction between Q and D
# diff(site) # all zeros and 1s meaning that unique(u_site) is ordered by site

# ...OK here I need to add 3 binary variables
# mb = bulk-processing experiment (mix of subsamples to make 3 reps) (run 13)
# ma = bulk processed run with many false positives (run 15)
# mc = bulk processed runs with many false negatives (runs 19 and 22)
# cases where all of the above are zero = bulk processed runs with moderate level of fps (runs 16,17,18,20)

u_obs <- as.matrix(data.frame(riff = samples$riffle, season = samples$season,
                              mb = as.numeric(samples$run == 13),
                              ma = as.numeric(samples$run %in% 15:16),
                              mc = as.numeric(samples$run %in% c(19,22))))

# Matrix for modelling latent variation unrelated to predictors
# Compute the Hat projection matrix onto the site covariate space
H_site <- u_site %*% solve(t(u_site) %*% u_site) %*% t(u_site)
# Compute the Residual Projection Matrix (The Null-Space Operator)
M_site_res <- diag(nrow(u_site)) - H_site

# Assemble biota and predictor data for post-processing purposes
species_set <- colnames(biota_1_ct)
species_set <- biot_prev[match(species_set, biot_prev$taxoncode),]
sample_set <- row.names(biota_1_ct)
sample_set <- samples[match(sample_set, samples$smpcode),c("smpcode","reach_v12","site_no")]

# Assemble stan data for second trial set
master_data <- list(n_obs = nrow(biota_1_ct),          # no. samples
                    n_site_pred = ncol(u_site),           # no. site-specific predictors
                    n_obs_pred = ncol(u_obs),             # no. sample-specific predictors
                    n_taxa = ncol(biota_1_ct),            # no. taxa
                    n_site = max(site),                   # no. sites with a site upstream
                    n_latent = 2,                         # no. latent factors for site effect
                    u_site = u_site,                      # site-specific predictor matrix
                    u_obs = u_obs,                        # sample-specific predictor matrix
                    y = as.matrix(biota_1_ct),
                    site  = site,
                    M_site_res = M_site_res,
                    phylo_cor = spp_vcv_ou)            

params_to_summarise <- c(
  "mu_beta_site", "scale_beta_site", "beta_site",  # Fixed site-level effects
  "mu_beta_obs", "scale_beta_obs", "beta_obs",  # Fixed sample-level effects
  "mu_taxon", "sigma_taxon",        # Species hyper-parameters
  "z_expanded", "lambda",          # latent factors - these do not converge as they are unconstrained
  "a_taxon", "alpha_phylo",           # Random species effect, strength of phylogeny effect
  #  "tjurs_r2_total",                 # Explanatory power using fixed predictors and latent random site variation
  #  "tjurs_r2_env",                   # Explanatory power solely using fixed predictors
  "y_rep"                           # For observed vs predicted plots
) # to reduce the RAM required to build the summary table

pred_site <- c("I","F","Q","Q2","C","D","T","T2","B","L","Q_D","F_I","Q_F","Q_I","F_I_Q")
pred_obs <- c("mb","mc","ma","riff","season") # included in all candidate models
proj_dir <- "~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/trial2_218spp_models/"

run_218spp_model <- function(pred_site, mod_code, 
                             proj_dir = "~/uomShare/wergStaff/ChrisW/git-data/mw_metabarcoding_model/trial2_218spp_models/",
                             mod_path = "stancode/jointspp_pa_mb.stan",
                             n_iterations = 500, n_chains = 4,
                             compile_summary = FALSE, 
                             compile_loo = FALSE){
  pred_set <- u_site[,match(pred_site, colnames(u_site))]
  mod_data <- master_data
  mod_data$n_site_pred <- length(pred_site)
  mod_data$u_site <- pred_set
  mod_bundle <- list(species_set = species_set, sample_set = sample_set, 
                     pred_set = pred_set, data = mod_data)
  mod <- cmdstan_model(stan_file = mod_path,
                       compile_model_methods = TRUE) #, 
  #  force_recompile = TRUE )
  mod_dir <-paste0(proj_dir, mod_code)
  if (!dir.exists(mod_dir)) dir.create(mod_dir)
  model_fit <- mod$sample(data = mod_data, seed = proj_seed,
                          output_dir = mod_dir,
                          chains = n_chains, 
                          parallel_chains = n_chains,
                          iter_warmup = n_iterations,  
                          iter_sampling = n_iterations)
  mod_bundle$fit <- model_fit
  mod_bundle$diagnostics <- model_fit$diagnostic_summary()
  saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  # The following is approach is taken instead of using cmdstanr's wrapper fit$loo()
  # to avoid spikes in RAM usage.
  # (ess_bulk is used because elpd is an expected mean. bulk ess is a measure of 
  # sampling efficiency of central tendenciees)
  if(compile_loo){
    log_lik_matrix <- model_fit$draws(variables = "log_lik_env", format = "draws_matrix")
    n_chains <- posterior::nchains(log_lik_matrix)
    n_draws <- posterior::ndraws(log_lik_matrix)
    r_eff <- posterior::ess_bulk(log_lik_matrix) / (n_draws / n_chains)
    mod_loo <- loo::loo(log_lik_matrix, r_eff = r_eff, cores = 1) # 4 cores sails too close to the wind o 64 Gb RAM machine
    mod_bundle$loo <- mod_loo
    saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  }
  if(compile_summary){
    mod_summary <- model_fit$summary(variables = params_to_summarise); gc()
    raw_z <- model_fit$draws("z_expanded", format = "array")
    raw_lambda <- model_fit$draws("lambda_raw", format = "array")
    first_vector <- raw_z[dim(raw_z)[1], dim(raw_z)[2], ]
    target_matrix <- matrix(NA, nrow = mod_data$n_obs, ncol = mod_data$n_latent)
    target_matrix[, 1] <- first_vector[1:mod_data$n_obs]           # Elements 1 to 274 = Axis 1
    target_matrix[, 2] <- first_vector[(mod_data$n_obs + 1):(mod_data$n_obs*2)]    # Elements 275 to 548= Axis 2
    target_matrix <- scale(target_matrix)
    # Create empty, 4D arrays to store the aligned results
    # dim =iterations, chains, observations, latent axes
    z_unpacked_aligned <- array(NA, dim = c(n_iterations, n_chains, mod_data$n_obs, mod_data$n_latent))
    lambda_unpacked_aligned <- array(NA, dim = c(n_iterations, n_chains, mod_data$n_taxa, mod_data$n_latent))
    for (c in 1:n_chains) {
      for (i in 1:n_iterations) {
        # current 550-element vector
        current_vector <- raw_z[i ,c ,]
        # Create an empty  275 x 2 matrix to be populated
        current_matrix <- matrix(NA, nrow = mod_data$n_obs, ncol = mod_data$n_latent)
        current_matrix[, 1] <- current_vector[1:mod_data$n_obs]            # Dimension 1 (Axis 1)
        current_matrix[, 2] <- current_vector[(mod_data$n_obs + 1):(mod_data$n_obs*2)]    # Dimension 2 (Axis 2)
        # Standardise to remove cross-chain scale discrepancies
        current_matrix_scaled <- scale(current_matrix)
        if (any(is.na(current_matrix_scaled)) || any(is.infinite(current_matrix_scaled))) {
          next
        }
        # Rotate the iteration onto the master template: dimensions [274,2]
        rotation <- vegan::procrustes(X = target_matrix, Y = current_matrix_scaled, symmetric = FALSE)
        z_unpacked_aligned[i, c, , ] <- rotation$Yrot
        # Extract the 2 x 2 rotation matrix calculated for this iteration for aligning the species loadings (lambda)
        rotation_matrix <- rotation$rotation
        current_lambda_vector <- raw_lambda[i, c, ]
        current_lambda_matrix <- matrix(current_lambda_vector, nrow = mod_data$n_taxa, ncol = mod_data$n_latent)
        # multiply lambda by the same rotation matrix
        lambda_unpacked_aligned[i, c, , ] <- current_lambda_matrix %*% rotation_matrix
      }
    }
    unique_site_ids <- 1:mod_data$n_site
    first_obs_indices <- match(unique_site_ids, mod_data$site)
    all_sites_samples_raw <- z_unpacked_aligned[, , first_obs_indices, ] # dim [1200, 4, 223, 2]
    # Reshape from a 4-D array into a standard 3-D MCMC array: [1200, 4, 446]
    n_parameters <- unique_site_ids <- 1:mod_data$n_site
    first_obs_indices <- match(unique_site_ids, mod_data$site)
    all_sites_samples_raw <- z_unpacked_aligned[, , first_obs_indices, ] # dim [1200, 4, 223, 2]
    # Reshape from a 4-D array into a standard 3-D MCMC array: [1200, 4, 446]
    n_parameters <- mod_data$n_site * mod_data$n_latent
    all_sites_samples_flat <- array(all_sites_samples_raw, dim = c(n_iterations, n_chains, n_parameters))
    parameter_names <- character(n_parameters)
    idx <- 1
    for (s in unique_site_ids) {
      parameter_names[idx]     <- paste0("site_", s, "_z1")
      parameter_names[idx + 1] <- paste0("site_", s, "_z2")
      idx <- idx + 2
    }# Apply the naming dimensions to the reshaped array
    dimnames(all_sites_samples_flat) <- list(
      Iteration = 1:n_iterations,
      Chain     = 1:n_chains,
      Parameter = parameter_names
    )
    # The formal draws object
    all_sites_aligned_draws <- as_draws_array(all_sites_samples_flat)
    # 2. Recalculate draws for aligned lambda (latent factor scores for species)
    n_lambda_params <- mod_data$n_taxa * mod_data$n_latent
    lambda_flat <- array(lambda_unpacked_aligned, dim = c(n_iterations, n_chains, n_lambda_params))
    species_names <- colnames(mod_data$y)
    lambda_names <- character(n_lambda_params)
    idx <- 1
    for (j in 1:mod_data$n_taxa) {
      lambda_names[idx]     <- paste0("Species_", species_names[j], "_Load_Axis1")
      lambda_names[idx + 1] <- paste0("Species_", species_names[j], "_Load_Axis2")
      idx <- idx + 2
    }
    dimnames(lambda_flat) <- list(
      Iteration = 1:n_iterations,
      Chain     = 1:n_chains,
      Parameter = lambda_names
    )
    # Cast into a formal MCMC draws object
    lambda_aligned_draws <- as_draws_array(lambda_flat)  # dim [1200, 4, 1762] 1762 = 881 spp by 2 latent factors
    mod_summary <- mod_summary[-grep("z_exp|lambda",mod_summary$variable),]
    lambda_diagnostic_summary <- summarise_draws(lambda_aligned_draws)
    z_diagnostic_summary <- summarise_draws(all_sites_aligned_draws)
    mod_summary <- rbind(mod_summary, lambda_diagnostic_summary,z_diagnostic_summary)
    mod_diagnostics <- model_fit$diagnostic_summary()
    # mod_summary$variable[which(mod_summary$ess_bulk == min(mod_summary$ess_bulk, na.rm = TRUE))]
    diag_summary <- data.frame(spec = c("No. sampling iterations","No. divergences",
                                        "No. max_treedepth exceedences",
                                        "BFMI (Chain 1)","BFMI (Chain 2)","BFMI (Chain 3)","BFMI (Chain 4)",
                                        "Min. Bulk ESS", "Min. Tail ESS","Max R-hat",
                                        "Worst sampled parameter"),
                               value = c(n_iterations,
                                         sum(mod_diagnostics$num_divergent),
                                         sum(mod_diagnostics$num_max_treedepth),
                                         round(mod_diagnostics$ebfmi,2),
                                         round(min(mod_summary$ess_bulk, na.rm = TRUE)),
                                         round(min(mod_summary$ess_tail, na.rm = TRUE)),
                                         round(max(mod_summary$rhat, na.rm = TRUE),2),
                                         mod_summary$variable[which(mod_summary$ess_bulk == min(mod_summary$ess_bulk, na.rm = TRUE))])
    )
    if(grepl("z|axis", tolower(diag_summary$value[diag_summary$spec == "Worst sampled parameter"]))){
      non_latent <- mod_summary[-grep("z|axis", tolower(mod_summary$variable)),]
      w1 <- non_latent[order(non_latent$ess_bulk),][1,]
      diag_summary <- rbind(diag_summary, data.frame(spec = c("Min. Bulk ESS (non-latent params)",
                                                              "Worst non-latent param"),
                                                     value = c(round(w1$ess_bulk), w1$variable)))
    }
    mod_bundle$summary <- mod_summary
    mod_bundle$diag_summary <- diag_summary
    saveRDS(mod_bundle, file = paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  }
  # rm(mod_bundle,mod_summary, model_fit, mod_loo, log_lik_matrix, r_eff); gc()
}

