mwstr <- DBI::dbConnect(RPostgres::Postgres(),"mwstr")
db_curr_data <- DBI::dbGetQuery(mwstr, 
                                paste0("SELECT subcs.reach, subcs.site, ei_l44_sw1942_2022, ",
                                       "af_l5_w633_2022, meanq_mm, carea_km2, meant_y30_2022, dor_2010, ",
                                       "c_basalt, elev_asl_m FROM subcs JOIN cat_env ON subcs.site = ",
                                       "cat_env.site JOIN subc_env ON subcs.site = subc_env.site;"))
sampleable <- DBI::dbGetQuery(mwstr, "SELECT DISTINCT(site) FROM streams WHERE sampleable = 1;")
streams_to_plot <- sf::st_read(mwstr, query = paste0("SELECT streams.site, streams.reach, ",
                                                     "type, streams.geom FROM streams JOIN subcs ", 
                                                     "ON streams.site = subcs.site WHERE ", 
                                                     "sampleable = 1 OR (type IN ('pipe', ",
                                                     "'connecting line through waterbody') AND carea_km2 > 10);"))
coast <- sf::st_read(mwstr, "coast")
db_curr_data <- db_curr_data[db_curr_data$site %in% sampleable$site,]
db_ref_data <- db_curr_data
db_ref_data$ei_l44_sw1942_2022 <- 0
db_ref_data$dor_2010 <- 0
db_1750_af <- DBI::dbGetQuery(mwstr, "SELECT site, af_l5_w633_1750 FROM cat_env;")
db_ref_data$af_l5_w633_2022 <- db_1750_af$af_l5_w633_1750[match(db_ref_data$site, db_1750_af$site)]

# scaling parameters (from MW_metabarcoding_data_assembly.qmd)
load("~/uomShare/wergStaff/ChrisW/temp/model_allvars_data.rda")
scale_params <- data.frame(pred = c("I", "F", "Q", "Q2", "C", "T", "D", "B", 
                                    "L","F_I","Q_F","Q_I","B_I","Q_D"),
                           scale = NA, center = NA)

mwbugs <- DBI::dbConnect(RPostgres::Postgres(),"mwbugs")
taxon_all <- DBI::dbReadTable(mwbugs, "taxon_all")
species <- data.frame(species_no = 1:33, 
                      taxon_code = names(as.data.frame(biota_ct)),
                      species = taxon_all$taxon[match(names(as.data.frame(biota_ct)), taxon_all$shortcode)])


for(i in 1:nrow(scale_params)){
  scale_params$scale[i] <- ifelse(scale_params$pred[i] == "L", NA,
                                  attr(get(scale_params$pred[i]), "scaled:scale"))
  scale_params$center[i] <- ifelse(scale_params$pred[i] == "L", NA,
                                   attr(get(scale_params$pred[i]), "scaled:center"))
}

trnsfm_scale_preds <- function(raw_pred_data){
  if(sum(c("reach","ei_l44_sw1942_2022", "af_l5_w633_2022", "meanq_mm", "carea_km2", 
           "meant_y30_2022", "dor_2010", "c_basalt", "elev_asl_m") %in% 
         names(raw_pred_data)) != 9){
    stop(paste0("raw_pred_data must contain fields with names reach, ei_l44_sw1942_2022 ",
                "af_l5_w633_2022, meanq_mm, carea_km2, meant_y30_2022, and dor_2010 - most ",
                "easily by extracting from the mwstr database."))
    }
  x <- raw_pred_data
  x$lowness <- 1
  x$lowness[samples$elev < 30] <- -0.5
  x$lowness[samples$elev < 15] <- 0
  x$lowness[samples$elev < 7.5] <- 0.5
  x$lowness[samples$elev < 2.5] <- 1
  out_unscaled <- data.frame(I = log(x$ei_l44_sw1942_2022 + 0.001),
                             F = x$af_l5_w633_2022, 
                             Q = sqrt(x$meanq_mm),
                             Q2 = x$meanq_mm,
                             C = log(x$carea_km2),
                             T = x$meant_y30_2022,
                             D = pmin(x$dor_2010,1),
                             B = x$c_basalt,
                             L = x$lowness,
                             F_I = log(x$ei_l44_sw1942_2022 + 0.001) * x$af_l5_w633_2022,
                             Q_F = x$af_l5_w633_2022 * sqrt(x$meanq_mm), 
                             Q_I = log(x$ei_l44_sw1942_2022 + 0.001) * sqrt(x$meanq_mm),
                             B_I = log(x$ei_l44_sw1942_2022 + 0.001) * x$c_basalt,
                             Q_D = log(sqrt(x$meanq_mm) * pmin(x$dor_2010,1) + 0.1))
  out_scaled <- out_unscaled
  for(i in 1:ncol(out_scaled)){
    if(names(out_scaled)[i] != "L"){
      cent_i <- scale_params$center[scale_params$pred == names(out_scaled)[i]]
      scale_i <- scale_params$scale[scale_params$pred == names(out_scaled)[i]]
      out_scaled[,i] <- scale(out_unscaled[,i], 
                              center = cent_i,
                              scale = scale_i)
    }
  }
  list(transformed_pred_data = out_unscaled,
       scaled_pred_data = out_scaled)
}

curr_data <- trnsfm_scale_preds(db_curr_data)
ref_data <- trnsfm_scale_preds(db_ref_data)

ft_to_word <- function(ft, pgwidth = 7){
  # Set as autofit to make width parameters adjustable
  ft_out <- flextable::autofit(ft)
  # Set width as function of page width
  ft_out <- flextable::width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable::flextable_dim(ft_out)$widths))
  return(ft_out)
}

load_model_component_bundle <- function(mod_code, 
                                        mod_dir = paste0("~/uomShare/wergStaff/ChrisW/git-data/",
                                                         "mw_metabarcoding_model/trial_33spp_models/")){
  mod_dir <-paste0(mod_dir, mod_code)
  mod_files <- paste0(mod_dir, "/", dir(mod_dir))
  mod_files <- mod_files[grepl(".csv", mod_files)]
  mod_bundle <- readRDS(paste0(mod_dir, "/", mod_code, "model_bundle.rds"))
  mod_bundle
}

compile_model_diagnostics <- function(mod_code, 
                                      mod_dir = paste0("~/uomShare/wergStaff/ChrisW/git-data/",
                                                       "mw_metabarcoding_model/trial_33spp_models/")){
  mod_path <-paste0(mod_dir, mod_code)
  mod_files <- paste0(mod_path, "/", dir(mod_path))
  mod_files <- mod_files[grepl(".csv", mod_files)]
  mod_bundle <- load_model_component_bundle(mod_code, mod_dir)
  diagnostic_summary <- mod_bundle$diagnostics
  out_tab <- data.frame(spec = NA, value = NA)[0,]
  ## Metadata
  metadata <- suppressWarnings(data.table::fread(
    cmd = paste("head -n 50", shQuote(path.expand(mod_files[1])))))
  out_tab <- rbind(out_tab, data.frame(spec = metadata$V2, value = metadata$V4))
  # Number of iterations
  cmd <- "grep"; args <- c("-vc", "^#", mod_files[1])
  output_string <- system2(cmd, args, stdout = TRUE)
  out_tab <- rbind(out_tab, 
                   data.frame(spec = c("No. iterations per chain", "No. chains"),
                              value = c(round(as.numeric(output_string),-2),
                                        length(mod_files))))
  ## Time taken
  t_taken <- vector("numeric")
  for(i in 1:length(mod_files)){
    time_taken_i <- suppressWarnings(data.table::fread(
      cmd = paste("tail -n 15", shQuote(path.expand(mod_files[i]))),
      skip = "#" ))
    time_taken_i <- time_taken_i[!is.na(time_taken_i$Elapsed),]
    t_taken <- c(t_taken, as.numeric(tail(time_taken_i$Elapsed,1))/60) # time taken in minutes
  }
  t_tab <- data.frame(spec = "Modelling time (min)", value = round(max(t_taken)))
  out_tab <- rbind(out_tab, t_tab)
  # stan diagnostics
  s_tab <- data.frame(spec = c("No. divergences", "No. exceeding max_treedepth",
                               "BFMI (chain 1)","BFMI (chain 1)",
                               "BFMI (chain 3)","BFMI (chain 4)"),
                      value = c(sum(diagnostic_summary$num_divergent),
                                sum(diagnostic_summary$num_max_treedepth),
                                round(diagnostic_summary$ebfmi,2)))
  out_tab <- rbind(out_tab, s_tab)
  # ESS
  if("summary" %in% names(mod_bundle)){
    mod_summary <- mod_bundle$summary
    s_tab1 <- data.frame(spec = c("Min. Bulk ESS","Min. Tail ESS",
                                  "Least well-sampled parameter"),
                         value = c(round(min(mod_summary$ess_bulk, na.rm = TRUE)),
                                   round(min(mod_summary$ess_tail, na.rm = TRUE)),
                                   mod_summary$variable[
                                     which(mod_summary$ess_bulk == min(mod_summary$ess_bulk, na.rm = TRUE))]))
    out_tab <- rbind(out_tab, s_tab1)
  }
  out_tab
}

predict_to_new <- function(u, u_new, nriff = 1, season = 1){
n_preds <- ncol(u) # Must match the number of columns in your original 'u'
# predict to reaches assuming a riffle is present, and a spring sample
u_new$nriff <- nriff
u_new$season_mod <- season
u_new <- as.matrix(u_new[match(colnames(u),colnames(u_new))])
n_new_sites <- nrow(u_new)
# Extract the hyperparameters (standard deviations) to simulate new site/error effects
sigma_site  <- fit$draws(variables = "sigma_site",format = "draws_matrix") 
sigma_epsilon <- fit$draws(variables = "sigma_epsilon",format = "draws_matrix") 
a_taxon_draws <-  fit$draws(variables = "a_taxon",format = "draws_matrix") 
n_draws <- nrow(spno_beta_draws)
mu_new_marginal    <- matrix(NA, nrow = n_draws, ncol = n_new_sites)
prob_new_marginal  <- matrix(NA, nrow = n_draws, ncol = n_new_sites)
# binary_simulations <- matrix(NA, nrow = n_draws, ncol = n_new_sites)
system.time({
  for (i in 1:n_draws) {
    # Extract single iteration parameters
    a_tax  <- a_taxon_draws[i]
    b_tax  <- t(spno_beta_draws[i, ])
    s_site <- sigma_site[i]
    s_eps  <- sigma_epsilon[i]
    # Calculate the fixed structural effect: a_taxon + (u_new * beta)
    fixed_effect <- a_tax + as.vector(u_new %*% b_tax)
    # Scenario A: Marginal (Average Expected Probability)
    # Exclude site, spatial, and overdispersion terms by leaving them at 0
    mu_new_marginal[i, ]   <- fixed_effect
    prob_new_marginal[i, ] <- 1 / (1 + exp(-fixed_effect)) # Inverse logit link
    # # Scenario B: Conditional (Simulating Out-of-Sample Realizations)
    # # Generate new random effects for unmeasured sites
    # a_site_sim  <- rnorm(n_new_sites, mean = 0, sd = s_site)
    # epsilon_sim <- rnorm(n_new_sites, mean = 0, sd = s_eps)
    # mu_conditional <- fixed_effect + a_site_sim + epsilon_sim
    # prob_conditional <- 1 / (1 + exp(-mu_conditional))
    # # Generate actual 0/1 presence/absence data
    # binary_simulations[i, ] <- rbinom(n_new_sites, size = 1, prob = prob_conditional)
  }
})  # 24 s for 1 species
prob_new_marginal
}

# # Marginal trend plots for significant effects.
# preds[spno_beta_quantiles$X97.5. < 0 | spno_beta_quantiles$X2.5. > 0]
# # example code for unscaling (to transformed original values)
# # I * scales$scale[scales$pred == 'I'] + scales$center[scales$pred == 'I']
# # Q * scales$scale[scales$pred == 'Q'] + scales$center[scales$pred == 'Q']
# # Q2 * scales$scale[scales$pred == 'Q2'] + scales$center[scales$pred == 'Q2']
# # B * scales$scale[scales$pred == 'B'] + scales$center[scales$pred == 'B']
# # T * scales$scale[scales$pred == 'T'] + scales$center[scales$pred == 'T']
# 
# # and to re-scale
# # (log(0.1 + 0.001) - scales$center[scales$pred == 'I'])/scales$scale[scales$pred == 'I']
# 
# # mean conditions for marginal plots
# I_m <- c(min(I), (log(0.101) -scales$center[scales$pred == 'F_I'])/
#            scales$scale[scales$pred == 'F_I'] )  # 0 and 0.1 EI
# F_m <- range(F) 
# Q_m <- median(Q)  
# # (range(Q) * scales$scale[scales$pred == 'Q']) + scales$center[scales$pred == 'Q']
# # = c(3.46566, 30.19750)
# Q2_m <- median(Q2) # in both cases median is equivalent to 264.26 mm (mean is 323, which is too high)
# I_Qm <- (sqrt(264) * log(c(0, 0.1) + 0.001) -  
#            scales$center[scales$pred == 'Q_I'])/scales$scale[scales$pred == 'Q_I']
# Q_Fm <- (sqrt(264) * c(0,1) -  
#            scales$center[scales$pred == 'Q_F'])/scales$scale[scales$pred == 'Q_F']
# # if EI = 0
# I_Fm <- (log(0.001) * c(0,1) -  
#            scales$center[scales$pred == 'F_I'])/scales$scale[scales$pred == 'F_I']
# # # if EI = 0.1
# # I_F <- (log(0.101) * c(0,1) -  
# #           scales$center[scales$pred == 'F_I'])/scales$scale[scales$pred == 'F_I']
# C_m <- mean(C)
# # exp(mean(C) * scales$scale[scales$pred == 'C'] + scales$center[scales$pred == 'C']) 
# #23 km2 - fair enough 
# T_m <- mean(T) # 13.6 # quantile(T, 0.2);  
# D_m <- median(D) # 0.03 (mean is 0.144, which is too high)
# B_m <- median(B) # 0, so eastern streams (median is 0.15, which is not representative of many sites at all)


# predicted_cond_probs <- data.frame(
#   Site        = 1:n_new_sites,
#   Mean_Prob   = colMeans(binary_simulations),
#   Lower_95_CI = apply(binary_simulations, 2, quantile, probs = 0.025),
#   Upper_95_CI = apply(binary_simulations, 2, quantile, probs = 0.975)
# )
# 
# plot(seq(3.5,30,length = 20)^2, predicted_probabilities$Mean_Prob[1:20], 
#      type = "l", ylim = c(0,1), axes = FALSE,
#      xlab = "Mean annual runoff (mm)", ylab = "Probability of occurrence")
# lines(seq(3.5,30,length = 20)^2, predicted_probabilities$Lower_95_CI[1:20], lty = 3)
# lines(seq(3.5,30,length = 20)^2, predicted_probabilities$Upper_95_CI[1:20], lty = 3)
# 
# lines(seq(3.5,30,length = 20)^2, predicted_probabilities$Mean_Prob[21:40], col = "red")
# lines(seq(3.5,30,length = 20)^2, predicted_probabilities$Lower_95_CI[21:40], lty = 3, col = "red")
# lines(seq(3.5,30,length = 20)^2, predicted_probabilities$Upper_95_CI[21:40], lty = 3, col = "red")
# axis(1); axis(2, las = 1); box(bty = "l")
# # Red: EI = 0.1, Black = EI = 0 
# # Highly uncertain prob of occurrence at low Q, low EI sites.
# lines(seq(3.5,30,length = 20)^2,colMeans(binary_simulations[,1:20]), col = "black", lwd = 2)
# lines(seq(3.5,30,length = 20)^2,colMeans(binary_simulations[,21:40]), col = "red", lwd = 2)

