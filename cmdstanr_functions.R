suppressPackageStartupMessages(library(cmdstanr))
suppressPackageStartupMessages(library(posterior))
suppressPackageStartupMessages(library(data.table))

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
    s_tab1 <- data.frame(spec = c("Min. Bulk ESS","Min. Tail ESS","Max Rhat",
                                  "Least well-sampled parameter"),
                         value = c(round(min(mod_summary$ess_bulk, na.rm = TRUE)),
                                   round(min(mod_summary$ess_tail, na.rm = TRUE)),
                                   round(max(mod_summary$rhat, na.rm = TRUE),2),
                                   mod_summary$variable[
                                     which(mod_summary$ess_bulk == min(mod_summary$ess_bulk, na.rm = TRUE))]))
    out_tab <- rbind(out_tab, s_tab1)
  }
  out_tab
}

base_traceplot <- function(param_draws, param_name = "Unspecified parameter"){
  colpal <- RColorBrewer::brewer.pal(4,"PuOr")
  plot(param_draws[,1,], type = "l", col = colpal[1], axes = FALSE,
       ylim = range(param_draws,na.rm=TRUE), xlab = "Iteration", 
       ylab = param_name)
  for(i in 2:dim(param_draws)[2]){
    lines(param_draws[,i,], col = colpal[i])
  }
  axis(1);axis(2, las = 1); box(bty = "l")
}
