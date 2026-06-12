library(cmdstanr)
library(posterior)
library(knitr)     # Required for kable() markdown formatting
library(data.table)

save_model_summary <- function(mod_code, 
                               mod_dir = "~/uomShare/wergStaff/ChrisW/git-data/=mw_metabarcoding_model/trial_33spp_models/"){
  mod_dir <-paste0(mod_dir, mod_code)
  mod_files <- paste0(mod_dir, "/", dir(mod_dir))
  mod_files <- mod_files[!grepl("diagnostic_summary", mod_files)]
  cmdstan_path <- cmdstanr::cmdstan_path()
  stansummary_binary <- file.path(cmdstan_path, "bin", "stansummary")
  output_csv_summary <- paste0("temp_data/", mod_code, "_diagnostic_summary.csv")
  cmd <- paste0('"', stansummary_binary, '" --csv_filename="', output_csv_summary, '" ', paste(mod_files, collapse = " "))
  system(cmd)
  system(paste0("mv temp_data/", mod_code,  output_csv_summary, " ", mod_dir))
}

compile_model_diagnostics <- function(mod_code,
                                      mod_dir = "~/uomShare/wergStaff/ChrisW/git-data/=mw_metabarcoding_model/trial_33spp_models/"){
  mod_dir <-paste0(mod_dir, mod_code)
  mod_files <- paste0(mod_dir, "/", dir(mod_dir))
  mod_files <- mod_files[!grepl("diagnostic_summary", mod_files)]
  mod_summary <- paste0(mod_dir, "/", mod_code, "_diagnostic_summary.csv")
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
  # ESS
  mod_summary <- read.csv(mod_summary)
  s_tab <- data.frame(spec = c("Min. Bulk ESS","Min. Tail ESS","Least well-sampled parameter"),
                      value = c(min(mod_summary$ESS_bulk, na.rm = TRUE),
                                min(mod_summary$ESS_tail, na.rm = TRUE),
                                mod_summary$name[which(mod_summary$ESS_bulk.s == min(mod_summary$ESS_bulk.s, na.rm = TRUE))]))
  out_tab <- rbind(out_tab, s_tab)
}

