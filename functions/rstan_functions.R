### Some useful functions that help with rstan model diagnostics in SPARTAN:

# Calc EBFMI:
check_energy <- function(stanfit){
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  EBFMI <- rep(0, times = length(sampler_params))
  for (n in 1:length(sampler_params)) {
    energies <- sampler_params[n][[1]][,'energy__']
    numer <- sum(diff(energies)**2) / length(energies)
    denom <- var(energies)
    EBFMI[n] <- numer / denom
  }
  return(EBFMI)
}

# Number of divergences:
num_divergent <- function(stanfit){
  sampler_params <- get_sampler_params(stanfit, inc_warmup=FALSE)
  num_divergent_per_chain <- sapply(sampler_params, function(x) sum(x[, "divergent__"]))
  num_divergent_all_chains <- sum(num_divergent_per_chain)
  return(num_divergent_all_chains)
}

# Generate diagnostic report:
convergence_report <- function(stanfit, mon){
  if (!exists("mon")){
    mon <- as.data.frame(monitor(stanfit, print = FALSE))
  }
  energy <- check_energy(stanfit)
  divergences <- sum(num_divergent(stanfit))
  energy <- check_energy(stanfit)
  report <- list(
    max_Rhat = max(mon$Rhat),
    min_Bulk_ESS = min(mon$Bulk_ESS),
    min_Tail_ESS = min(mon$Tail_ESS),
    num_divergent = divergences,
    energy = round(energy, 3),
    EBFMI_min = min(energy),
    EBFMI_mean = mean(energy),
    EBFMI_max = max(energy)
  )
  return(report)
}

# Print diagnostic report:
print_convergence_report <- function(report){
  cat("Max Rhat:", report$max_Rhat, "\n")
  cat("Min Bulk ESS:", report$min_Bulk_ESS, "\n")
  cat("Min Tail ESS:", report$min_Tail_ESS, "\n")
  cat("Number of Divergences:", report$num_divergent, "\n")
  cat("EBFMI:", report$energy[1], report$energy[2], report$energy[3], report$energy[4], "\n")
}
