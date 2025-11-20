model {
  # Grand mean (prior is fine for scaled data where mean is ~0)
  mu ~ dnorm(0, 1)
  
  # --- Priors for Standard Deviations (Half-Normal) ---
  # More robust and theoretically sound than uniform priors.
  sigma_species ~ dnorm(0, 1) T(0,)
  sigma_ind     ~ dnorm(0, 1) T(0,)
  sigma_error   ~ dnorm(0, 1) T(0,)
  
  # --- Non-Centered Parameterization for ALL Random Effects ---
  
  # Species-level effects
  for (s in 1:N_species) {
    z_species[s] ~ dnorm(0, 1)
    a_species[s] <- z_species[s] * sigma_species
  }
  
  # Individual-level effects
  for (j in 1:N_individuals) {
    z_ind[j] ~ dnorm(0, 1)
    b_ind[j] <- z_ind[j] * sigma_ind
  }
  
  # --- Likelihood ---
  tau_error <- 1 / (sigma_error * sigma_error)
  for (n in 1:N) {
    mu_n[n] <- mu + a_species[ species_idx[n] ] + b_ind[ ind_idx[n] ]
    y[n] ~ dnorm(mu_n[n], tau_error)
  }
  
  # --- Derived Variances for Monitoring ---
  var_species <- sigma_species * sigma_species
  var_ind     <- sigma_ind * sigma_ind
  var_error   <- sigma_error * sigma_error
}