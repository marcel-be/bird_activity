library(dplyr)
library(data.table)
library(glmmTMB)
library(rjags)
library(coda)
library(brms)

path<- "E:/Uni_Arbeit/rts_activity/"
path<- "D:/rts_activity/"
df<- fread(file.path(path, "data/bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"), stringsAsFactors = TRUE)
df$date_f <- as.factor(as.character(df$date_f))


# simulation truths
var_species_true  <- 1
var_ind_true      <- 0.5
var_residual_true <- 0.3
mu_true           <- 10
N_simulations     <- 100

# prepare results frames
results_lmm <- data.frame(sim_id = 1:N_simulations,
                          var_species_est = NA,
                          var_ind_est = NA,
                          var_residual_est = NA)
results_brms <- results_lmm

results_jags <- results_lmm

# path to model file - update to your path
model_file <- file.path(path, "scripts/model_jags.r")

for (i in seq_len(N_simulations)) {
  
  # --------------------------
  # simulate on your real df
  # --------------------------
  df_temp <- df %>%
    mutate(species = factor(species_en),
           ring_ID = factor(ring_ID))
  
  # simulate random effects respecting existing levels
  species_eff <- data.frame(species = levels(df_temp$species),
                            species_eff = rnorm(n = nlevels(df_temp$species), 0, sqrt(var_species_true)))
  ind_eff <- data.frame(ring_ID = levels(df_temp$ring_ID),
                        ind_eff = rnorm(n = nlevels(df_temp$ring_ID), 0, sqrt(var_ind_true)))
  
  df_sim <- df_temp %>%
    left_join(species_eff, by = "species") %>%
    left_join(ind_eff, by = "ring_ID")
  
  mu_vec <- mu_true + df_sim$species_eff + df_sim$ind_eff
  df_sim$y_sim <- rnorm(nrow(df_sim), mu_vec, sqrt(var_residual_true))
  
  
  
  # --------------------------
  # frequentist REML LMM
  # --------------------------
  
  fit <- suppressMessages(glmmTMB(scale(y_sim) ~ (1|species) + (1|species_en:ring_ID), data = df_sim))
  
  # extract variance components (conditional)
  vc <- VarCorr(fit)$cond
  # if any are missing, fallback
  ring_var <- if(!is.null(vc$ring_ID)) vc$ring_ID[1,1] else 0.1
  species_var <- if(!is.null(vc$species)) vc$species[1,1] else 0.1
  resid_var <- sigma(fit)^2
  
  results_lmm$var_species_est[i] <- species_var
  results_lmm$var_ind_est[i] <- ring_var
  results_lmm$var_residual_est[i] <- resid_var
  
  
  
  # --------------------------
  # brms
  # --------------------------
  
  # Priors 
  priors <- c(
    prior(normal(0, 1), class = "Intercept"),
    prior(student_t(3, 0, 0.5), class = "sd", group = "species"),  # species (only 8 groups)
    prior(student_t(3, 0, 1), class = "sd", group = "ring_ID"),      # individual
    prior(student_t(3, 0, 1), class = "sigma")
  )
  
  # Model
  fit <- brm(
    formula = scale(y_sim) ~ 1 + (1 | species) + (1 | species:ring_ID),
    data = df_sim,
    prior = priors,
    iter = 4000, warmup = 2000, chains = 4,
    control = list(adapt_delta = 0.95, max_treedepth = 15)
  )
  
  # Diagnostics & extraction
  vc<- VarCorr(fit)            # variance components
  results_brms$var_species_est[i] <- vc$species$sd[1]^2
  results_brms$var_ind_est[i] <- vc$ring_ID$sd[1]^2
  results_brms$var_residual_est[i] <-vc$residual__$sd[1]^2

  
  
  
  
  # --------------------------
  # Jags
  # --------------------------
  df_j <- df_sim %>%
    mutate(species_idx = as.integer(factor(species)),
           ind_idx = as.integer(factor(paste0(species, "_", ring_ID)))) # unique across species
  
  y_numeric <- as.numeric(scale(df_j$y_sim))  # <--- IMPORTANT: numeric vector
  
  jags_data <- list(
    N = length(y_numeric),
    N_species = length(unique(df_j$species_idx)),
    N_individuals = length(unique(df_j$ind_idx)),
    y = y_numeric,
    species_idx = df_j$species_idx,
    ind_idx = df_j$ind_idx
  )
  
  # inits for JAGS
  n_chains <- 4
  set.seed(100 + i)
  inits_list <- vector("list", n_chains)
  sd_species_est <- sqrt(species_var)
  sd_ind_est <- sqrt(ring_var)
  sd_resid_est <- sqrt(resid_var)
  for (ch in seq_len(n_chains)) {
    inits_list[[ch]] <- list(
      mu = mean(y_numeric, na.rm = TRUE) + rnorm(1,0,0.1),
      sigma_species = abs(sd_species_est * runif(1, 0.5, 1.5)) + 1e-3,
      sigma_ind     = abs(sd_ind_est * runif(1, 0.5, 1.5)) + 1e-3,
      sigma_error   = abs(sd_resid_est * runif(1, 0.8, 1.2)) + 1e-3,
      z_ind = rnorm(jags_data$N_individuals, 0, 0.1)
    )
  }
  
  # run JAGS
  jm <- jags.model(file = model_file,
                   data = jags_data,
                   inits = inits_list,
                   n.chains = n_chains,
                   n.adapt = 2000)
  
  update(jm, n.iter = 20000)  # burn-in (increase for final analyses)
  
  params <- c("mu", "sigma_species", "sigma_ind", "sigma_error",
              "var_species", "var_ind", "var_error")
  
  samples_coda <- coda.samples(jm, variable.names = params, n.iter = 8000, thin = 5)

  # extract posterior medians and store
  m <- as.matrix(samples_coda)
  results_jags$var_species_est[i] <- median(m[ , "var_species"])
  results_jags$var_ind_est[i]     <- median(m[ , "var_ind"])
  results_jags$var_residual_est[i]<- median(m[ , "var_error"])
  
}

results_brms$method<- "brms"
results_lmm$method <- "lmm"
results_jags$method<- "jags"
results_sim<- rbind(results_brms, results_jags, results_lmm)

fwrite(results_sim, paste0(path, "output/model_output/VCA/simulation_results.csv"))




results_old<- fread(file.path(path, "output/model_output/VCA/simulation_results.csv"), stringsAsFactors = TRUE)
#results_sim<- rbind(results_old, results_sim) 



## true values: species=1, individual=0.5, residual=0.3
rescaling_factor<- sum(var_species_true + var_ind_true+ var_residual_true)
total<- results_sim %>% 
  filter(!is.na(var_species_est)) %>% 
  group_by(method) %>% 
  summarise(n= n(),
            mean_var_species_est=mean(var_species_est)*rescaling_factor,
            mean_var_ind_est = mean(var_ind_est)*rescaling_factor,
            mean_var_residual_est=mean(var_residual_est)*rescaling_factor,
            mae_species = sum(abs(mean_var_species_est - var_species_est*rescaling_factor))/n,
            mae_ind = sum(abs(mean_var_ind_est - var_ind_est*rescaling_factor))/n,
            mae_residual = sum(abs(mean_var_residual_est - var_residual_est*rescaling_factor))/n)
# Calculate how much larger the estimated variance-component values are compared to the true values
((total$mean_var_species_est-var_species_true)/var_species_true)*100
((total$mean_var_ind_est-var_ind_true)/var_ind_true)*100
((total$mean_var_residual_est-var_residual_true)/var_residual_true)*100
