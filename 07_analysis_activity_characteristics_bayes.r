# packages
library(data.table)
library(dplyr)
library(glmmTMB)
library(rjags)
library(coda)

# path to your project (adjust if needed)
path<- "/Users/pandadamda/rts_activity/"
path<- "E:/Uni_Arbeit/rts_activity/"
path<- "D:/rts_activity/"

## load df
df<- fread(file.path(path, "data/bird_data_storage/activity_characteristics/activity_characteristics_individual.csv"), stringsAsFactors = TRUE)
df$date_f <- as.factor(as.character(df$date_f))


variables<- c("length_act",
              "act_at_sunset_rel",
              "act_at_sunrise_rel",
              "auc_sunrise",
              "steepest_ascend",
              "steepest_descend",
              "auc" )


results<- data.frame(variable=character(0),
                     level=character(0),
                     VC_lmm=numeric(0),
                     VCprop_lmm=numeric(0),
                     VC_jags=numeric(0),
                     VCprop_jags=numeric(0),
                     CI_low=numeric(0),
                     CI_high=numeric(0))

# loop through all variables of interest
for(i in 1:length(variables)){
  
  ## create empty result table
  results_i<- data.frame(variable=rep(variables[i], 3),
                         level=c("interspecific", "intraspecific", "intraindividual"),
                         VC_lmm=rep(NA,3),
                         VCprop_lmm=rep(NA,3),
                         VC_jags=rep(NA,3),
                         VCprop_jags=rep(NA,3),
                         CI_low=rep(NA,3),
                         CI_high=rep(NA,3))
  
  ##############################################################################
  # Variance Component Analysis LMM style ####
  formula_string <- paste("scale(",variables[i], ")~ 1 + (1|species_en) + (1|species_en:ring_ID)")
  lmm_formula <- as.formula(formula_string)
  lmm <- glmmTMB(lmm_formula, data = df)
  summary(lmm)

# extract SDs (for plausible initial values)
sd_species_est <- sqrt(VarCorr(lmm)$cond$species_en[1,1])
sd_ind_est     <- sqrt(VarCorr(lmm)$cond$ring_ID[1,1])
sd_resid_est   <- sigma(lmm)# residual SD


results_i$VC_lmm[results_i$level == "interspecific"] <- VarCorr(lmm)$cond$species_en[1,1]
results_i$VC_lmm[results_i$level == "intraspecific"] <- VarCorr(lmm)$cond$ring_ID[1,1]
results_i$VC_lmm[results_i$level == "intraindividual"]  <-sigma(lmm)^2

total<- results_i$VC_lmm[results_i$level == "interspecific"]+
         results_i$VC_lmm[results_i$level == "intraspecific"]+
  results_i$VC_lmm[results_i$level == "intraindividual"]

results_i$VCprop_lmm[results_i$level == "interspecific"] <- results_i$VC_lmm[results_i$level == "interspecific"] /total
results_i$VCprop_lmm[results_i$level == "intraspecific"]<-results_i$VC_lmm[results_i$level == "intraspecific"] /total
results_i$VCprop_lmm[results_i$level == "intraindividual"] <- results_i$VC_lmm[results_i$level == "intraindividual"]/ total


# ---------------------------
# Prepare data for JAGS
# ---------------------------
# Create indices:
df_j <- df %>%
  filter(!is.na(!!sym(variables[i]))) %>% # Remove NAs for the current variable
  mutate(
    species_idx = as.integer(factor(species_en)),
    ind_global_char = paste(species_en, ring_ID, sep = "_"),
    ind_idx = as.integer(factor(ind_global_char))
  ) %>%
  arrange(species_idx, ind_idx)

# Observations (and ensure it is a vector)
y <- as.vector(scale(df_j[[variables[i]]]))
N <- nrow(df_j)
species_idx <- df_j$species_idx
ind_idx <- df_j$ind_idx
N_species <- length(unique(species_idx))

N_individuals <- length(unique(ind_idx))

# data list for JAGS
jags_data <- list(
  N = N,
  N_species = N_species,
  N_individuals = N_individuals,
  y = y,
  species_idx = species_idx,
  ind_idx = ind_idx
)

# ---------------------------
# Write JAGS model file (non-centered additive random effects)
# ---------------------------
model_file <- file.path(path, "scripts/model_jags.r")


# ---------------------------
# Initial values for chains
# ---------------------------
set.seed(123)
n_chains <- 4
inits_list <- list()
for (ch in 1:n_chains) {
  inits_list[[ch]] <- list(
    mu = mean(y, na.rm = TRUE) + rnorm(1, 0, 0.1),
    sigma_species = abs(sd_species_est * runif(1, 0.5, 1.5)) + 1e-3,
    sigma_ind = abs(sd_ind_est * runif(1, 0.5, 1.5)) + 1e-3,
    sigma_error = abs(sd_resid_est * runif(1, 0.8, 1.2)) + 1e-3,
    z_species = rnorm(N_species, 0, 0.1),
    z_ind = rnorm(N_individuals, 0, 0.1)
  )
}

# ---------------------------
# Run JAGS
# ---------------------------
jm <- jags.model(file = model_file,
                 data = jags_data,
                 inits = inits_list,
                 n.chains = n_chains,
                 n.adapt = 2000)

# burn-in
update(jm, n.iter = 20000)

# parameters to monitor
params <- c("mu", "sigma_species", "sigma_ind", "sigma_error",
            "var_species", "var_ind", "var_error")

# draw samples
samples_coda <- coda.samples(jm, variable.names = params, n.iter = 8000, thin = 5)

# ---------------------------
# Diagnostics & summaries
# ---------------------------
# trace plots and densities
plot(samples_coda)

# Gelman-Rubin (Rhat)
print(gelman.diag(samples_coda, multivariate = FALSE))

# effective sample size
print(effectiveSize(samples_coda))

# Posterior summaries
m <- as.matrix(samples_coda)

# Extract posterior draws for variances
var_species_post <- m[ , "var_species"]
var_ind_post     <- m[ , "var_ind"]
var_error_post   <- m[ , "var_error"]

# compute proportions per draw (preserve uncertainty)
total_post <- var_species_post + var_ind_post + var_error_post
prop_species_post <- var_species_post / total_post
prop_ind_post     <- var_ind_post / total_post
prop_error_post   <- var_error_post / total_post

# Summaries (median and 95% credible intervals)
# Calculate variance components (posterior means)
results_i$VC_jags[results_i$level == "interspecific"] <- median(var_species_post)
results_i$VC_jags[results_i$level == "intraspecific"] <- median(var_ind_post)
results_i$VC_jags[results_i$level == "intraindividual"]  <- median(var_error_post)

results_i$CI_low[results_i$level == "interspecific"]<- quantile(var_species_post, 0.025)
results_i$CI_high[results_i$level == "interspecific"] <- quantile(var_species_post, 0.975)
results_i$CI_low[results_i$level == "intraspecific"]<- quantile(var_ind_post, 0.025)
results_i$CI_high[results_i$level == "intraspecific"] <- quantile(var_ind_post, 0.975)
results_i$CI_low[results_i$level == "intraindividual"] <- quantile(var_error_post, 0.025)
results_i$CI_high[results_i$level == "intraindividual"]  <- quantile(var_error_post, 0.975)

results_i$VCprop_jags[results_i$level == "interspecific"] <- median(prop_species_post)
results_i$VCprop_jags[results_i$level == "intraspecific"]<- median(prop_ind_post)
results_i$VCprop_jags[results_i$level == "intraindividual"] <- median(prop_error_post)

results<- rbind(results, results_i)
}


results$percent_lmm <- results$VCprop_lmm * 100
results$percent_jags <- results$VCprop_jags * 100

fwrite(results, paste0(path, "output/model_output/VCA/jags_results.csv"))
openxlsx::write.xlsx(results, paste0(path, "output/model_output/VCA/jags_results.xlsx"))

results[3:8]<- round(results[3:8], digits = 2)

library(ggplot2)
ggplot(results,aes(y=variable, x=percent_jags, fill=level, label = round(percent_jags,digits=2))) +
         geom_bar(stat="identity")+
  scale_fill_brewer(palette="Paired")+
  #geom_text(size = 5, position = position_stack(vjust = 1.025))+
  xlab("Variance Partition Coefficients [%]")+
  ylab("Activity characteristic")+
  theme_minimal()+
  theme(text = element_text(size=25),
        axis.text = element_text(face="bold"),
        legend.title=element_blank(),
        legend.position = c(0.84, 0.88),
        legend.text = element_text(size=25),
        legend.box.background = element_rect(colour = "black"))
ggsave(filename = paste0(path, "output/model_output/VCA/jags_results.png") , width = 15, height = 9)

