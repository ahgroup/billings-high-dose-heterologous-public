###
# Fitting brms models to data
# Zane Billings
# 2024-01-29
# This script fits the finalized brms models to the data. There are four
# models, each with the same structure. The only difference is the outcome and
# the outcome distribution.
# titerincrease and prevactiter both use Gaussian likelihoods
# seroconversion and seroprotection both use Binomial likelihoods
# In contrast to the previous modeling approach, all strains are fitted within
# the same model, and we allow the random effects to handle differences across
# the strains.
# Models rerun on: 2024-04-10
###

# Setup ========================================================================
## Declare dependencies ####
library(readr, include.only = NULL)
library(here, include.only = NULL)
library(cmdstanr)
library(brms)

## Load data ####
dat <- readr::read_rds(here::here("data", "processed", "model-data.Rds"))

## Load settings ####
# This script sets up variables for brms settings that are reused, so they can
# all be adjusted in the same place.
source(here::here("code", "common-functions", "brms-settings.R"))

# This script sets up convenient project paths
source(here::here("code", "common-functions", "path-setting.R"))

# This is for chunk saving the models
source(here::here("code", "common-functions", "model-manip.R"))

# Fix mixed effects models =====================================================
# In this section, we'll fit the four different mixed effects models. Each
# model has a different outcome, but they otherwise have the same functional
# form (wrt the independent variables) and the same priors, where applicable.
fit_me_models <- function(dat) {
	## titer increase outcome model ================================================
	message("ME 1/4: titer increase outcome")
	mod_titer_increase <- brms::brm(
		titer_increase ~ dose +
			s(birth_year_c, k = 5) + s(age_c, k = 5) +
			s(log_pretiter, k = 5) + s(year_c, k = 5, by = study) +
			(1 | id) + (1 | study) +
			(1 + dose | strain_type) + (1 + dose | strain_type:strain_name) +
			(1 + dose | vaccine_name) + (1 + dose | vaccine_name:strain_type),
		data = dat,
		family = gaussian(),
		prior = priors_me_gaussian_outcome,
		backend = "cmdstanr",
		silent = 0,
		refresh = BRMS_ITERS / 5,
		iter = BRMS_ITERS,
		warmup = BRMS_WARMUPS,
		chains = BRMS_CHAINS,
		cores = BRMS_CORES,
		seed = BRMS_SEED,
		control = list(
			adapt_delta = BRMS_ADAPT_DELTA,
			max_depth = BRMS_MAX_TREEDEPTH
		)
	)
	
	readr::write_rds(
		mod_titer_increase,
		here::here(largefile_path, "fits", "mod-ti.Rds"),
		compress = "gz"
	)
	
	## Post-titer outcome model ====================================================
	message("ME 2/4: posttiter outcome")
	mod_posttiter <- brms::brm(
		log_posttiter ~ dose +
			s(birth_year_c, k = 5) + s(age_c, k = 5) +
			s(log_pretiter, k = 5) + s(year_c, k = 5, by = study) +
			(1 | id) + (1 | study) +
			(1 + dose | strain_type) + (1 + dose | strain_type:strain_name) +
			(1 + dose | vaccine_name) + (1 + dose | vaccine_name:strain_type),
		data = dat,
		family = gaussian(),
		prior = priors_me_gaussian_outcome,
		backend = "cmdstanr",
		silent = 0,
		refresh = BRMS_ITERS / 5,
		iter = BRMS_ITERS,
		warmup = BRMS_WARMUPS,
		chains = BRMS_CHAINS,
		cores = BRMS_CORES,
		seed = BRMS_SEED,
		control = list(
			adapt_delta = BRMS_ADAPT_DELTA,
			max_depth = BRMS_MAX_TREEDEPTH
		)
	)
	
	readr::write_rds(
		mod_posttiter,
		here::here(largefile_path, "fits", "mod-pt.Rds"),
		compress = "gz"
	)
	
	## Seroconversion outcome model ================================================
	message("ME 3/4: seroconversion outcome")
	mod_seroconversion <- brms::brm(
		seroconversion ~ dose +
			s(birth_year_c, k = 5) + s(age_c, k = 5) +
			s(log_pretiter, k = 5) + s(year_c, k = 5, by = study) +
			(1 | id) + (1 | study) +
			(1 + dose | strain_type) + (1 + dose | strain_type:strain_name) +
			(1 + dose | vaccine_name) + (1 + dose | vaccine_name:strain_type),
		data = dat,
		family = bernoulli("logit"),
		prior = priors_me_binomial_outcome,
		backend = "cmdstanr",
		silent = 0,
		refresh = BRMS_ITERS / 5,
		iter = BRMS_ITERS,
		warmup = BRMS_WARMUPS,
		chains = BRMS_CHAINS,
		cores = BRMS_CORES,
		seed = BRMS_SEED,
		control = list(
			adapt_delta = BRMS_ADAPT_DELTA,
			max_depth = BRMS_MAX_TREEDEPTH
		)
	)
	
	readr::write_rds(
		mod_seroconversion,
		here::here(largefile_path, "fits", "mod-sc.Rds"),
		compress = "gz"
	)
	
	## Seroprotection outcome model ================================================
	message("ME 4/4: seroprotection outcome")
	mod_seroprotection <- brms::brm(
		seroprotection ~ dose +
			s(birth_year_c, k = 5) + s(age_c, k = 5) +
			s(log_pretiter, k = 5) + s(year_c, k = 5, by = study) +
			(1 | id) + (1 | study) +
			(1 + dose | strain_type) + (1 + dose | strain_type:strain_name) +
			(1 + dose | vaccine_name) + (1 + dose | vaccine_name:strain_type),
		data = dat,
		family = bernoulli("logit"),
		prior = priors_me_binomial_outcome,
		backend = "cmdstanr",
		silent = 0,
		refresh = BRMS_ITERS / 5,
		iter = BRMS_ITERS,
		warmup = BRMS_WARMUPS,
		chains = BRMS_CHAINS,
		cores = BRMS_CORES,
		seed = BRMS_SEED,
		control = list(
			adapt_delta = BRMS_ADAPT_DELTA,
			max_depth = BRMS_MAX_TREEDEPTH
		)
	)
	
	readr::write_rds(
		mod_seroprotection,
		here::here(largefile_path, "fits", "mod-sp.Rds"),
		compress = "gz"
	)
	
	invisible(TRUE)
}

# Run model fits ===============================================================
main <- function(ME = TRUE) {
	# Fit the models
	message("Starting to fit mixed effects models.")
	fit_me_models(dat)
	message("Done fitting.")
	
	# Split up the models for saving
	message("Saving models now.")
	model_split()
	
	invisible(TRUE)
}

# TODO fix parallel processing setup to dynamically adjust # of cores passed
# to brms
main(
	ME = TRUE
)

# END OF FILE ==================================================================
