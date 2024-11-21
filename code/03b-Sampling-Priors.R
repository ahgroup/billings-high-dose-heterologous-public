###
# Sampling from priors of brms models
# Zane Billings
# 2024-10-31
# This script takes the models as described in script 3 and samples from
# the priors without exposing the models to any data.
###

# Setup ========================================================================
## Declare dependencies ####
library(readr, include.only = NULL)
library(here, include.only = NULL)
library(cli, include.only = NULL)
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
sample_me_priors <- function(dat) {
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
			max_treedepth = BRMS_MAX_TREEDEPTH
		),
		sample_prior = "only"
	)
	
	readr::write_rds(
		mod_titer_increase,
		here::here(largefile_path, "priors", "mod-ti.Rds"),
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
			max_treedepth = BRMS_MAX_TREEDEPTH
		),
		sample_prior = "only"
	)
	
	readr::write_rds(
		mod_posttiter,
		here::here(largefile_path, "priors", "mod-pt.Rds"),
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
			max_treedepth = BRMS_MAX_TREEDEPTH
		),
		sample_prior = "only"
	)
	
	readr::write_rds(
		mod_seroconversion,
		here::here(largefile_path, "priors", "mod-sc.Rds"),
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
			max_treedepth = BRMS_MAX_TREEDEPTH
		),
		sample_prior = "only"
	)
	
	readr::write_rds(
		mod_seroprotection,
		here::here(largefile_path, "priors", "mod-sp.Rds"),
		compress = "gz"
	)
	
	invisible(TRUE)
}

# Run model fits ===============================================================
main <- function(ME = TRUE) {
	dir.create(
		here::here(largefile_path, "priors"),
		showWarnings = FALSE,
		recursive = TRUE
	)
	cli::cli_alert_info(
		"Results will be stored at: {.path {here::here(largefile_path, 'priors')}}"
	)
	# Fit the models
	message("Starting to fit mixed effects models.")
	sample_me_priors(dat)
	message("Done fitting.")
	
	# Split up the models for saving
	message("Saving models now.")
	model_split(priors = TRUE)
	
	invisible(TRUE)
}

# TODO fix parallel processing setup to dynamically adjust # of cores passed
# to brms
main(
	ME = TRUE
)

# END OF FILE ==================================================================
