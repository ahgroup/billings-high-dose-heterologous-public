###
# Posterior Summaries
# Zane Billings
# 2024-02-01
# This script takes in the fitted brms model objects, and summarizes the
# posterior distributions.
###

# Setup ========================================================================
## Package requirements ####
# We need to attach brms to prevent any funny business.
library(brms)
library(ggdist)
library(here, include.only = NULL)
library(readr, include.only = NULL)
library(rlang, include.only = NULL)
library(crayon, include.only = NULL)
library(dplyr, include.only = NULL)
library(forcats, include.only = NULL)
library(purrr, include.only = NULL)
library(stringr, include.only = NULL)
library(tibble, include.only = NULL)
library(tidyr, include.only = NULL)

## Source scripts as needed ####
# This script sets up convenient project paths
source(here::here("code", "common-functions", "path-setting.R"))

# This script recreates the model objects
source(here::here("code", "common-functions", "model-manip.R"))

# Posterior summaries ==========================================================

# Function to get the counterfactual model predictions for each person
get_contrast_preds <- function(.data, mod, raw = TRUE) {
	
	# Get counterfactual HD predictions
	hd_preds <-
		.data |>
		dplyr::mutate(dose = "HD") |>
		brms::posterior_epred(
			object = mod,
			newdata = _,
			re_formula = NULL
		)
	
	# Get counterfactual SD predictions
	sd_preds <-
		.data |>
		dplyr::mutate(dose = "SD") |>
		brms::posterior_epred(
			object = mod,
			newdata = _,
			re_formula = NULL
		)
	
	# Get the contrast by subtracting those matrices elementwise
	contrast_preds <- hd_preds - sd_preds
	
	# If the raw argument is TRUE, just return that. Otherwise, we need to tidy
	# the results.
	if (isTRUE(raw)) {
		return (contrast_preds)
	} else {
		
		# Now we just need to tidy up the results
		contrast_predictions_wide <-
			contrast_preds |>
			`colnames<-`(1:ncol(contrast_preds)) |>
			tibble::as_tibble() |>
			tibble::rowid_to_column(var = "sample") |>
			tidyr::pivot_longer(cols = -sample, names_to = "subject") |>
			tidyr::pivot_wider(names_from = sample, values_from = value)
		
		# Bind the predictions to the dataset and then pivot the samples longer.
		contrast_augmented <-
			.data |>
			# Since we already have the contrasts over the dose, specify that column
			# so that we don't do anything silly later
			dplyr::mutate(dose = "HD - SD") |>
			dplyr::bind_cols(contrast_predictions_wide) |>
			tidyr::pivot_longer(
				cols = dplyr::num_range("", 1:nrow(contrast_preds)),
				names_to = "sample",
				values_to = "epred"
			)
		
		return(contrast_augmented)
	}
}

# This function extracts the information on the total number of samples that
# are in a fitted brms object. This is so we can automatically know how
# many samples there are for counting the individuals in each stratum.
count_samples <- function(mod) {
	# Code copied from print.brmssummary and summary.brmsfit to get the
	# number of total samples for a given brms model.
	n_samples <- ceiling(
		(mod$fit@sim$iter - mod$fit@sim$warmup) /
			brms:::nthin(mod) * brms::nchains(mod)
	)
	return(n_samples)
}

# Function for summarizing posterior ITEs. If we define a function here, it is
# easier to change the way ITEs are summarized into CATEs.
summarize_post <- function(x, ...) {
	return(
		ggdist::mean_hdci(
			x,
			...,
			.width = 0.95,
			n = 4096
		)
	)
}

# Function to calculate the ATE or conditional ATE, which is defined as the
# average of the counterfactual effects when treatment = 1 minus the effects
# when treatment = 0 for each person. For the CATE, only average within strata.
cate <- function(ca, .data, ..., label = NULL, summary_fun = summarize_post) {
	# This function is memory intensive, so when we exit make sure to run
	# garbage collection. Lots of people online say this doesn't do anything, but
	# I notice a substantial improvement in vRAM usage.
	on.exit(invisible(gc()), add = TRUE)
	
	# The input is a numeric matrix of counterfactual predictions (HD - SD) for
	# each person, with dimensions n_samples by n_observations.
	# So for each cATE, we need to figure out which person indices go into each
	# group, and do the matrix operations on those. We need this part to accept
	# ... as an argument which is a list of bare names, so this part needs to be
	# tidyverse. Fortunately dplyr provides a really nice and convenient way to
	# get this.
	gdata <- .data |>
		dplyr::group_by(...) |>
		dplyr::group_data() |>
		dplyr::mutate(
			# Get the size of each group as well, unfortunately there isn't a slick
			# way to use dplyr::group_size() here
			n = lengths(.rows),
			# And add the specified label, which is recycled to repeat over all rows
			cate_label = label
		)
	
	# Next we map over that group data frame. For each row, we want to summarize
	# the posterior samples using some summary function, over only the
	# observations (which are the columns of the contrast matrix)
	# specified in the integer vector in the .rows column. Fortunately we can
	# easily do that by subsetting the contrast matrix.
	res <-
		purrr::map(
			gdata$.rows,
			\(idx) {ca[, idx] |> c() |> summary_fun()}
		) |>
		purrr::list_rbind()
	
	out <- dplyr::bind_cols(gdata, res)
	
	return(out)
}

# Function to get the specific cATEs we want for the paper -- can add more
# strata to this as necessary.
get_estimate_set <- function(ctrs, .data, data_label = NULL) {
	
	# Calculate the various CATE estimates for strata that we want
	cate_list <- list(
		# Overall effect
		cate(ctrs, .data, label = "overall"),
		# Conditional on vaccine
		cate(ctrs, .data, vaccine_name, label = "vaccine"),
		# Conditional on vaccine and strain
		cate(
			ctrs, .data, vaccine_name, strain_name,
			label = "vaccine-strain"
		),
		# Conditional on season
		cate(ctrs, .data, season, label = "season")
	)
	
	# Combine the estimates into one data frame.
	out <-
		cate_list |>
		dplyr::bind_rows() |>
		dplyr::mutate(
			# Fill in the missing factor levels with "overall"
			dplyr::across(
				c(vaccine_name, strain_name, season),
				\(x) forcats::fct_na_value_to_level(x, level = "overall")
			),
			# Add the data label as a column
			data_label = data_label
		)
	
	return(out)
}

# Function to get the CATEs on all strains, heterologous only, and homologous
# only. This calls the CATE set estimating function on all three of those sets. 
get_multiple_cates <- function(ctrs, .data, model_label = NULL) {
	
	# Get the estimates combined on all strains
	all_strains <- .data |>
		get_estimate_set(ctrs, .data = _, data_label = "all")
	
	# Get the estimates for only the homologous strains included
	homologous <- .data |>
		dplyr::filter(as.character(strain_name) == as.character(vaccine_name)) |>
		get_estimate_set(ctrs, .data = _, data_label = "homologous")
	
	# Get the estimates with the homologous strains EXCLUDED
	heterologous <- .data |>
		dplyr::filter(as.character(strain_name) != as.character(vaccine_name)) |>
		get_estimate_set(ctrs, .data = _, data_label = "heterologous")
	
	# Combine the estimates into one data frame
	out <-
		list(all_strains, homologous, heterologous) |>
		dplyr::bind_rows() |>
		dplyr::mutate(model_label = model_label)
	
	return(out)
}

# Function to get all of the CATE estimates for the homologous, heterologous,
# and combined strain sets, for multiple models.
get_cates_for_all_models <- function(.data, model_list, model_names) {
	
	# Start the nested CATE cascade -- this function calls the data set splitting
	# CATE function for a given model; that function splits the data based on
	# strain sets and calls the CATE set function; the CATE set function calls
	# the CATE estimator function for each of the strata we're interested in.
	cate_set_list <-
		purrr::map2(
			model_list, model_names,
			\(model, model_name) {
				# First calculate the HD/SD counterfactual contrasts for every observation
				# in the data set. If we do this here, we only have to do it once per model
				# because it's augmented to the actual data.
				#n_samps <- count_samples(model)
				ctrs <- get_contrast_preds(.data, model, raw = TRUE)
				
				# Now do the CATEs for that model.
				cates <- get_multiple_cates(ctrs, .data, model_label = model_name)
				
				# Save the entire posterior sample for supplementary figures
				# Commented out because we didn't end up needing them and the files
				# will be at least 2GB each.
				# readr::write_rds(
				# 	x = ctrs,
				# 	file = here::here(largefile_path, "post", paste0(model_name, ".Rds")),
				# 	compress = "gz"
				# )
				
				return(cates)
			},
			.progress = "Getting CATEs."
		)
	
	# Bind them together into one data frame
	out <- dplyr::bind_rows(cate_set_list) |>
		# Make the label variables into factors
		dplyr::mutate(
			dplyr::across(c(cate_label, data_label, model_label), as.factor)
		)
	return(out)
}

# Code execution and saving results ============================================
main <- function() {
	rlang::inform(crayon::green("Loading data and models."))
	
	# Run the model de-chunker if necessary
	mod_base <- here::here(largefile_path, "mods")
	if (length(list.files(mod_base)) != 4) {
		crayon::yellow(
			"Reassembling large model objects! They will be saved at:\n\t",
			mod_base
		) |> rlang::inform()
		model_cat()
	} else {
		crayon::green(
			"Large model objects found at:\n\t",
			mod_base
		) |> rlang::inform()
	}
	
	# Load the models
	save_pth <- here::here(res_path, "files", "all-cates-combined.Rds")
	mod_paths <- list.files(mod_base, full.names = TRUE)
	mod_names <- stringr::str_extract(mod_paths, "^.*mod-(.{2})\\.Rds", group = 1)
	model_list <- purrr::map(
		mod_paths,
		readr::read_rds
	) |>
		rlang::set_names(mod_names)
	
	# Load the data set
	dat <- readr::read_rds(here::here("data", "processed", "model-data.Rds"))

	# Get the cate estimates
	rlang::inform(crayon::green("Starting CATE calculations."))
	all_cates_combined <- get_cates_for_all_models(dat, model_list, mod_names)
	
	# Save the cate estimates to file
	rlang::inform(c(
		crayon::green("Saving to file."),
		"i" = crayon::blue(save_pth) |> crayon::underline()
	))
	readr::write_rds(
		all_cates_combined, save_pth
	)
	
	rlang::inform(crayon::bold("All done!") |> crayon::green())
}

main()

# END OF FILE ==================================================================