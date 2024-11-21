####
# Diagnostic Plots for Model Parameters
# Zane Billings
# 2024-11-02
# Creates diagnostic plots for some of the model parameters to include in the
# supplement. The code is complicated and memory intensive so gets its own
# script instead of going in supplementary analyses.
####

# Setup ####
## Packages and script dependencies ####
library(ggplot2)
library(brms)
library(cmdstanr)
library(bayesplot)

source(here::here("code", "common-functions", "path-setting.R"))
source(here::here("code", "common-functions", "model-manip.R"))

## Set up file paths ####
# Run the model de-chunker if necessary
# First do it for the posterior models
post_base_pth <- here::here(largefile_path, "mods")
if (length(list.files(post_base_pth)) != 4) {
	crayon::yellow(
		"Reassembling large posterior model objects! They will be saved at:\n\t",
		post_base_pth
	) |> rlang::inform()
	model_cat()
} else {
	crayon::green(
		"Large model objects found at:\n\t",
		post_base_pth
	) |> rlang::inform()
}

# Now do it for the prior models
pre_base_pth <- here::here(largefile_path, "priors")
if (length(list.files(pre_base_pth)) != 4) {
	crayon::yellow(
		"Reassembling large prior model objects! They will be saved at:\n\t",
		pre_base_pth
	) |> rlang::inform()
	model_cat(priors = TRUE)
} else {
	crayon::green(
		"Large model objects found at:\n\t",
		pre_base_pth
	) |> rlang::inform()
}

## Load the models ####
post_paths <- list.files(here::here(largefile_path, "mods"), full.names = TRUE)
post_models <- purrr::map(post_paths, \(f) readr::read_rds(f))

pre_paths <- list.files(here::here(largefile_path, "priors"), full.names = TRUE)
pre_models <- purrr::map(pre_paths, \(f) readr::read_rds(f))

# Printing model diagnostics ####
# Get the total number of divergent transitions (formatted as fraction of
# total samples), the minimum E-BFMI, tail ESS, and bulk ESS across all model
# parameters, and the maximum R_hat4 across all model parameters.
# Gets those diagnostics for each of the four models.
model_diagnostics_raw <-
	purrr::map(
		post_models,
		\(p) {
			cur_summary <- posterior::summarise_draws(p$fit)
			out <- tibble::tibble(
				"Num. Divergences" = paste0(
					rstan:::get_num_divergent(p$fit), " / ",
					brms::ndraws(p)
				),
				"min E-BFMI" = rstan::get_bfmi(p$fit) |> min(),
				"min ESS (tail)" = min(cur_summary$ess_tail),
				"min ESS (bulk)" = min(cur_summary$ess_bulk),
				"max R_hat" = max(cur_summary$rhat)
			)
		}
	)

# Create a table by binding that list together, and then do some formatting
# to make it look presentable
model_names_df <- data.frame(
	full = c(
		"Post-vaccination titer",
		"Seroconversion",
		"Seroprotection",
		"Titer increase"
	),
	short = c("pt", "sc", "sp", "ti")
)
model_diagnostics <-
	model_diagnostics_raw |>
	purrr::list_rbind() |>
	dplyr::mutate(
		dplyr::across(tidyselect::contains("ESS"), \(x) sprintf("%.0f", x)),
		dplyr::across(c(`min E-BFMI`, `max R_hat`), \(x) sprintf("%.3f", x))
	) |>
	dplyr::mutate(
		Model = model_names_df$full,
		.before = dplyr::everything()
	)

# Change the table format and export so we can load the table in the
# manuscript document
model_diagnostics_table <- flextable::flextable(model_diagnostics)
readr::write_rds(
	model_diagnostics_table,
	here::here("results", "tables", "mod-diag.Rds")
)

# Select parameters for diagnostics ####
# These are the parameters that we want to extract from the brms model. You can
# get the full list by brms::variables(post_models[[1]]) and add or remove
# as many as you want.
parms_to_select <- c(
	"Intercept", # The population-level intercept term
	"b_doseHD", # Population-level effect of HD
	"bs_sbirth_year_c_1", # First spline coefficient term
	"sd_id__Intercept", # First random intercept variance term
	"sd_strain_type__Intercept", # First random intercept variance w/ correlation
	"sd_strain_type__doseHD", # First random slope
	"cor_strain_type__Intercept__doseHD", # First random effect correlation
	"sds_sage_c_1" # First sds-class term.
)

# Go through all of the posterior sample models (all 4 of them), and extract
# the draws in a tidy format.
post_draws <-
	purrr::map(
		post_models,
		\(m) m |>
			tidybayes::tidy_draws() |>
			dplyr::select(
				dplyr::starts_with("."),
				dplyr::all_of(parms_to_select)
			)
	) |>
	rlang::set_names(model_names_df$full) |>
	purrr::list_rbind(names_to = ".model")

# Same thing as above but for the prior sample models
pre_draws <-
	purrr::map(
		pre_models,
		\(m) m |>
			tidybayes::tidy_draws() |>
			dplyr::select(
				dplyr::starts_with("."),
				dplyr::all_of(parms_to_select)
			)
	) |>
	rlang::set_names(model_names_df$full) |>
	purrr::list_rbind(names_to = ".model")

# Put the prior and posterior sample models into the same data frame for
# easier plotting later
combined_draws_wide <-
	dplyr::bind_rows(
		"prior" = pre_draws,
		"posterior" = post_draws,
		.id = ".samples"
	) |>
	dplyr::mutate(
		.samples = factor(
			.samples,
			levels = c("prior", "posterior"),
			labels = c("Prior", "Posterior")
		)
	)

# Create a long form for ggplot so we can facet on the different parameters
combined_draws_long <-
	combined_draws_wide |>
	tidyr::pivot_longer(
		cols = dplyr::all_of(parms_to_select),
		names_to = "parameter",
		values_to = "value"
	)

# Free up some memory -- we could avoid this by doing the combined draws in a
# local scope or just not saving the intermediates, but the way the code is
# written now is easier to read.
# On a computer with enough (v)RAM we don't need to do this, but we won't use
# these intermediate objects again so might as well get rid of them. Manual
# invocation of the garbage collector is not necessary but ensures those things
# release memory now instead of in the future.
rm(post_draws, pre_draws)
invisible(gc())

# Trank plots ####
# First we need to create the trace ranks.
BIN_WIDTH <- 500
trace_rank_data <-
	combined_draws_long |>
	dplyr::filter(.samples == "Posterior") |>
	dplyr::group_by(.model, parameter) |>
	dplyr::mutate(
		value_rank = dplyr::min_rank(value),
		value_bin = cut(
			value_rank,
			breaks = seq(0, 20000 + BIN_WIDTH, BIN_WIDTH),
			labels = seq(0, 20000, BIN_WIDTH)
		),
		value_int = as.numeric(as.character(value_bin))
	) |>
	dplyr::group_by(.model, parameter, .chain) |>
	dplyr::add_count(value_bin) |>
	dplyr::ungroup() |>
	dplyr::select(-.samples) |>
	tidyr::nest(plt_dat = -.model)

# Function to make a trank plot for a single model
make_trank_plot <- function(data, model, dir = NULL) {
	out <- data |>
	ggplot2::ggplot() +
		ggplot2::aes(
			x = value_int, y = n,
			color = factor(.chain), group = factor(.chain)
		) +
		ggplot2::geom_step(alpha = 0.85) +
		ggplot2::facet_wrap(~parameter, ncol = 4, nrow = 2) +
		ggplot2::scale_color_viridis_d(option = "mako", begin = 0.1, end = 0.9) +
		ggplot2::coord_cartesian(
			xlim = c(1, 20000),
			ylim = c(0, 60)
		) +
		ggplot2::labs(
			x = "Rank (across chains)",
			y = "Count"
		) +
		ggplot2::guides(
			color = ggplot2::guide_legend(
				title = "Chain",
				nrow = 2
			)
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			strip.text = ggplot2::element_text(size = 12),
			legend.text = ggplot2::element_text(size = 18),
			legend.title = ggplot2::element_text(size = 18),
			legend.key.size = ggplot2::unit(0.25, "in")
		)
	
	# Set up saving directory
	if (!is.null(dir)) {
		file_path <- here::here(dir, paste0("mcmc-trank-", model, ".png"))
		ggplot2::ggsave(
			file_path,
			plot = out,
			width = 13,
			height = 6.5,
			dpi = 300
		)
	}
	
	invisible(out)
}

trank_plots <-
	purrr::map2(
		trace_rank_data$plt_dat,
		model_names_df$short,
		\(x, y) make_trank_plot(x, y, dir = here::here("results", "figures"))
	)

# Density plots ####
dens_plot_data <- 
	combined_draws_long |>
	tidyr::nest(plt_dat = -.model)

make_pp_density_plot <- function(data, model, dir = NULL) {
	out <- data |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = value,
			fill = .samples
		) +
		ggplot2::geom_density(
			position = ggplot2::position_identity(),
			alpha = 0.65
		) +
		ggplot2::facet_wrap(
			ggplot2::vars(parameter),
			ncol = 4, nrow = 2,
			scales = "free"
		) +
		labs(
			x = "Parameter value",
			y = "Estimated density",
			fill = NULL
		) +
		ggplot2::scale_fill_manual(
			values = c("#56B4E9", "#CC79A7")
		) +
		hgp::theme_ms() +
		ggplot2::theme(
			strip.text = ggplot2::element_text(size = 11),
			legend.text = ggplot2::element_text(size = 18),
			legend.title = ggplot2::element_text(size = 18),
			legend.key.size = ggplot2::unit(0.25, "in")
		)
	
	# Set up saving directory
	if (!is.null(dir)) {
		file_path <- here::here(dir, paste0("mcmc-dens-", model, ".png"))
		ggplot2::ggsave(
			file_path,
			plot = out,
			width = 13,
			height = 6.5,
			dpi = 300
		)
	}
}

dens_plots <-
	purrr::map2(
		dens_plot_data$plt_dat,
		model_names_df$short,
		\(x, y) make_pp_density_plot(x, y, dir = here::here("results", "figures"))
	)

# END OF SCRIPT ####
