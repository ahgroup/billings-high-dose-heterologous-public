###
# Supplementary Analyses
# Zane
# 2024-06-06
# Produces figures and tables for the supplement
###

# This script sets up convenient project paths
source(here::here("code", "common-functions", "path-setting.R"))

# Load the virus info data and renaming function
source(here::here("code", "common-functions", "virus-info.R"))

# Load the data
dat <- readr::read_rds(here::here("data", "processed", "model-data.Rds")) |>
	# Replace the virus names with the short names
	dplyr::mutate(
		dplyr::across(
			c(vaccine_name, strain_name),
			replace_strain_names
		)
	)
 
# DAG ####
library(dagitty)
library(ggdag)
library(ggplot2)

dag <- dagitty::dagitty(
	"dag {
	d <- {a t U i b} -> y
	{d p sv sa} -> y
	d [exposure]
	y [outcome]
	}")
dagitty::coordinates(dag) <- list(
	x = c(y = 5, d = 0, a = 1, t = 2, i = 3, U = 1, b = 2, p = 2.5,
				sv = 3, sa = 3.5),
	y = c(y = 0, d = 0, a = 1, t = 1, i = 1, U = -1, b = -1, p = -1,
				sv = -1, sa = -1)
)

tidy_dag <- ggdag::tidy_dagitty(dag)
tidy_dag$data <-
	tidy_dag$data |>
	dplyr::mutate(lab = dplyr::case_match(
		name,
		"U" ~ "U",
		"a" ~ "Age",
		"b" ~ "Birth year",
		"d" ~ "Dose",
		"i" ~ "Individual effect",
		"p" ~ "Pre-vaccination titer",
		"sv" ~ "Vaccine strain",
		"sa" ~ "Assay strain",
		"t" ~ "Calendar time",
		"y" ~ "Outcome"
	))

my_ggdag <-
	tidy_dag |>
	ggdag::node_status() |>
	dplyr::mutate(
		observed = factor(
			ifelse(name %in% c("U", "i"), 0, 1),
			levels = c(1, 0),
			labels = c("Observed", "Unobserved")
		),
		status = factor(
			ifelse(is.na(status), "bleh", status),
			levels = c("1", "2", "bleh"),
			labels = c("Exposure", "Outcome", "Covariate")
		)
	) |>
	ggdag::ggdag() +
	ggdag::geom_dag_point(aes(color = status, shape = observed)) +
	ggdag::geom_dag_text() +
	ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00", "black")) +
	ggplot2::scale_shape_manual(values = c(16, 15)) +
	ggdag::theme_dag()
ggsave(
	here::here("results", "figures", "dag.tiff"),
	plot = my_ggdag,
	width = 10,
	height = 5
)

ggsave(
	here::here("results", "figures", "dag.png"),
	plot = my_ggdag,
	width = 10,
	height = 5
)

# Table of assay strains by seasons ####
season_counts_dat <-
	dat |>
	dplyr::group_by(season, strain_name, strain_type) |>
	dplyr::summarize(
		n_assays = dplyr::n(),
		.groups = "drop"
	) |>
	# Format season labels to short names
	dplyr::mutate(
		season = gsub("20([0-9][0-9])", "\\1", gsub(" - ", "/", season))
	) |>
	tidyr::pivot_wider(
		names_from = season,
		values_from = n_assays
	) |>
	# Replace NAs with empty strings
	dplyr::mutate(
		dplyr::across(-c(strain_name, strain_type), \(x) ifelse(is.na(x), 0, x))
	) |>
	dplyr::relocate(strain_type, .before = strain_name) |>
	dplyr::arrange(strain_type, forcats::fct_inorder(strain_name)) |>
	# Add the total column at the end
	dplyr::rowwise() |>
	dplyr::mutate(Total = sum(dplyr::c_across(`13/14`:`21/22`))) |>
	dplyr::ungroup()

season_counts_tab <-
	season_counts_dat |>
	flextable::flextable() |>
	flextable::set_header_labels(
		values = list(
			strain_type = "Subtype",
			strain_name = "Strain"
		)
	) |>
	flextable::merge_v(j = 1) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues()

readr::write_rds(
	season_counts_tab,
	here::here("results", "tables", "assay-seasons-table.Rds")
)

# Printing model diagnostics ####
mod_paths <- list.files(here::here(largefile_path, "mods"), full.names = TRUE)
model_diagnostics <-
	purrr::map(
		mod_paths,
		\(p) {
			cur_mod <- readr::read_rds(p)
			cur_summary <- posterior::summarise_draws(cur_mod$fit)
			out <- tibble::tibble(
				Model = p,
				"Num. Divergences" = paste0(
					rstan:::get_num_divergent(cur_mod$fit), " / 20000"
					),
				"min E-BFMI" = rstan::get_bfmi(cur_mod$fit) |> min(),
				"min ESS (tail)" = min(cur_summary$ess_tail),
				"min ESS (bulk)" = min(cur_summary$ess_bulk),
				"max R_hat" = max(cur_summary$rhat)
			)
		}
	) |>
	purrr::list_rbind() |>
	dplyr::mutate(
		dplyr::across(tidyselect::contains("ESS"), \(x) sprintf("%.0f", x)),
		dplyr::across(c(`min E-BFMI`, `max R_hat`), \(x) sprintf("%.3f", x))
	)

model_diagnostics <- model_diagnostics[c(4, 1, 3, 2), ]

model_diagnostics$Model <- c(
	"Titer increase", "Post-vaccination titer", "Seroprotection", "Seroconversion"
)

model_diagnostics <- flextable::flextable(model_diagnostics)

readr::write_rds(model_diagnostics, here::here("results", "tables", "mod-diag.Rds"))

# END OF FILE ####
