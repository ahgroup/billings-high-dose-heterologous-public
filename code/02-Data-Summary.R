###
# Data summary tables and/or figures
# Zane Billings
# 2024-01-29
# This script makes Table 1 for the paper, in addition to any plots and tables
# that summarize the data but not the model fits.
###

# ---- Setup ----
# Declare package dependencies
library(readr, include.only = NULL)
library(here, include.only = NULL)
library(dplyr, include.only = NULL)
library(tidyr, include.only = NULL)
library(purrr, include.only = NULL)
library(forcats, include.only = NULL)
library(stringr, include.only = NULL)
library(tibble, include.only = NULL)
library(smd, include.only = NULL)
library(gtsummary, include.only = NULL)
library(flextable, include.only = NULL)

# Load the virus info data and renaming function
source(here::here("code", "common-functions", "virus-info.R"))

# Load data
dat <- readr::read_rds(here::here("data", "processed", "model-data.Rds")) |>
	# Replace the virus names with the short names
	dplyr::mutate(
		dplyr::across(
			c(vaccine_name, strain_name),
			replace_strain_names
		),
		ind_id = stringr::str_extract(id, "[0-9]{1,4}")
	)

# Helper functions for geometric mean and SD
geo_mean <- function(x) {
	return(x |> log() |> mean() |> exp())
}

geo_sd <- function(x) {
	return(x |> log() |> sd() |> exp())
}

# Set up gtsummary theme
gtsummary::theme_gtsummary_journal(
	journal = "jama",
	set_theme = TRUE
)
gtsummary::theme_gtsummary_language(language = "en", big.mark = "")

# ---- Study sample demographics ----
# Make a table that summarizes the features of the study population
tbl_one <-
	dat |> 
	dplyr::select(
		ind_id, study, season, dose
		#pretiter, posttiter, fold_change, seroprotection, seroconversion
	) |>
	dplyr::distinct() |>
	gtsummary::tbl_summary(
		by = dose,
		include = -ind_id,
		label = list(
			season ~ "Season",
			study ~ "Study"
		),
		statistic = list(
			gtsummary::all_continuous() ~ "{median} ({min} - {max})"
		),
		digits = list(
			gtsummary::all_continuous() ~ 0,
			gtsummary::all_categorical() ~ 0
		)
	) |>
	gtsummary::add_overall(last = TRUE) |>
	gtsummary::as_flex_table()

table_studies_seasons <-
	dat |> 
	dplyr::select(
		id, study, season, dose
		#pretiter, posttiter, fold_change, seroprotection, seroconversion
	) |>
	dplyr::distinct() |>
	gtsummary::tbl_summary(
		by = study,
		include = -id,
		label = list(
			season ~ "Season",
			dose ~ "Dose"
		),
		statistic = list(
			gtsummary::all_continuous() ~ "{median} ({min} - {max})"
		),
		digits = list(
			gtsummary::all_continuous() ~ 0,
			gtsummary::all_categorical() ~ 0
		)
	) |>
	gtsummary::add_overall(last = TRUE) |>
	gtsummary::as_flex_table()

# Make the same thing but for individuals instead of measurements
n_indiv <- dplyr::n_distinct(dat$id)
tab_indiv <-
	dat |>
	dplyr::select(id, study, age, birth_year) |>
	dplyr::summarise(
		age = min(age),
		birth_year = max(birth_year),
		.by = c(id, study)
	) |>
	gtsummary::tbl_summary(
		by = study,
		include = -id,
		label = list(
			age ~ "Age at first enrollment",
			birth_year ~ "Birth year"
		),
		statistic = list(
			gtsummary::all_continuous() ~ "{median} ({min} - {max})"
		),
		digits = list(
			gtsummary::all_continuous() ~ 0,
			gtsummary::all_categorical() ~ 0
		)
	) |>
	gtsummary::add_overall(last = TRUE) |>
	gtsummary::as_flex_table()

# Save everything as a file
readr::write_rds(
	n_indiv,
	here::here("results", "files", "n-indiv.Rds"),
	compress = "gz"
)

readr::write_rds(
	tab_indiv,
	here::here("results", "tables", "tab-indiv.Rds"),
	compress = "gz"
)

readr::write_rds(
	tbl_one,
	here::here("results", "tables", "demographics.Rds"),
	compress = "gz"
)

readr::write_rds(
	table_studies_seasons,
	here::here("results", "tables", "demographics-studies.Rds"),
	compress = "gz"
)

# ---- Strain count tables ----
# Make a table showing the number of total HAI assays performed for each
# vaccine strain/assay strain combo.
unique_counts <-
	dat |>
	dplyr::group_by(
		vaccine_name, strain_name, strain_type
	) |>
	dplyr::summarize(
		n_assays = dplyr::n(),
		n_people = dplyr::n_distinct(id),
		n_season = dplyr::n_distinct(season),
		.groups = "drop"
	) |>
	dplyr::arrange(vaccine_name, strain_name)

# Do it again, but this time ignore the assay strains to get an "overall" for
# each vaccine
count_vaccines <-
	dat |>
	dplyr::group_by(
		vaccine_name, strain_type
	) |>
	dplyr::summarize(
		strain_name = "Overall",
		n_assays = dplyr::n(),
		n_people = dplyr::n_distinct(id),
		n_season = dplyr::n_distinct(season),
		.groups = "drop"
	) |>
	dplyr::relocate(
		strain_name,
		.after = vaccine_name
	) |>
	dplyr::arrange(vaccine_name)

# Now just make an overall one
count_overall <-
	dat |>
	dplyr::summarize(
		vaccine_name = "Overall",
		strain_name = "Overall",
		n_assays = dplyr::n(),
		n_people = dplyr::n_distinct(id),
		n_season = dplyr::n_distinct(season),
		.groups = "drop"
	) |>
	dplyr::relocate(
		vaccine_name,
		strain_name,
		.before = dplyr::everything()
	)

# Take a look here to make sure that all of the strain names match correctly.
# Then we make a nice table by dropping some columns.
count_vaccines_strains <-
	# Put the by strain and overall strain counts together
	dplyr::bind_rows(
		unique_counts,
		count_vaccines,
		count_overall
	) |>
	# Mutate strain name so that overall always comes first
	dplyr::mutate(
		strain_name = relevel(
			forcats::fct_inorder(strain_name), ref = "Overall"
		),
		vaccine_name = relevel(
			forcats::fct_inorder(vaccine_name), ref = "Overall"
		)
	) |>
	dplyr::arrange(vaccine_name, strain_name) |>
	# Table formatting stuff
	dplyr::select(
		`Vaccine` = vaccine_name,
		`Assay strain` = strain_name,
		`n (assays)` = n_assays,
		`n (individuals)` = n_people
	) |>
	flextable::flextable() |>
	flextable::merge_v(j = 1) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues() |>
	flextable::autofit()

readr::write_rds(
	count_vaccines_strains,
	here::here("results", "tables", "unique-counts-vs.Rds"),
	compress = "gz"
)

# Make the same table but just for vaccines
count_vaccines_only <-
	dat |>
	dplyr::group_by(
		vaccine_name
	) |>
	dplyr::summarize(
		n_assays = dplyr::n(),
		n_people = dplyr::n_distinct(id),
		n_season = dplyr::n_distinct(season),
		.groups = "drop"
	) |>
	dplyr::arrange(vaccine_name) |>
	dplyr::select(
		`Vaccine` = vaccine_name,
		`n (assays)` = n_assays,
		`n (individuals)` = n_people,
		`n (seasons used)` = n_season
	) |>
	flextable::flextable() |>
	flextable::merge_v(j = 1) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues() |>
	flextable::autofit()

readr::write_rds(
	count_vaccines_only,
	here::here("results", "tables", "unique-counts-v.Rds"),
	compress = "gz"
)

# Make the same table but just for strains
count_strains_only <-
	dat |>
	dplyr::group_by(
		strain_name, strain_type
	) |>
	dplyr::summarize(
		n_assays = dplyr::n(),
		n_people = dplyr::n_distinct(id),
		n_season = dplyr::n_distinct(season),
		.groups = "drop"
	) |>
	dplyr::arrange(strain_name) |>
	dplyr::select(
		`Assay strain` = strain_name,
		`n (assays)` = n_assays,
		`n (individuals)` = n_people,
		`n (seasons used)` = n_season
	) |>
	flextable::flextable() |>
	flextable::fix_border_issues() |>
	flextable::autofit()

readr::write_rds(
	count_strains_only,
	here::here("results", "tables", "unique-counts-s.Rds"),
	compress = "gz"
)

# ---- Strain names and short names ----
strain_names <-
	dat |>
	dplyr::mutate(
		`Strain name` = replace_strain_names(
			strain_name,
			from = "short",
			to = "full"
		),
		`Short name` = strain_name,
		`Subtype` = strain_type,
		.keep = "none"
	) |>
	dplyr::arrange(`Strain name`) |>
	dplyr::distinct()

strain_names_h1n1 <-
	strain_names |>
	dplyr::filter(Subtype == "H1N1") |>
	dplyr::select(-Subtype) |>
	flextable::flextable()

strain_names_h3n2 <-
	strain_names |>
	dplyr::filter(Subtype == "H3N2") |>
	dplyr::select(-Subtype) |>
	flextable::flextable()

readr::write_rds(
	strain_names_h1n1,
	here::here("results", "tables", "h1n1-strain-names-table.Rds"),
	compress = "gz"
)

readr::write_rds(
	strain_names_h3n2,
	here::here("results", "tables", "h3n2-strain-names-table.Rds"),
	compress = "gz"
)

#### Vaccine Strains for Each Year ####
season_vaccines <-
	dat |>
	dplyr::select(season, vaccine_name, strain_type) |>
	dplyr::distinct() |>
	tidyr::pivot_wider(
		names_from = strain_type,
		values_from = vaccine_name
	) |>
	flextable::flextable()

readr::write_rds(
	season_vaccines,
	here::here("results", "tables", "vaccines-by-season.Rds"),
	compress = "gz"
)

# ---- Outcome summary stats ----
outcomes_continuous <-
	dat |>
	# Select only vaccine name and three outcome variables
	dplyr::select(strain_type, vaccine_name, strain_name,
								pretiter, posttiter, fold_change) |>
	# Pivot the three outcomes to long form so we can summarize all of them at
	# one time
	tidyr::pivot_longer(
		-c(strain_type, vaccine_name, strain_name),
		names_to = "outcome",
		values_to = "value"
	)

# Use this trick to get an "overall" column in the summary -- create an entire
# replicate data frame where all the vaccine names are set to "overall"
outcomes_continuous_combined <-
	dplyr::bind_rows(
		# Stratified by vaccines AND assay strains
		outcomes_continuous,
		# Stratified by just vaccine strains
		outcomes_continuous |> dplyr::mutate(strain_name = "Overall"),
		# Stratified by subtype only
		outcomes_continuous |>
			dplyr::mutate(vaccine_name = "Overall", strain_name = "Overall")
	)

outcomes_continuous_summary <-
	outcomes_continuous_combined |>
	# Compute geometric mean, geometric SD, and number of observations for each
	# outcome by vaccine strain
	dplyr::summarize(
		gm = geo_mean(value),
		gs = geo_sd(value),
		.by = c(strain_type, strain_name, vaccine_name, outcome)
	) |>
	# Paste together the GM and GS
	dplyr::mutate(
		gm = sprintf("%.2f", gm) |>
			stringr::str_pad(6, "left"),
		gs = sprintf("%.2f", gs),
		res = paste0(gm, " (±", gs, ")")
	) |>
	dplyr::select(-c(gm, gs)) |>
	# Pivot outcomes wider to have 1 col per outcome
	tidyr::pivot_wider(names_from = outcome, values_from = res)

# Do the same thing for binary outcomes
outcomes_binary <-
	dat |>
	# Select only vaccine name and three outcome variables
	dplyr::select(strain_type, vaccine_name, strain_name,
								seroprotection, seroconversion) |>
	# Pivot the three outcomes to long form so we can summarize all of them at
	# one time
	tidyr::pivot_longer(
		-c(strain_type, vaccine_name, strain_name),
		names_to = "outcome",
		values_to = "value"
	)

# Use this trick to get an "overall" column in the summary -- create an entire
# replicate data frame where all the vaccine names are set to "overall"
outcomes_binary_combined <-
	dplyr::bind_rows(
		# Stratified by vaccines AND assay strains
		outcomes_binary,
		# Stratified by just vaccine strains
		outcomes_binary |> dplyr::mutate(strain_name = "Overall"),
		# Stratified by subtype
		outcomes_binary |>
			dplyr::mutate(vaccine_name = "Overall", strain_name = "Overall")
	)

outcomes_binary_summary <-
	outcomes_binary_combined |>
	# Compute geometric mean, geometric SD, and number of observations for each
	# outcome by vaccine strain
	dplyr::summarize(
		gm = sum(value),
		gs = sum(value) / dplyr::n(),
		.by = c(strain_type, strain_name, vaccine_name, outcome)
	) |>
	# Paste together the GM and GS
	dplyr::mutate(
		gm = gm |>
			stringr::str_pad(4, "left"),
		gs = sprintf("%.1f", gs * 100) |>
			stringr::str_pad(4, "left", pad = "0"),
		res = paste0(gm, " (", gs, "%)")
	) |>
	dplyr::select(-c(gm, gs)) |>
	# Pivot outcomes wider to have 1 col per outcome
	tidyr::pivot_wider(names_from = outcome, values_from = res)

# Put the continuous and binary outcomes together
outcome_summary_combined <-
	dplyr::bind_cols(
		outcomes_continuous_summary,
		outcomes_binary_summary[, c("seroprotection", "seroconversion")]
	)

outcome_summary_clean <-
	outcome_summary_combined |>
	dplyr::mutate(
		# Have to replace the strain names to get them to sort right
		vaccine_name = replace_strain_names(
			vaccine_name,
			from = "short",
			to = "short"
		) |>
			relevel(ref = "Overall"),
		strain_name = replace_strain_names(
			strain_name,
			from = "short",
			to = "short"
		) |>
			relevel(ref = "Overall")
	) |>
	dplyr::arrange(strain_type, vaccine_name, strain_name) |>
	dplyr::select(
		`Subtype` = strain_type,
		`Vaccine strain` = vaccine_name,
		`Assay strain` = strain_name,
		`Pre-vaccination titer` = pretiter,
		`Post-vaccination titer` = posttiter,
		`Fold change` = fold_change,
		`Seroprotection` = seroprotection,
		`Seroconversion` = seroconversion
	)

outcomes_summary_table <-
	outcome_summary_clean |>
	flextable::flextable() |>
	flextable::align(j = 4:8, align = "right", part = "all") |>
	flextable::merge_v(j = 1:2) |>
	flextable::valign(j = 1:2, valign = "top") |>
	flextable::fix_border_issues() |>
	#flextable::hline(i = c(5, 12)) |>
	flextable::footnote(
		i = 1,
		j = 3,
		part = "header",
		ref_symbol = "1",
		value = flextable::as_paragraph(c("Number of total HAI assays"))
	) |>
	flextable::footnote(
		i = 1,
		j = 4:6,
		part = "header",
		ref_symbol = "2",
		value = flextable::as_paragraph(c("Geometric mean (±geometric SD)"))
	) |>
	flextable::autofit()

readr::write_rds(
	outcomes_summary_table,
	here::here("results", "tables", "vaccine-outcomes.Rds"),
	compress = "gz"
)

# ---- Crude analyses ----
# For each of the four outcomes, and also pre-vaccination titer, we examined
# crude summary statistics and SMDs.
# Summaries for pre-titer, post-titer, and fold change are based on geometric
# mean and SD, while summaries for seroprotection and seroconversion are based
# on counts and proportions.

# This function automates the crude analysis.
# Stat arg: "am" = arithmetic mean, "gm" = geometric mean, "n" = n/% for binary
# Format arg: "df" = raw calculations, "d2" = cleaned up / formatted calcs,
#    "ft" = flextable object
# Outcome should be passed as a bare name since tidyselection is used
# "Label" will be a header above the statistics in the output table and
# "cpt" will become a footnote for those statistics.
# "alpha" is the alpha level for the two-sided CI for the SMD.
crude_analysis <- function(.data, outcome, label, cpt, stat = c("am", "gm", "n"),
													 format = c("df", "d2", "ft"), alpha = 0.05) {
	subset_data <- .data |>
		dplyr::select(strain_type, vaccine_name, dose, {{outcome}})
	
	# Select the correct summary functions
	if (stat == "am") {
		mean_fun = mean
		sd_fun = sd
	} else if (stat == "gm") {
		mean_fun = geo_mean
		sd_fun = geo_sd
	} else if (stat == "n") {
		mean_fun = sum
		sd_fun = mean
	} else {
		stop("invalid 'stat' argument")
	}
	
	# Sometimes we need to take the log of the outcome before it gets passed to
	# the SMD function. The maybe_log() function will evaluate as log() if that
	# is specified in the arguments, but will evaluated as identity() if the
	# log is not specified in the arguments.
	if (stat == "gm") {
		maybe_log <- log
	} else {
		maybe_log <- identity
	}
	
	# Annoying trick to add an "overall" dose group to the summary
	dose_data <-
		dplyr::bind_rows(
			subset_data |> dplyr::mutate(
				vaccine_name = "Overall",
				strain_type = "Overall"
			),
			subset_data |> dplyr::mutate(
				vaccine_name = "Overall"
			),
			subset_data
		)

	# Calculate summary stats 	
	summary_data <-
		dose_data |>
		dplyr::group_by(strain_type, vaccine_name, dose) |>
		dplyr::summarize(
			n = dplyr::n(),
			m = mean_fun({{outcome}}),
			s = sd_fun({{outcome}}),
			.groups = "drop"
		) |>
		tidyr::pivot_wider(
			names_from = dose,
			values_from = c(n, m, s),
			names_glue = "{dose}_{.value}"
		)
	
	# Calculate the SMD -- we have to use nesting
	dose_smd <-
		dose_data |>
		tidyr::nest(smd_data = -c(strain_type, vaccine_name)) |>
		dplyr::mutate(
			smd = purrr::map(
				smd_data,
				\(d) smd::smd(
					x = maybe_log(dplyr::pull(d, {{outcome}})),
					g = dplyr::pull(d, dose),
					std.error = TRUE,
					gref = 2L
				)
			)
		)
	
	# Calculate critical value for SMD CI
	cz <- qnorm(1 - alpha / 2)
	
	# Add the SMDs to the summary data
	comb_data <- dplyr::left_join(
		dose_smd,
		summary_data,
		by = c("strain_type", "vaccine_name")
	) |>
		tidyr::unnest(smd) |>
		dplyr::select(-term, -smd_data)
	
	if (format == "df") {return(comb_data)}
	
	# If stat == "n" then these cols are integers and don't need to be rounded
	# And we need to format the _s differnetly.
	if (stat != "n") {
		comb_data2 <-
			comb_data |>
			dplyr::mutate(
				dplyr::across(
					c(SD_m, HD_m),
					\(x) sprintf("%.2f", x) |> stringr::str_pad(width = 4)
				),
				dplyr::across(
					c(SD_s, HD_s),
					\(x) sprintf("%.2f", x) |> stringr::str_pad(width = 4)
				),
				SD_s = paste0("±", SD_s),
				HD_s = paste0("±", HD_s)
			)
	} else {
		comb_data2 <-
			comb_data |>
			dplyr::mutate(
				dplyr::across(
					c(SD_s, HD_s),
					\(x) sprintf("%.0f", x * 100) |>
						stringr::str_pad(width = 2, pad = "0")
				),
				SD_s = paste0(SD_s, "%"),
				HD_s = paste0(HD_s, "%")
			)
	}
	
	# Clean up the data
	tbl_data <-
		comb_data2 |>
		dplyr::mutate(
			lwr = sprintf("%.2f", estimate - cz * std.error),
			upr = sprintf("%.2f", estimate + cz * std.error),
			estimate = sprintf("%.2f", estimate),
			dplyr::across(
				c(SD_m, HD_m, estimate, lwr, upr),
				\(x) ifelse(x == "-0.00", "0.00", x)
			)
		) |>
		dplyr::mutate(
			SD_v = paste0(SD_m, " (", SD_s, ")"),
			HD_v = paste0(HD_m, " (", HD_s, ")"),
			smd = paste0(estimate, " (", lwr, ", ", upr, ")"),
			.keep = "unused"
		) |>
		dplyr::select(-std.error) |>
		dplyr::arrange(
			forcats::fct_inorder(strain_type),
			forcats::fct_inorder(vaccine_name)
		)
	
	if (format == "d2") {return(tbl_data)}	
	
	tbl <-
		tbl_data |>
		dplyr::rename(
			"Subtype" = strain_type,
			"Vaccine strain" = vaccine_name,
			"SMD" = smd,
			"n_SD" = SD_n,
			"n_HD" = HD_n,
			"v_SD" = SD_v,
			"v_HD" = HD_v
		) |>
		dplyr::rename_with(
			.fn = \(x) stringr::str_replace(x, "^v_", paste0(label, "_")),
			.cols = dplyr::starts_with("v_")
		) |>
		flextable::flextable() |>
		flextable::separate_header(
			split = "_", fixed = TRUE,
			opts = c("span-top", "center-hspan", "bottom-vspan")
		) |>
		flextable::align(i = 2, j = 3:7, align = "center", part = "header") |>
		flextable::align(j = 3:6, align = "right", part = "body") |>
		flextable::merge_v(j = 1) |>
		flextable::valign(j = 1:2, valign = "top") |>
		flextable::fix_border_issues() |>
		flextable::hline(i = c(1, 7)) |>
		flextable::footnote(
			i = 1,
			j = 3,
			part = "header",
			ref_symbol = "1",
			value = flextable::as_paragraph(paste0(
				"Total number of HAI assays across all assay strains and seasons."
			))
		) |>
		flextable::footnote(
			i = 1,
			j = 5,
			part = "header",
			ref_symbol = "2",
			value = flextable::as_paragraph(cpt)
		) |>
		flextable::footnote(
			i = 1,
			j = 7,
			part = "header",
			ref_symbol = "3",
			value = flextable::as_paragraph(paste0(
				"Standardized mean difference (HD - SD); SMD (95% CI)."
			))
		) |>
		flextable::autofit()

	if (format == "ft") {return(tbl)} else (stop("invalid 'format' argument"))	
}

## ---- Prevaccination titer ----
pret_dose_tbl <- crude_analysis(
	dat, outcome = pretiter,
	label = "Pre-vaccination titer",
	stat = "gm", format = "ft",
	cpt = paste0(
		"Pre-vaccination HAI titer. Geometric mean ",
		"(± geometric standard deviation)."
	)
)

readr::write_rds(
	pret_dose_tbl,
	here::here("results", "tables", "dose-pret-comparison.Rds"),
	compress = "gz"
)

## ---- Postvaccination titer ----
postt_dose_tbl <- crude_analysis(
	dat, outcome = posttiter,
	label = "Post-vaccination titer",
	stat = "gm", format = "ft",
	cpt = paste0(
		"Post-vaccination HAI titer. Geometric mean ",
		"(± geometric standard deviation)."
	)
)

readr::write_rds(
	postt_dose_tbl,
	here::here("results", "tables", "dose-postt-comparison.Rds"),
	compress = "gz"
)

# ---- Fold change by dose table ----
fc_dose_tbl <- crude_analysis(
	dat, outcome = fold_change,
	label = "Fold change",
	stat = "gm", format = "ft",
	cpt = paste0(
		"Fold change (post-vaccination titer divided by pre-vaccination titer. ",
		"Geometric mean (± geometric standard deviation)."
	)
)

readr::write_rds(
	fc_dose_tbl,
	here::here("results", "tables", "dose-fc-comparison.Rds"),
	compress = "gz"
)

# ---- Seroprotection by dose table ----
sp_dose_tbl <- crude_analysis(
	dat, outcome = seroprotection,
	label = "Seroprotection",
	stat = "n", format = "ft",
	cpt = paste0(
		"Seroprotection (indicator for post-titer ≥1:40); ",
		"n (%)."
	)
)

readr::write_rds(
	sp_dose_tbl,
	here::here("results", "tables", "dose-sp-comparison.Rds"),
	compress = "gz"
)

# ---- Seroconversion by dose table ----
sc_dose_tbl <- crude_analysis(
	dat, outcome = seroconversion,
	label = "Seroconversion",
	stat = "n", format = "ft",
	cpt = paste0(
		"Seroconversion (indicator for post-titer ≥1:40 only after vaccination); ",
		"n (%)."
	)
)

readr::write_rds(
	sc_dose_tbl,
	here::here("results", "tables", "dose-sc-comparison.Rds"),
	compress = "gz"
)

# ---- Raw data figure -----
library(ggplot2)
# A little bit of preprocessing before plotting
set.seed(101)
plot_dat <-
	dat |>
	dplyr::select(
		id, season, vaccine_name, strain_name, strain_type,
		pretiter, posttiter, dose
	) |>
	# Get all HAI values as one column so pre/post can go on x axis
	tidyr::pivot_longer(
		c(pretiter, posttiter),
		names_to = "timing",
		values_to = "HAI"
	) |>
	dplyr::mutate(
		# Clean up x axis labels to save space
		timing = factor(
			timing,
			levels = c("pretiter", "posttiter"),
			labels = c("Pre", "Post")
		),
		dose_timing = factor(
			paste(dose, timing),
			levels = c(
				"SD Pre",
				"SD Post",
				"HD Pre",
				"HD Post"
			),
			labels = c(
				"Pre (SD)",
				"Post (SD)",
				"Pre (HD)",
				"Post (HD)"
			)
		),
		# Add manual jitter -- ggplot converts factor x axis to integer anyways
		# so this is how it would do it. We do it manually to ensure lines and
		# points get the same amount of jitter.
		# For this x axis points around 0 will be Pre and around 1 will be Post,
		# we can manually change the plot labels to words.
		x_j = as.integer(dose_timing) + runif(dplyr::n(), -0.25, +0.25),
		# Also jitter the HAIs, but we have to do it on the log scale since that's
		# what we're plotting, and then put it back on natural scale to pass to
		# ggplot or else it can't create the log scale correctly.
		HAI_j = 5 * 2 ^ (log2(HAI / 5) + runif(dplyr::n(), -0.15, +0.15)),
		# Create a grouping variable to make sure person-years all get their own
		# point and line
		unique_id = paste(id, season, sep = "_")
	)

# Function to create one raw data plot
plot_fcn <- function(data, title) {
	plt <-
		data |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = x_j, y = HAI_j, group = unique_id,
			shape = dose, color = dose
		) +
		ggplot2::geom_line(
			color = "black",
			alpha = 0.1
		) +
		ggplot2::geom_point(
			size = 3,
			alpha = 0.25
		) +
		ggplot2::facet_wrap(~strain_name, ncol = 4, scales = "free_x") +
		ggplot2::scale_x_continuous(
			breaks = c(1, 2, 3, 4),
			minor_breaks = NULL,
			expand = c(0, 0),
			limits = c(0.5, 4.5),
			labels = c("Pre\nSD", "Post\nSD", "Pre\nHD", "Post\nHD")
		) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = 5 * 2 ^ seq(0, 10, 2)
		) +
		ggplot2::scale_color_brewer(palette = "Dark2") +
		hgp::theme_ms() +
		ggplot2::labs(
			x = NULL,
			y = "HAI titer",
			title = paste0("Vaccine: ", title)
		) +
		ggplot2::guides(
			color = ggplot2::guide_legend(override.aes = list(alpha = 1))
		)
	
	return(plt)
}

# Create a separate raw data plot for each vaccine strain
raw_data_plot_df <-
	plot_dat |>
	tidyr::nest(plt_data = -vaccine_name) |>
	dplyr::mutate(plt = purrr::map2(plt_data, vaccine_name, plot_fcn))

# Save raw data plots as png
purrr::walk2(
	raw_data_plot_df$plt, raw_data_plot_df$vaccine_name,
	\(plt, y) ggplot2::ggsave(
		here::here(
			"results", "figures", "raw-data",
			paste0(gsub("/", "-", y), ".png")
		),
		plot = plt,
		width = 13,
		height = 14,
		units = "in"
	)
)

# # Save raw data plots as tiff
# purrr::walk2(
# 	raw_data_plot_df$plt, raw_data_plot_df$vaccine_name,
# 	\(plt, y) ggplot2::ggsave(
# 		here::here(
# 			"results", "figures", "raw-data",
# 			paste0(gsub("/", "-", y), ".tiff")
# 		),
# 		plot = plt,
# 		width = 12,
# 		height = 8,
# 		units = "in"
# 	)
# )

# Raw fold change plots ----
# A little bit of preprocessing before plotting
fc_plot_dat <-
	dat |>
	dplyr::select(
		id, season, vaccine_name, strain_name, strain_type,
		fold_change, dose
	)

# Function to create one raw data plot
make_fc_plot <- function(data, title) {
	plt <-
		data |>
		ggplot2::ggplot() +
		ggplot2::aes(x = dose, y = fold_change, color = dose, shape = dose) +
		ggplot2::geom_point(
			size = 3,
			alpha = 0.25,
			show.legend = FALSE,
			position = ggplot2::position_jitter(0.125, 0.125, 123123)
		) +
		ggplot2::scale_y_continuous(
			trans = "log2",
			breaks = 2 ^ seq(-10, 10, 2),
			limits = 2 ^ c(-7, 11),
			minor_breaks = 2 ^ seq(-10, 10, 1),
			labels = MASS::fractions
		) +
		ggplot2::facet_wrap(~strain_name, ncol = 7, scales = "free_x") +
		ggplot2::scale_color_brewer(palette = "Dark2") +
		hgp::theme_ms(text_size_axis_text = 12) +
		ggplot2::labs(
			x = NULL,
			y = "Fold change (post titer / pre titer)",
			title = paste0("Vaccine: ", title)
		) +
		ggplot2::guides(
			color = ggplot2::guide_legend(override.aes = list(alpha = 1))
		)
	
	return(plt)
}

# Create a separate raw data plot for each vaccine strain
fc_data_plot_df <-
	fc_plot_dat |>
	tidyr::nest(plt_data = -vaccine_name) |>
	dplyr::mutate(plt = purrr::map2(plt_data, vaccine_name, make_fc_plot))

# Save fc data plots as png
purrr::walk2(
	fc_data_plot_df$plt, fc_data_plot_df$vaccine_name,
	\(plt, y) ggplot2::ggsave(
		here::here(
			"results", "figures", "raw-fc",
			paste0(gsub("/", "-", y), ".png")
		),
		plot = plt,
		width = 13,
		height = 14,
		units = "in"
	)
)

# # Save fc data plots as tiff
# purrr::walk2(
# 	fc_data_plot_df$plt, fc_data_plot_df$vaccine_name,
# 	\(plt, y) ggplot2::ggsave(
# 		here::here(
# 			"results", "figures", "raw-fc",
# 			paste0(gsub("/", "-", y), ".tiff")
# 		),
# 		plot = plt,
# 		width = 12,
# 		height = 8,
# 		units = "in"
# 	)
# )



# END OF FILE ==================================================================
