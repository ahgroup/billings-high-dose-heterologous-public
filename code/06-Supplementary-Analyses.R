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

# Set up gtsummary theme
gtsummary::theme_gtsummary_journal(
	journal = "jama",
	set_theme = TRUE
)
gtsummary::theme_gtsummary_language(language = "en", big.mark = "")

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
requireNamespace("dagitty")
requireNamespace("ggdag")
requireNamespace("ggplot2")

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
	ggdag::geom_dag_point(ggplot2::aes(color = status, shape = observed)) +
	ggdag::geom_dag_text() +
	ggplot2::scale_color_manual(values = c("#56B4E9", "#E69F00", "black")) +
	ggplot2::scale_shape_manual(values = c(16, 15)) +
	ggdag::theme_dag() +
	ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white"))

ggplot2::ggsave(
	here::here("results", "figures", "dag.tiff"),
	plot = my_ggdag,
	width = 10,
	height = 5
)

ggplot2::ggsave(
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

season_counts_dat_h1n1 <-
	season_counts_dat |>
	dplyr::filter(strain_type == "H1N1") |>
	dplyr::select(-strain_type)

season_counts_dat_h3n2 <-
	season_counts_dat |>
	dplyr::filter(strain_type == "H3N2") |>
	dplyr::select(-strain_type)

season_counts_tab_h1n1 <-
	season_counts_dat_h1n1 |>
	flextable::flextable() |>
	flextable::set_header_labels(
		values = list(
			strain_name = "Strain"
		)
	) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues()

readr::write_rds(
	season_counts_tab_h1n1,
	here::here("results", "tables", "h1n1-assay-seasons-table.Rds")
)

season_counts_tab_h3n2 <-
	season_counts_dat_h3n2 |>
	flextable::flextable() |>
	flextable::set_header_labels(
		values = list(
			strain_name = "Strain"
		)
	) |>
	flextable::valign(j = 1, valign = "top") |>
	flextable::fix_border_issues()

readr::write_rds(
	season_counts_tab_h3n2,
	here::here("results", "tables", "h3n2-assay-seasons-table.Rds")
)

# Dose as outcome analysis ####
# Analyze whether any patient characteristics were more likely to get HD.
requireNamespace("brms")
requireNamespace("broom.mixed")
requireNamespace("cmdstanr")

source(here::here("code", "common-functions", "brms-settings.R"))

# Get the extra race/ethnicity/etc variables
dat_hsdi <- readr::read_rds(here::here("data", "processed", "reporting-data.Rds"))
dose_analysis_data <-
	dat_hsdi |>
	dplyr::transmute(
		"study" = factor(
			gsub("_.*$", "", Cohort_ID),
			levels = c("UGA", "PA", "FL")
		),
		"season" = stringr::str_extract(Cohort_ID, "[0-9]{4}"),
		season = factor(
			paste0(season, " - ", as.integer(season) + 1)
		),
		"sex" = Sex_Assigned_at_Birth,
		"age" = cut(
			Min_Age,
			breaks = seq(65, 85, 5),
			include.lowest = TRUE,
			labels = paste0(
				c(65, 71, 76, 81), " - ", c(70, 75, 80, 85)
			)
		),
		"by" = cut(
			Birth_Year,
			breaks = seq(1930, 1960, 5),
			include.lowest = TRUE,
			labels = paste0(
				c("1930", seq(1935, 1955, 5) + 1), " - ",
				seq(1935, 1960, 5)
			)
		),
		"raceeth" = factor(
			dplyr::case_when(
				Race == "White" & startsWith(as.character(Ethnicity), "Not") ~ "W",
				startsWith(as.character(Race), "Black") &
					startsWith(as.character(Ethnicity), "Not") ~ "B",
				TRUE ~ "O"
			),
			levels = c("W", "B", "O"),
			labels = c(
				"White or Caucasian",
				"Black or African American",
				"Other"
			)
		),
		"Dose" = ifelse(grepl("SD", Cohort_ID), 0, 1)
	)

if (interactive() |> isTRUE()) {
	summary(dose_analysis_data)
}

regression_tbl_dose_fast <-
	dose_analysis_data |>
	gtsummary::tbl_uvregression(
		method = glm,
		method.args = list(
			family = "binomial"
		),
		exponentiate = TRUE,
		y = Dose,
		label = list(
			study ~ "Study",
			season ~ "Season",
			sex ~ "Sex assigned at birth",
			age ~ "Age at enrollment",
			by ~ "Birth year",
			raceeth ~ "Race/Ethnicity"
		)
	)

regression_tbl_dose_bayes <-
	dose_analysis_data |>
	gtsummary::tbl_uvregression(
		method = brms::brm,
		method.args = list(
			family = "bernoulli",
			cores = 4,
			backend = "cmdstanr"
		),
		exponentiate = TRUE,
		y = Dose,
		label = list(
			study ~ "Study",
			season ~ "Season",
			sex ~ "Sex assigned at birth",
			age ~ "Age at enrollment",
			by ~ "Birth year",
			raceeth ~ "Race/Ethnicity"
		),
		tidy_fun = \(...) broom.mixed:::tidy.brmsfit(
			...,
			conf.method = "HPDinterval"
		)
	)

N_dat <- nrow(dose_analysis_data)
N_indiv <- dplyr::n_distinct(dat_hsdi$Subject_ID)

regression_tbl_dose_fmt <-
	regression_tbl_dose_bayes |>
	gtsummary::modify_header(
		label = "**Variable**",
		estimate = "**OR**"
	) |>
	gtsummary::modify_column_hide(stat_n) |>
	gtsummary::modify_footnote(
		estimate ~ "OR = Odds Ratio"
	)

regression_tbl_dose_flextable <-
	regression_tbl_dose_fmt |>
	gtsummary::as_flex_table() |>
	# Long annoying process of manually fixing cells since apparently we can't
	# format brms stuff correctly
	# First is the seasons
	flextable::compose(
		i = 7, j = 1,
		flextable::as_paragraph("2014 - 2015")
	) |>
	flextable::compose(
		i = 8, j = 1,
		flextable::as_paragraph("2015 - 2016")
	) |>
	flextable::compose(
		i = 9, j = 1,
		flextable::as_paragraph("2016 - 2017")
	) |>
	flextable::compose(
		i = 10, j = 1,
		flextable::as_paragraph("2017 - 2018")
	) |>
	flextable::compose(
		i = 11, j = 1,
		flextable::as_paragraph("2018 - 2019")
	) |>
	flextable::compose(
		i = 12, j = 1,
		flextable::as_paragraph("2019 - 2020")
	) |>
	flextable::compose(
		i = 13, j = 1,
		flextable::as_paragraph("2020 - 2021")
	) |>
	flextable::compose(
		i = 14, j = 1,
		flextable::as_paragraph("2021 - 2022")
	) |>
	# Next is the ages
	flextable::compose(
		i = 20, j = 1,
		flextable::as_paragraph("71 - 75")
	) |>
	flextable::compose(
		i = 21, j = 1,
		flextable::as_paragraph("76 - 80")
	) |>
	flextable::compose(
		i = 22, j = 1,
		flextable::as_paragraph("81 - 85")
	) |>
	# Next is the birth years
	flextable::compose(
		i = 25, j = 1,
		flextable::as_paragraph("1936 - 1940")
	) |>
	flextable::compose(
		i = 26, j = 1,
		flextable::as_paragraph("1941 - 1945")
	) |>
	flextable::compose(
		i = 27, j = 1,
		flextable::as_paragraph("1946 - 1950")
	) |>
	flextable::compose(
		i = 28, j = 1,
		flextable::as_paragraph("1951 - 1955")
	) |>
	flextable::compose(
		i = 29, j = 1,
		flextable::as_paragraph("1956 - 1960")
	) |>
	# Finally race and ethnicity, just one fix
	flextable::compose(
		i = 32, j = 1,
		flextable::as_paragraph("Black or Afircan American")
	) |>
	# Additionally set the font size since the default is big
	flextable::fontsize(part = "all", size = 10) |>
	flextable::set_table_properties(opts_pdf = list(tabcolsep = 0.75))

readr::write_rds(
	regression_tbl_dose_flextable,
	here::here("results", "tables", "dose-outcome-analysis-table.Rds")
)

# Table of number of repeats ####
# First count how many times each individual was observed
id_counts <-
	dat |>
	dplyr::distinct(study, id, season) |>
	dplyr::count(id, name = "visits")

# Now count how many times each individual/dose combo was observed
dose_counts <-
	dat |>
	dplyr::distinct(study, id, dose, season) |>
	dplyr::count(id, dose) |>
	# Add zero counts for people who never got one of the two treatments
	tidyr::complete(id, dose, fill = list(n = 0)) |>
	# Make one column for SD counts and one column for HD counts
	tidyr::pivot_wider(names_from = dose, values_from = n)

# Get the number of times someone switched doses
switch_counts <-
	dat |>
	dplyr::distinct(study, id, dose, season) |>
	dplyr::group_by(id) |>
	dplyr::mutate(last_dose = dplyr::lag(dose)) |>
	dplyr::ungroup() |>
	dplyr::mutate(
		dose_switch = dplyr::case_when(
			is.na(last_dose) ~ "no_last_dose",
			dose == "HD" & last_dose == "SD" ~ "HD_to_SD",
			dose == "SD" & last_dose == "HD" ~ "SD_to_HD",
			dose == last_dose ~ "same_dose",
			TRUE ~ NA_character_
		)
	) |>
	dplyr::count(id, dose_switch) |>
	# Add zero counts for people who never got one of the two treatments
	tidyr::complete(id, dose_switch, fill = list(n = 0)) |>
	# Make one column for SD counts and one column for HD counts
	tidyr::pivot_wider(names_from = dose_switch, values_from = n) |>
	# Select only the interesting columns
	dplyr::select(id, HD_to_SD, SD_to_HD)

# Combine the visit counts and the dose counts into one DF. This could be done
# with column binding instead but since the data frames are small the join is
# safer in case something has gone awry since this code was run last.
repeat_counts_df <-
	dplyr::left_join(
		id_counts,
		dose_counts,
		by = "id"
	) |>
	dplyr::left_join(
		switch_counts,
		by = "id"
	) |>
	dplyr::mutate(
		# Check to make sure the counting was consistent
		chk = (visits == (SD + HD)),
		# Recreate the study site variable -- this is the easiest way to get it
		# back rather than changing any of the counting commands
		study = id |>
			stringr::str_extract("^([a-z]{2,3})_.*", group = 1) |>
			factor(
				levels = c("uga", "fl", "pa"),
				labels = c("UGA", "FL", "PA")
			),
		# Force treatment of integer valued counts as factors
		dplyr::across(c(visits, SD, HD, SD_to_HD, HD_to_SD), as.ordered)
	)

# Check to make sure that sum of SD and HD equals number of repeats
if (isFALSE(all(repeat_counts_df$chk))) rlang::abort("SD + HD != Visits")

# Now make the table to show this info
repeat_counts_table <-
	repeat_counts_df |>
	gtsummary::tbl_summary(
		by = study,
		label = list(
			visits ~ "Number of total visits by same individual",
			SD ~ "Number of SD vaccinations for same individual",
			HD ~ "Number of HD vaccinations for same individual",
			HD_to_SD ~ "Number of times an individual switched from HD to SD",
			SD_to_HD ~ "Number of times an individual switched from SD to HD"
		),
		include = c(visits, SD, HD, HD_to_SD, SD_to_HD)
	) |>
	gtsummary::add_overall(last = TRUE)

repeat_counts_flextable <-
	repeat_counts_table |>
	gtsummary::as_flex_table() |>
	flextable::fontsize(size = 10, part = "all")

readr::write_rds(
	repeat_counts_flextable,
	here::here("results", "tables", "repeat-counts-table.Rds")
)

# END OF FILE ####
