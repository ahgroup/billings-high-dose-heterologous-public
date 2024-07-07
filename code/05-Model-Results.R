#
# Model Results: Figures and Maybe Tables
# Zane Billings
# 2024-02-01
# This script takes the CATEs calculated from the previous script and makes the
# figures and tables we need for the manuscript and supplement
#
# ---- Setup ----

## Package reqs ####
library(here, include.only = NULL)
library(readr, include.only = NULL)
library(dplyr, include.only = NULL)
library(ggplot2)
library(patchwork)
library(cowplot, include.only = NULL)
library(gtsummary, include.only = NULL)
library(gt, include.only = NULL)
library(hgp, include.only = NULL)
library(Hmisc, include.only = NULL)

# This script sets up convenient project paths
source(here::here("code", "common-functions", "path-setting.R"))

# Load the virus info for renaming strains
source(here::here("code", "common-functions", "virus-info.R"))

# Script with plotting functions
source(here::here("code", "common-functions", "plotting-functions.R"))

## Set the ggplot theme and gtsummary theme ####
ggplot2::theme_set(hgp::theme_ms())

# ---- Data Processing ----

## Data loading ####
# Load the main dataset
dat <- readr::read_rds(here::here("data", "processed", "model-data.Rds")) |>
	# Replace the vaccine and strain names with short form
	dplyr::mutate(
		dplyr::across(
			c(vaccine_name, strain_name),
			replace_strain_names
		)
	)

# Load the calculated CATEs
raw_cate <- readr::read_rds(
	here::here(res_path, "files", "all-cates-combined.Rds")
)

# Clean up the CATE data and format nicely for plots
cate_formatted <-
	raw_cate |>
	# Remove the Sich/89 erroneous entries
	# TODO REMOVE THIS WHEN THE MODELS ARE RERUN
	dplyr::filter(strain_name != "H3N2-Sichuan-1989") |>
	# Format the strain and vaccine names nicely
	dplyr::mutate(
		# Capitalize all the "overall" entries in the strain name
		dplyr::across(
			c(strain_name, vaccine_name),
			\(x) {
				y <- as.character(x)
				out <- ifelse(y == "overall", "Overall", y)
				return(out)
			}
		),
		# Make a variable with the subtypes for plotting
		subtype = replace_strain_names(vaccine_name, to = "subtype"),
		# Replace the vaccine and strain names with short form
		dplyr::across(
			c(vaccine_name, strain_name),
			replace_strain_names
		)
	) |>
	# Format the different factors nicely
	dplyr::mutate(
		data_label = forcats::fct_recode(data_label, "all strains" = "all"),
		model_label = factor(
			as.character(model_label),
			levels = c("pt", "ti", "sp", "sc"),
			labels = c(
				"Post-vaccination titer", "Titer increase",
				"Seroprotection", "Seroconversion"
			)
		),
		# Homologous indicator for plot colors
		hom = factor(
			as.character(vaccine_name) == as.character(strain_name),
			levels = c(FALSE, TRUE),
			labels = c("Heterologous", "Homologous")
		)
	) |>
	dplyr::select(-.rows)

## ---- Cate Transformation ----
# Define a function to transform the CATEs into a more understandable unit
# This is a function so it is easy to modify/add arguments if Andreas changes
# his mind about what units make the most sense.
# N.b. the transformation must be monotonic in order to preserve the properties
# of the credible interval (ymin, ymax).
transform_cates <- function(x, method = "fc") {
	if (method == "fc") {
		# Fold change method -- since the outcome is TI, just do 2 ^ value.
		out <- 2^x
	} else if (method == "ih") {
		# ih = Inverse HAI function transformation. Before modeling, we did the
		# transformation y = g(z) = log2(z / 5). This is the inverse transformation
		# g^{-1}(y) = 5 * (2 ^ y).
		out <- 5 * (2^y)
	}
	
	return(out)
}

# Transform the raw CATE outcomes into a more understandable unit
cate <-
	cate_formatted |>
	dplyr::mutate(
		dplyr::across(
			c(y, ymin, ymax),
			\(x) transform_cates(x, method = "fc")
		)
	)

readr::write_rds(cate, here::here("results", "files", "cate-formatted.Rds"))

# ---- MAIN RESULTS FOR TITER INCREASE ----
# homologous strains only
homologous_vaccine_cace_ti <-
	vaccine_cace_plot(
		outcome = "Titer increase",
		strains_to_plot = "homologous",
		img_path = here::here("Results", "Figures", "homologous-vaccine-cate")
	)

# heterologous stratified by vaccine
heterologous_vaccine_cace_ti <-
	vaccine_cace_plot(
		outcome = "Titer increase",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "heterologous-vaccine-cate")
	)

# Stratified by season
season_cace_ti <-
	season_cace_plot(
		outcome = "Titer increase",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "season-only-cate")
	)

# cACE plots for all assay strains ####
strains_cace_ti <-
	all_strains_cace_plot(
		outcome = "Titer increase",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "all-strains-cate"),
		img_width = 6.5 * 2.25,
		img_height = 8 * 2.25
	)

# Supplement results ===========================================================

## Post-vaccination titer ======================================================
# homologous only stratified by vaccine
homologous_vaccine_cace_pt <-
	vaccine_cace_plot(
		outcome = "Post-vaccination titer",
		strains_to_plot = "homologous",
		img_path = here::here("Results", "Figures", "homologous-vaccine-cate-post")
	)

# heterologous stratified by vaccine
heterologous_vaccine_cace_pt <-
	vaccine_cace_plot(
		outcome = "Post-vaccination titer",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "heterologous-vaccine-cate-post")
	)

# Stratified by season
season_cace_pt <-
	season_cace_plot(
		outcome = "Post-vaccination titer",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "season-only-cate-post")
	)

# cACE plots for all assay strains ####
strains_cace_pt <-
	all_strains_cace_plot(
		outcome = "Post-vaccination titer",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "all-strains-cate-post"),
		img_width = 6.5 * 2.25,
		img_height = 8 * 2.25
	)

## Seroprotection ==============================================================
# homologous only stratified by vaccine
homologous_vaccine_cace_sp <-
	vaccine_cace_plot(
		outcome = "Seroprotection",
		strains_to_plot = "homologous",
		img_path = here::here("Results", "Figures", "homologous-vaccine-cate-sp")
	)

# heterologous stratified by vaccine
heterologous_vaccine_cace_sp <-
	vaccine_cace_plot(
		outcome = "Seroprotection",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "heterologous-vaccine-cate-sp")
	)

# Stratified by season
season_cace_sp <-
	season_cace_plot(
		outcome = "Seroprotection",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "season-only-cate-sp")
	)

# cACE plots for all assay strains ####
strains_cace_sp <-
	all_strains_cace_plot(
		outcome = "Seroprotection",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "all-strains-cate-sp"),
		img_width = 6.5 * 2.25,
		img_height = 8 * 2.25
	)

## Seroconversion ==============================================================
# homologous only stratified by vaccine
homologous_vaccine_cace_sc <-
	vaccine_cace_plot(
		outcome = "Seroconversion",
		strains_to_plot = "homologous",
		img_path = here::here("Results", "Figures", "homologous-vaccine-cate-sc")
	)

# heterologous stratified by vaccine
heterologous_vaccine_cace_sc <-
	vaccine_cace_plot(
		outcome = "Seroconversion",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "heterologous-vaccine-cate-sc")
	)

# Stratified by season
season_cace_sc <-
	season_cace_plot(
		outcome = "Seroconversion",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "season-only-cate-sc")
	)

# cACE plots for all assay strains ####
strains_cace_sc <-
	all_strains_cace_plot(
		outcome = "Seroconversion",
		strains_to_plot = "all strains",
		img_path = here::here("Results", "Figures", "all-strains-cate-sc"),
		img_width = 6.5 * 2.25,
		img_height = 8 * 2.25
	)

# END OF FILE ==================================================================
