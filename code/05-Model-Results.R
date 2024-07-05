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
library(ggokabeito, include.only = NULL)
library(gtsummary, include.only = NULL)
library(gt, include.only = NULL)
library(hgp, include.only = NULL)
library(Hmisc, include.only = NULL)

# This script sets up convenient project paths
source(here::here("code", "common-functions", "path-setting.R"))

# Load the virus info for renaming strains
source(here::here("code", "common-functions", "virus-info.R"))

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

# ---- Raw data dose plot by strain ----
plt1_nested <-
	dat |>
	dplyr::select(dose, titer_increase, vaccine_name, strain_name, strain_type) |>
	tidyr::nest(dat = -strain_type)

plt1_data <-
	plt1_nested |>
	dplyr::mutate(
		plt = purrr::map(
			dat,
			\(d) ggplot(d) +
				aes(x = dose, y = titer_increase, color = dose, shape = dose) +
				geom_point(
					position = position_jitter(
						width = 0.25,
						height = 0.25,
						seed = 123124
					),
					size = 0.75,
					alpha = 0.25,
				) +
				coord_cartesian(ylim = c(-5, 10)) +
				labs(
					x = NULL,
					y = "Titer increase"
				) +
				facet_wrap(vars(vaccine_name), nrow = 2, scales = "free_x") +
				ggokabeito::scale_color_okabe_ito() +
				scale_shape_manual(values = c(1, 2)) +
				guides(
					color = guide_legend(override.aes = list(alpha = 1, size = 3))
				) +
				theme(
					strip.text = ggplot2::element_text(size = 14)
				)
		),
		plt = purrr::map2(
			plt, strain_type,
			\(p, st) p + ggtitle(st)
		)
	)

raw_dose_plt_vaccines <-
	plt1_data |>
	dplyr::pull(plt) |>
	purrr::reduce(`+`) +
	patchwork::plot_layout(
		guides = "collect",
		axes = "collect",
		axis_titles = "collect"
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "raw-dose-plot.tiff"),
	plot = raw_dose_plt_vaccines,
	width = 9.75,
	height = 9
)

# ---- Raw dose plot by vaccine/strain ----
plt2_summaries <-
	dat |>
	dplyr::select(dose, titer_increase, vaccine_name, strain_name, strain_type) |>
	# dplyr::mutate(
	# 	i = dplyr::row_number(),
	# 	.by = c(strain_type, vaccine_name, strain_name, dose)
	# ) |>
	# tidyr::pivot_wider(
	# 	names_from = dose,
	# 	values_from = titer_increase
	# ) |>
	dplyr::summarise(
		ggplot2::mean_cl_boot(titer_increase),
		.by = c(strain_type, vaccine_name, strain_name, dose)
	)

test <-
	plt2_summaries |>
	dplyr::filter(strain_type == "H1N1") |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(strain_name),
		color = dose
	) +
	ggplot2::geom_vline(xintercept = 0, color = "gray50", linewidth = 1) +
	# ggplot2::scale_y_discrete() +
	# ggplot2::geom_hline(
	# 	yintercept = 1:13,
	# 	linetype = "dashed",
	# 	linewidth = 0.25,
	# 	color = "gray70"
	# ) +
	ggplot2::geom_pointrange() +
	ggplot2::facet_wrap(
		facets = ggplot2::vars(vaccine_name),
		scales = "free_x"
	) +
	ggplot2::labs(
		x = "Expected change in titer increase (titer units)",
		y = NULL
	)

# ---- MAIN RESULTS FOR TITER INCREASE ----
# homologous strains only
homologous_vaccine_cate <-
	cate |>
	dplyr::filter(
		model_label == "Titer increase",
		data_label == "homologous",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (fold change for HD / fold change for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate.tiff"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate.png"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# heterologous stratified by vaccine
{
heterologous_vaccine_cate_nolines <-
	cate |>
	dplyr::filter(
		model_label == "Titer increase",
		data_label == "all strains",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (fold change for HD / fold change for SD)",
		y = NULL
	)

heterologous_vaccine_cate <-
	heterologous_vaccine_cate_nolines +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 1
	) +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 2
	) +
	theme(
		axis.text.x = element_text(size = 14),
		axis.text.y = element_text(size = 12)
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate.tiff"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate.png"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
}

# Stratified by season
cate_season_only <-
	cate |>
	dplyr::filter(
		model_label == "Titer increase",
		data_label == "all strains",
		cate_label %in% c("season", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(season)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		axis.text.x = element_text(size = 18),
		axis.text.y = element_text(size = 16)
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (fold change for HD / fold change for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate.tiff"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate.png"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)

# cACE plots for all assay strains ####
# First we need to extract the overall data
all_strains_cate_data <-
	dplyr::bind_rows(
		# This part has the "Overall" row for each vaccine
		cate |>
			dplyr::filter(
				model_label == "Titer increase",
				data_label == "all strains",
				cate_label %in% c("vaccine")
			),
		# This part has all of the assay-strain_specific CATES
		cate |>
			dplyr::filter(
				model_label == "Titer increase",
				data_label == "all strains",
				cate_label %in% c("vaccine-strain")
			)
	)

# H1N1 assay strains stratified
h1n1_all_strains <-
	all_strains_cate_data |>
	dplyr::filter(subtype == "H1N1") |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom, shape = hom, linetype = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dotted",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (fold change for HD / fold change for SD)",
		y = NULL,
		color = NULL,
		shape = NULL,
		linetype = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate.tiff"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate.png"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

# H3N2 assay strains stratified
h3n2_all_strains <-
	all_strains_cate_data |>
	dplyr::filter(subtype == "H3N2") |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom, shape = hom, linetype = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dotted",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (fold change for HD / fold change for SD)",
		y = NULL,
		color = NULL,
		shape = NULL,
		linetype = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate.tiff"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate.png"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

# Make combined figure
{
p1 <- h1n1_all_strains
p2 <- h3n2_all_strains
fig1comb <-
	p1 / p2 +
	patchwork::plot_annotation(tag_levels = "A", tag_suffix = ")") &
	theme(
		plot.tag = element_text(size = 36),
		strip.text = element_text(size = 24),
		axis.text.x = element_text(size = 16),
		axis.text.y = element_text(size = 10),
		axis.title = element_text(size = 20),
		legend.text = element_text(size = 20)
	)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "all-strains-cate.tiff"),
	plot = fig1comb,
	width = 6.5 * 2.25,
	height = 8 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "all-strains-cate.png"),
	plot = fig1comb,
	width = 6.5 * 2.25,
	height = 8 * 2.25
)
}

# Supplement results ===========================================================

## Post-vaccination titer ======================================================
# Homologous plot
homologous_vaccine_cate <-
	cate |>
	dplyr::filter(
		model_label == "Post-vaccination titer",
		data_label == "homologous",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (post-titer for HD / post-titer for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-post.tiff"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-post.png"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# Stratified by vaccine
heterologous_vaccine_cate_nolines <-
	cate |>
	dplyr::filter(
		model_label == "Post-vaccination titer",
		data_label == "all strains",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (post-titer for HD / post-titer for SD)",
		y = NULL
	)

heterologous_vaccine_cate <-
	heterologous_vaccine_cate_nolines +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 1
	) +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 2
	) +
	theme(
		axis.text.x = element_text(size = 14),
		axis.text.y = element_text(size = 12)
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-post.tiff"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-post.png"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# Stratified by season
cate_season_only <-
	cate |>
	dplyr::filter(
		model_label == "Post-vaccination titer",
		data_label == "all strains",
		cate_label %in% c("season", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(season)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		axis.text.x = element_text(size = 18),
		axis.text.y = element_text(size = 16)
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (post-titer for HD / post-titer for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-post.tiff"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-post.png"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)

# ---- H1N1 strain/vaccine cATE plot ----
h1n1_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Post-vaccination titer",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H1N1"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (post-titer for HD / post-titer for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-post.tiff"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-post.png"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

# ---- H3N2 strain/vaccine cATEs plot ----
h3n2_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Post-vaccination titer",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H3N2"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (post-titer for HD / post-titer for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-post.tiff"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-post.png"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

## Seroprotection ==============================================================
homologous_vaccine_cate <-
	cate |>
	dplyr::filter(
		model_label == "Seroprotection",
		data_label == "homologous",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroprotection odds for HD / Seroprotection odds for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-sp.tiff"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-sp.png"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# ---- Heterologous cATEs for each vaccine ----
heterologous_vaccine_cate_nolines <-
	cate |>
	dplyr::filter(
		model_label == "Seroprotection",
		data_label == "all strains",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroprotection odds for HD / Seroprotection odds for SD)",
		y = NULL
	)

heterologous_vaccine_cate <-
	heterologous_vaccine_cate_nolines +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 1
	) +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 2
	) +
	theme(
		axis.text.x = element_text(size = 14),
		axis.text.y = element_text(size = 12)
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-sp.tiff"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-sp.png"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# Titer increase by season
cate_season_only <-
	cate |>
	dplyr::filter(
		model_label == "Seroprotection",
		data_label == "all strains",
		cate_label %in% c("season", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(season)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		axis.text.x = element_text(size = 18),
		axis.text.y = element_text(size = 16)
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroprotection odds for HD / Seroprotection odds for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-sp.tiff"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-sp.png"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)

# ---- H1N1 strain/vaccine cATE plot ----
h1n1_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Seroprotection",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H1N1"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroprotection odds for HD / Seroprotection odds for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-sp.tiff"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-sp.png"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

# ---- H3N2 strain/vaccine cATEs plot ----
h3n2_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Seroprotection",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H3N2"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroprotection odds for HD / Seroprotection odds for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-sp.tiff"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-sp.png"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

## Seroconversion ==============================================================
homologous_vaccine_cate <-
	cate |>
	dplyr::filter(
		model_label == "Seroconversion",
		data_label == "homologous",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroconversion odds for HD / Seroconversion odds for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-sc.tiff"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "homologous-vaccine-cate-sc.png"),
	plot = homologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# ---- Heterologous cATEs for each vaccine ----
heterologous_vaccine_cate_nolines <-
	cate |>
	dplyr::filter(
		model_label == "Seroconversion",
		data_label == "all strains",
		cate_label %in% c("vaccine", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(vaccine_name)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange(aes(y = forcats::fct_rev(vaccine_name))) +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::facet_grid(
		rows = ggplot2::vars(subtype),
		space = "free_y",
		scales = "free_y",
		switch = "y"
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		strip.placement = "outside",
		strip.background = element_blank(),
		panel.spacing=unit(0,"cm")
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroconversion odds for HD / Seroconversion odds for SD)",
		y = NULL
	)

heterologous_vaccine_cate <-
	heterologous_vaccine_cate_nolines +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 1
	) +
	ggh4x::at_panel(
		annotate(
			"segment", y = 0.25, yend = 0.25, x = -Inf, xend = Inf,
			color = "black", alpha = 0.25, linewidth = 2
		),
		PANEL == 2
	) +
	theme(
		axis.text.x = element_text(size = 14),
		axis.text.y = element_text(size = 12)
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-sc.tiff"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "heterologous-vaccine-cate-sc.png"),
	plot = heterologous_vaccine_cate,
	width = 6.5 * 2,
	height = 4 * 2
)

# Titer increase by season
cate_season_only <-
	cate |>
	dplyr::filter(
		model_label == "Seroconversion",
		data_label == "all strains",
		cate_label %in% c("season", "overall")
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax,
		y = forcats::fct_rev(season)
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:13,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::scale_x_continuous(
		limits = c(0.5, 1.5)
	) +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
		axis.text.x = element_text(size = 18),
		axis.text.y = element_text(size = 16)
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroconversion odds for HD / Seroconversion odds for SD)",
		y = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-sc.tiff"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "season-only-cate-sc.png"),
	plot = cate_season_only,
	width = 6.5 * 2,
	height = 4 * 2
)

# ---- H1N1 strain/vaccine cATE plot ----
h1n1_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Seroconversion",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H1N1"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroconversion odds for HD / Seroconversion odds for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-sc.tiff"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h1n1-all-strains-cate-sc.png"),
	plot = h1n1_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)


# ---- H3N2 strain/vaccine cATEs plot ----
h3n2_all_strains <-
	cate |>
	dplyr::filter(
		model_label == "Seroconversion",
		data_label == "all strains",
		cate_label %in% c("vaccine-strain"),
		subtype == "H3N2"
	) |>
	ggplot2::ggplot() +
	ggplot2::aes(
		x = y, xmin = ymin, xmax = ymax, color = hom,
		y = forcats::fct_rev(strain_name),
	) +
	ggplot2::geom_vline(xintercept = 1, color = "gray50", linewidth = 1) +
	ggplot2::scale_y_discrete() +
	ggplot2::geom_hline(
		yintercept = 1:100,
		linetype = "dashed",
		linewidth = 0.25,
		color = "gray70"
	) +
	ggplot2::geom_pointrange() +
	ggplot2::coord_cartesian(
		xlim = c(0.5, 1.5)
	) +
	ggokabeito::scale_color_okabe_ito() +
	ggplot2::theme(
		panel.grid.major.y = ggplot2::element_blank(),
		panel.grid.minor.y = ggplot2::element_blank(),
	) +
	ggplot2::facet_wrap(
		ggplot2::vars(vaccine_name),
		scales = "free_x",
		nrow = 2
	) +
	ggplot2::labs(
		x = "Exponentiated cACE (Seroconversion odds for HD / Seroconversion odds for SD)",
		y = NULL,
		color = NULL
	)

ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-sc.tiff"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)
ggplot2::ggsave(
	filename = here::here("Results", "Figures", "h3n2-all-strains-cate-sc.png"),
	plot = h3n2_all_strains,
	width = 6.5 * 2.25,
	height = 4 * 2.25
)

# END OF FILE ==================================================================
