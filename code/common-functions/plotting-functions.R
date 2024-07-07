###
# Plotting functions
# Zane
# 2024-07-06
# Many of the plots have to be repeated for all four outcomes so making the
# plot code a function means we don't have to copy and paste
###

# Helper function to get the correct x-axis label for the plot
get_x_axis_label_from_outcome <- function(
		outcome = c("Post-vaccination titer", "Seroconversion", "Seroprotection",
								"Titer increase")
	) {
	if (outcome == "Post-vaccination titer") {
		label <- "Exponentiated cACE\n(post-titer for HD / post-titer for SD)"
	} else if (outcome == "Seroconversion") {
		label <- "Exponentiated cACE\n(Seroprotection odds for HD / Seroprotection odds for SD)"
	} else if (outcome == "Seroprotection") {
		label <- "Exponentiated cACE\n(Seroconversion odds for HD / Seroconversion odds for SD)"
	} else if (outcome == "Titer increase") {
		label <- "Exponentiated cACE\n(fold change for HD / fold change for SD)"
	} else {
		rlang::abort("Wrong outcome provided! Check for typos.")
	}
	
	return(label)
}

# Function to save the plots to png and tiff format
save_plots <- function(
		plot_to_save,
		img_path = NULL, save_png = TRUE, save_tiff = TRUE,
		img_width = 6.5 * 2, img_height = 4 * 2
	) {
	# Optionally save as tiff
	if (isTRUE(save_tiff) & !is.null(img_path)) {
		ggplot2::ggsave(
			filename = paste0(img_path, ".tiff"),
			plot = plot_to_save,
			width = img_width,
			height = img_height
		)
	}
	
	# Optionally save as png
	if (isTRUE(save_png) & !is.null(img_path)) {
		ggplot2::ggsave(
			filename = paste0(img_path, ".png"),
			plot = plot_to_save,
			width = img_width,
			height = img_height
		)
	}
	
	invisible(plot_to_save)
}

# Plot for homologous cACEs, stratified by vaccine
vaccine_cace_plot <- function(
		outcome = c("Post-vaccination titer", "Seroconversion", "Seroprotection",
								"Titer increase"),
		strains_to_plot = c("all strains", "homologous", "heterologous"),
		...
	) {
	# First make the base plot
	vaccine_cate_nolines <-
		cate |>
		# Filter the cate dataset to get the right values
		dplyr::filter(
			model_label == outcome,
			data_label == strains_to_plot,
			cate_label %in% c("vaccine", "overall")
		) |>
		ggplot2::ggplot() +
		ggplot2::aes(
			x = y, xmin = ymin, xmax = ymax,
			y = forcats::fct_rev(vaccine_name)
		) +
		ggplot2::geom_vline(
			xintercept = 1, color = "gray50", linewidth = 1,
			linetype = "dashed", alpha = 0.8
		) +
		ggplot2::scale_y_discrete() +
		ggplot2::geom_hline(
			yintercept = 1:13,
			linetype = "dashed",
			linewidth = 0.25,
			color = "gray70"
		) +
		ggplot2::geom_linerange(
			ggplot2::aes(
				y = forcats::fct_rev(vaccine_name)
			),
			linewidth = 1.25
		) +
		ggplot2::geom_point(
			ggplot2::aes(
				y = forcats::fct_rev(vaccine_name),
				shape = subtype
			),
			size = 7,
			color = "black",
			show.legend = FALSE
		) +
		ggplot2::coord_cartesian(
			xlim = c(0.6, 1.52)
		) +
		ggplot2::scale_shape_manual(
			values = c(15, 17, 16)
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
			panel.spacing=unit(0,"cm"),
			strip.text = element_text(size = 32, margin = ggplot2::margin(r = 8)),
			axis.text.y = element_text(size = 24),
			axis.text.x = element_text(size = 24),
			axis.title.x = element_text(size = 32, margin = ggplot2::margin(t = 12))
		) +
		ggplot2::labs(
			x = get_x_axis_label_from_outcome(outcome),
			y = NULL
		)
	
	# Add the horizontal line annotations
	vaccine_cate <-
		vaccine_cate_nolines +
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
		)
	
	# Save the plots with the specifications and return the plot object
	invisible(save_plots(vaccine_cate, ...))
}

season_cace_plot <- function(
		outcome = c("Post-vaccination titer", "Seroconversion", "Seroprotection",
								"Titer increase"),
		strains_to_plot = c("all strains", "homologous", "heterologous"),
		...
	) {
	cate_season_only <-
		cate |>
		dplyr::filter(
			model_label == outcome,
			data_label == strains_to_plot,
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
		ggplot2::geom_linerange(
			linewidth = 1.25
		) +
		ggplot2::geom_point(
			size = 7,
			color = "black",
			show.legend = FALSE
		) +
		ggplot2::coord_cartesian(
			xlim = c(0.6, 1.52)
		) +
		ggplot2::theme(
			panel.grid.major.y = ggplot2::element_blank(),
			panel.grid.minor.y = ggplot2::element_blank(),
			strip.placement = "outside",
			strip.background = element_blank(),
			panel.spacing=unit(0,"cm"),
			strip.text = element_text(size = 32, margin = ggplot2::margin(r = 8)),
			axis.text.y = element_text(size = 24),
			axis.text.x = element_text(size = 24),
			axis.title.x = element_text(size = 32, margin = ggplot2::margin(t = 12))
		) +
		ggplot2::labs(
			x = get_x_axis_label_from_outcome(outcome),
			y = NULL
		)
	
	invisible(save_plots(cate_season_only, ...))
}

all_strains_cace_plot <- function(
		outcome = c("Post-vaccination titer", "Seroconversion", "Seroprotection",
								"Titer increase"),
		strains_to_plot = c("all strains", "homologous", "heterologous"),
		...
) {
	# Make the H1N1 subplot
	h1n1_all_strains <-
		cate |>
		dplyr::filter(
			model_label == outcome,
			data_label == strains_to_plot,
			cate_label %in% c("vaccine", "vaccine-strain"),
			subtype == "H1N1"
		) |>
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
		ggplot2::geom_linerange(
			linewidth = 0.75
		) +
		ggplot2::geom_point(
			size = 2.5,
			color = "black"
		) +
		ggplot2::coord_cartesian(
			xlim = c(0.6, 1.52)
		) +
		ggplot2::scale_color_manual(
			values = c("black", "gray50")
		) +
		ggplot2::scale_shape_manual(
			values = c(16, 17)
		) +
		ggplot2::facet_wrap(
			ggplot2::vars(vaccine_name),
			scales = "free_x",
			ncol = 3
		) +
		ggplot2::labs(
			x = NULL,
			y = NULL,
			color = NULL,
			shape = NULL,
			linetype = NULL
		) +
		ggplot2::theme(legend.position = "none")
	
	# Make the H3N2 subplot
	h3n2_all_strains <-
		cate |>
		dplyr::filter(
			model_label == outcome,
			data_label == strains_to_plot,
			cate_label %in% c("vaccine", "vaccine-strain"),
			subtype == "H3N2"
		) |>
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
		ggplot2::geom_linerange(
			linewidth = 0.75
		) +
		ggplot2::geom_point(
			size = 2.5,
			color = "black"
		) +
		ggplot2::coord_cartesian(
			xlim = c(0.6, 1.52)
		) +
		ggplot2::scale_color_manual(
			values = c("black", "gray50")
		) +
		ggplot2::scale_shape_manual(
			values = c(16, 17)
		) +
		ggplot2::facet_wrap(
			ggplot2::vars(vaccine_name),
			scales = "free_x",
			ncol = 4
		) +
		ggplot2::labs(
			x = get_x_axis_label_from_outcome(outcome),
			y = NULL,
			color = NULL,
			shape = NULL,
			linetype = NULL
		) +
		ggplot2::guides(
			shape = ggplot2::guide_legend(override.aes = list(size = 4))
		)
	
	# Combine the two subplots together and change theme elements
	fig1comb_noanno <-
		(h1n1_all_strains / h3n2_all_strains) +
		patchwork::plot_layout(heights = c(0.9, 1)) &
		theme(
			panel.grid.major.y = ggplot2::element_blank(),
			panel.grid.minor.y = ggplot2::element_blank(),
			legend.text = element_text(size = 28),
			axis.text.x = element_text(size = 20),
			axis.text.y = element_text(size = 12),
			strip.text = element_text(size = 32, margin = ggplot2::margin(b = 4)),
			axis.title.x = element_text(size = 32, margin = ggplot2::margin(t = 2)),
			panel.spacing=unit(0.2,"cm"),
			legend.key.size = unit(1, "cm"),
			legend.box.spacing = unit(0.1, "cm")
		)
	
	# Add the H1N1 and H3N2 labels
	fig1comb <-
		cowplot::ggdraw(fig1comb_noanno) +
		cowplot::draw_label(
			"H1N1:",
			x = 0.055, y = 0.981,
			size = 36,
			fontface = "bold"
		) +
		cowplot::draw_label(
			"H3N2:",
			x = 0.055, y = 0.55,
			size = 36,
			fontface = "bold"
		)
	
	# Save and return
	invisible(save_plots(fig1comb, ...))
}
