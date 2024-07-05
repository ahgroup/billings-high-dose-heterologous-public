###
# Virus info and name replacer
# Zane
# 2024-04-09
# Loads the virus info database along with a function that allows for
# replacing virus names with one of the other forms.
# Needed to be in two places, so refactored it into this file.
###

# Declare dependencies ----
requireNamespace("here")
requireNamespace("readr")
requireNamespace("dplyr")
requireNamespace("tibble")
requireNamespace("forcats")

# Load the virus info data ----
virus_info <- readr::read_csv(
	here::here("data", "UGAFluVac-virus-names.csv"),
	col_types = "fcccil"
) |>
	# Remove the useless columns
	dplyr::select(-c(vaccine_strain)) |>
	# Append a row so sorting the overall entry for CATEs is easy
	tibble::add_row(
		subtype = "",
		analysis_name = "Overall",
		short_name = "Overall",
		genbank_strain_name = "Overall",
		factor_order = 9999L
	) |>
	# Make all of the name variables ordered factors and clean up the subtypes
	dplyr::mutate(
		subtype = factor(
			as.character(subtype),
			levels = c("h1", "h3", ""),
			labels = c("H1N1", "H3N2", "")
		),
		# Put the different name factors in order
		dplyr::across(
			c(analysis_name, genbank_strain_name, short_name),
			\(x) forcats::fct_reorder(x, factor_order)
		),
	)

# Function for replacing strain names ----
replace_strain_names <- function(x, from = "analysis", to = "short",
																 drop = TRUE) {
	# Find the right column for selecting names from
	if (from == "analysis") {
		from_vec <- virus_info$analysis_name
	} else if (from == "full") {
		from_vec <- virus_info$genbank_strain_name
	} else if (from == "short") {
		from_vec <- virus_info$short_name
	} else {
		stop("'from' should be 'analysis', 'full', or 'short'.")
	}
	
	# Make sure all values of x exist in the virus info table
	if (!(all(x %in% from_vec))) {
		stop(paste0(
			"'x' should be a vector of ", from, " names that exist in the",
			'virus-info sheet.'
			))
	}
	
	# Now get the location in the virus info table for each element of x
	locs <- match(x, from_vec)
	
	# Based on the names argument, get the correct names to return.
	if (to == "analysis") {
		vals <- virus_info$analysis_name[locs]
	} else if (to == "full") {
		vals <- virus_info$genbank_strain_name[locs]
	} else if (to == "short") {
		vals <- virus_info$short_name[locs]
	} else if (to == "subtype") {
		vals <- virus_info$subtype[locs]
	}
	else {
		stop("'to' should be 'analysis', 'full', 'short', or 'subtype'.")
	}
	
	# If requested, remove unseen factor levels
	if (isTRUE(drop)) {
		vals <- forcats::fct_drop(vals)
	}
	
	return(vals)
}

# End of file ----
