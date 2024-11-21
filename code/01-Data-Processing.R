###
# Data Processing
# Zane Billings
# 2024-01-29
# Starting with the cleaned data from ahgroup/UGAFluVac repo, we need to do
# some data processing to get the data ready for our analyses.
# This script processes and filters the data.
###

# Setup ========================================================================
# Declare necessary dependencies
library(readr, include.only = NULL)
library(here, include.only = NULL)
library(dplyr, include.only = NULL)
library(forcats, include.only = NULL)

# Load the "raw" data
# This data is the output of the cleaning process from the repo
# https://github.com/ahgroup/UGAFluVac-data/, which includes the raw excel
# data files. We pulled the cleaned data from that repo on 2024-01-29.
raw_data <- readr::read_rds(here::here("data", "raw", "clean-data.Rds"))

# Data cleaning ================================================================
# First filter out the records we don't want
dat_filtered <-
	raw_data |>
	dplyr::filter(
		# Remove individuals with intradermal vaccines
		dose %in% c("SD", "HD"),
		# Remove individuals with infinite titerincrease (1 person who has a
		# pretiter of 0, probably a data entry error)
		is.finite(titerincrease),
		# Remove flu B assays
		strain_type %in% c("H1N1", "H3N2"),
		# Only individuals >= age 65 were allowed to get the HD vaccine, so we
		# want to filter out any records from individuals below that age.
		age >= 65
	) |>
	# Drop the missing factor levels that remain after filtering
	dplyr::mutate(
		# Change the sichuan 1989 (incorrect) records to 1987 so they are all
		# listed as the same strain, which they are, due to the insane data
		# cleaning error. SHould eventually be moved to the UGAFluVac repo.
		strains_fullname = forcats::fct_recode(
			strains_fullname,
			"H3N2-Sichuan-1987" = "H3N2-Sichuan-1989"
		),
		# Now drop the factor levels.
		dplyr::across(
			c(dose, strain_type, strains_fullname),
			forcats::fct_drop
		)
	)

# Now we need to make edits to certain variables and drop the variables that
# we don't need
dat_clean <-
	dat_filtered |>
	# First we need to "condense" the vaccine name variables. We just need one
	# vaccine name variable that has the name of the H1 vaccine strain for H1
	# assays, and the name of the H3 vaccine strain for H3 assays.
	dplyr::mutate(
		vaccine_name = dplyr::if_else(
			strain_type == "H1N1",
			h1n1_vaccine_fullname,
			h3n2_vaccine_fullname
		),
		.after = dose
	) |>
	# Now drop the *_vaccine_fullname columns -- we'll also drop out the
	# variables we don't plan to use at all in this analysis.
	dplyr::select(
		-dplyr::ends_with("vaccine_fullname"),
		-c(uniq_id, id, bmi, days_since_vac, date_vaccinated, race, gender),
		-dplyr::contains("particip")
	) |>
	# Next we need to clean up the date of birth variable. We just need a variable
	# for the birth year.
	dplyr::mutate(
		birth_year = substring(dateofbirth, 1, 4) |> as.integer(),
		.keep = "unused"
	) |>
	# Now we'll create a fold change variable which is 2 ^ (titerincrease)
	# which can make plotting easier
	dplyr::mutate(fold_change = 2 ^ titerincrease, .before = titerincrease) |>
	# Next let's transform the season variable into an ordered factor to ensure
	# it behaves the way we want during plotting. Since the levels are already
	# alphabetical it should be fine but better to go ahead and do it.
	dplyr::mutate(season = factor(season, ordered = TRUE)) |>
	# Finally let's clean up the names of the some of the variables.
	dplyr::rename(
		strain_name = strains_fullname,
		posttiter = postiter,
		log_pretiter = prevactiter,
		log_posttiter = postvactiter,
		titer_increase = titerincrease,
		id = subject_id
	)

# Now we need to make a few additional changes before we can pass the
# data to our models.
# First get the minimum birth year for making a centered/scaled/whatever version
# of the birth year variable.
MIN_BY <- raw_data$dateofbirth |> substring(1, 4) |> as.integer() |> min()
	
# Now we make those changes in the dataset.
dat_model <-
	dat_clean |>
	dplyr::mutate(
		# We need to make versions of the time and birth year variables that are
		# closer to being scale-free than the current versions -- models that
		# have those numbers in the thousands have worse conditioning problems that
		# can lead to numerical issues in an otherwise fine model. So we'll scale
		# the year variable by subtracting 2013, the first year, and we'll scale
		# the birth_year variable by subtracting the minimum birth year.
		year_c = year - 2013,
		birth_year_c = birth_year - MIN_BY,
		# We'll also center the age variable by subtracting 65.
		age_c = age - 65
	) |>
	# Finally let's reorder the variables. This is for no practical purpose
	# but it makes me happier.
	dplyr::select(
		id, study, season, year, year_c, age, age_c, birth_year, birth_year_c,
		dose, vaccine_name, strain_name, strain_type, pretiter, log_pretiter,
		posttiter, log_posttiter, fold_change, titer_increase, seroconversion,
		seroprotection
	)

# Save data to file ============================================================
readr::write_rds(
	dat_model,
	file = here::here("data", "processed", "model-data.Rds"),
	compress = "gz"
)
readr::write_csv(
	dat_model,
	file = here::here("data", "processed", "model-data.csv")
)

# ADDENDUM: Create CIVICs reporting data =======================================
rm(list = ls())
dat_used <- readr::read_rds(here::here("data", "processed", "model-data.Rds"))
dat <- readr::read_rds(here::here("data", "raw", "clean-data.Rds"))

# Filter the complete data to get just the people we used
dat_used$longid <- paste0(dat_used$id, dat_used$season)
dat$longid <- paste0(dat$subject_id, dat$season)

dat_is <-
	dat |>
	dplyr::filter(longid %in% dat_used$longid) |>
	dplyr::distinct(longid, .keep_all = TRUE)

dat_ages <- dat_is |>
	dplyr::group_by(subject_id) |>
	dplyr::mutate(
		Min_Age = min(age),
		Max_Age = max(age)
	) |>
	dplyr::ungroup()

dat_hsdi <-
	dat_is |>
	dplyr::transmute(
		Study_Code = "Study-248_HD_IIV",
		Subject_ID = subject_id,
		Cohort_ID = paste(study, year, dose, sep = "_"),
		Sex_Assigned_at_Birth = ifelse(
			is.na(gender), "Unknown", as.character(gender)
		),
		Gender = "Unknown",
		Min_Age = dat_ages$Min_Age,
		Max_Age = dat_ages$Max_Age,
		Subject_Age_Unit = "Years",
		Birth_Year = as.integer(substr(dateofbirth, 1, 4)),
		Subject_Age_Event = "Age at enrollment",
		Subject_Phenotype = "Not Collected",
		Subject_Location = ifelse(
			study == "UGA", "GA", study
		),
		Ethnicity = ifelse(
			race == "Hispanic",
			"Hispanic or Latino",
			"Not Hispanic or Latino"
		),
		Ethnicity = ifelse(
			is.na(Ethnicity),
			"Unknown",
			Ethnicity
		),
		Race = dplyr::case_match(
			race,
			"Black" ~ "Black or African American",
			"White" ~ "White",
			"Hispanic" ~ "Not Specified",
			"Other" ~ "OTH-Other",
			NA_character_ ~ "Unknown"
		),
		Subject_Description = "Not Provided"
	) |>
	dplyr::mutate(dplyr::across(dplyr::where(is.character), as.factor))

readr::write_rds(
	dat_hsdi,
	here::here("data", "processed", "reporting-data.Rds")
)

# END OF FILE ==================================================================
