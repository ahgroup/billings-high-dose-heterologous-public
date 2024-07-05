###
# convenient file path variables for saving plots and rds
# Yang Ge, edits by Zane
# Define all the paths in one place so we only have to do it once
###

library(purrr, include.only = NULL)

project_path <- here::here()
res_path <- here::here(project_path, "results")
rds_path <- here::here(project_path, "results", "files")
figure_path <- here::here(project_path, "results", "figures")
largefile_path <- here::here(project_path, "results", "largefiles")

# Create those paths if they don't exist (they may not depending on how the
# repo is set up)
purrr::walk(
	list(res_path, rds_path, figure_path, largefile_path),
	\(f) if (!dir.exists(f)) {dir.create(f)}
)
