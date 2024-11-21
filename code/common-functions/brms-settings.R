###
# brms model settings
# Yang Ge, edits by Zane
# 2024-01-29
# We want to use the same control parameters and priors across all of the
# models with different outcomes, so we specify them all in one place.
###

# We need to attach brms to make these things load correctly.
require(brms, quietly = TRUE)
# Used for determining number of cores for parallel
library(parallelly, include.only = NULL)

# MCMC control parameters ======================================================

BRMS_SEED <- 28127343
# Determine the number of cores for parallel processing
# This avoids many common issues such as NA, 0, or all cores being passed.
# THIS WILL LEAVE TWO OF YOUR CORES OPEN.
# You should ALWAYS leave at least one core open or you will experience
# substantial slow down. Leaving a second core open makes it easier for you
# to do something like look at emails while the code runs.
available_cores <-
	parallelly::availableCores(
		constraints = "connections",
		omit = 2
	)
# 500 warmup iterations/core should be sufficient for this model, we had a
# target of 20k total iterations though. So we need to set the number of
# iterations per chain correctly to match this target.
BRMS_WARMUPS <- 500
# The maximum configuration is 16 cores, which runs 1250 sampling iterations
# per chain (one chain is run per core). Since sampling iterations are typically
# a good bit faster than warmup iterations, it doesn't usually result in a
# noticeable increase to run across more cores, since they all need a
# non-trivial of warmup before sampling can begin. If you want to use more
# cores you can edit this code.
if (available_cores >= 16) {
	BRMS_CORES <- 16
	BRMS_CHAINS <- 16
	BRMS_SAMPLING <- 1250
} else if (available_cores >= 8) {
	# Step down from 16 to 8 cores. You could of course run 12 or any number in
	# the middle, but this was arbitrary to avoid integer division of the number
	# of sampling iterations.
	BRMS_CORES <- 8
	BRMS_CHAINS <- 8
	BRMS_SAMPLING <- 2500
} else if (available_cores >= 4) {
	# Next stepdown is 4 cores, which is most common.
	BRMS_CORES <- 4
	BRMS_CHAINS <- 4
	BRMS_SAMPLING <- 5000
} else {
	# In the case where < 4 cores are available, we'll run sequentially. However,
	# we still need to run multiple chains to ensure convergence. So we'll run
	# 4 chains in sequence (this will take a LONG time).
	BRMS_CORES <- 1
	BRMS_CHAINS <- 4
	BRMS_SAMPLING <- 5000
}

BRMS_ITERS <- BRMS_SAMPLING + BRMS_WARMUPS
BRMS_ADAPT_DELTA <- 0.99
BRMS_MAX_TREEDEPTH <- 15

# priors for mixed effect models ===============================================
priors_me_gaussian_outcome <- c(
	prior(normal(0,5), class = "Intercept"),
	prior(normal(0,5), class = "b"),
	prior(student_t(3, 0, 1), class = "sd", lb = 0),
	prior(lkj_corr_cholesky(2), class = "cor"),
	prior(student_t(3, 0, 3), class = "sigma", lb = 0),
	prior(student_t(3, 0, 3), class = "sds", lb = 0)
)

# The binomial models don't have a sigma parameter, which will make brms mad
# at us if we try to pass one.
priors_me_binomial_outcome <- c(
	prior(normal(0,1), class = "Intercept"),
	prior(normal(0,1), class = "b"),
	prior(lkj_corr_cholesky(2), class = "cor"),
	prior(student_t(3, 0, 1), class = "sd", lb = 0),
	prior(student_t(3, 0, 1), class = "sds", lb = 0)
)


# priors for glm (fixed effects) models=========================================
priors_fe_gaussian_outcome <- c(
	prior(normal(0,5), class = "Intercept"),
	prior(normal(0,5), class = "b"),
	prior(student_t(3, 0, 1), class = "sigma", lb = 0)
)

priors_fe_binomial_outcome <- c(
	prior(normal(0,1), class = "Intercept"),
	prior(normal(0,1), class = "b"),
	prior(student_t(3, 0, 1), class = "sd", lb = 0)
)

# End of File ==================================================================
