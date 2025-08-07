rm(list=ls())

# remotes::install_github("ampl-psych/EMC2",ref="ss-plots")
library(EMC2)
#setwd("MBNCSS/Day4-AdvancedCognition")
# ------------------------------------------------------------------------------
# Script: BEESTS Stop-Signal Task Analysis using EMC2
#
# Overview:
# This script demonstrates how to analyze stop-signal task data using the
# Bayesian Estimation of Ex-Gaussian Stop-signal Reaction Time distributions 
# (BEESTS) implemented in the EMC2 package.
#
# Sections:
# - Simulate single-subject stop-signal task data
# - Fit single-subject model
# - Perform posterior predictive checks and SSRT estimation
# - Analyze real data with hierarchical modeling
#
# NOTE: BEESTS functionality is a recent addition to EMC2. It is still experimental,
# undergoing robustness checks and validation, and not yet available
# on CRAN. Until officially released on CRAN, please contact the package authors
# about research applications.
#
# ------------------------------------------------------------------------------

# Load pre-computed analysis results to avoid re-running time-consuming sampling
print(load("objects_ss.RData"))

# ------------------------------------------------------------------------------
# Simulate Single-Subject Stop-Signal Data
# ------------------------------------------------------------------------------

# Simulate stop-signal task with a simple left/right go response
ADmat <- cbind(d = c(-0.5, 0.5))

# Define experimental design
designSS <- design(
  model = SSexG,
  factors = list(subjects = 1, S = c("left", "right")),
  Rlevels = c("left", "right"),
  matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
  contrasts = list(lM = ADmat),
  formula = list(
    mu ~ lM, sigma ~ 1, tau ~ 1,
    muS ~ 1, sigmaS ~ 1, tauS ~ 1,
    gf ~ 1, tf ~ 1
  )
)

# mu is the intercept
# mu_lMd reflects the difference between the matching and mismatching go runners.


# Model parameters:
# - Go process (ex-Gaussian): mu, sigma, tau
# - Stop process (ex-Gaussian): muS, sigmaS, tauS
# - Go failure probability: gf (probit scale)
# - Stop trigger failure probability: tf (probit scale)

# - mu, sigma, tau, muS, sigmaS, tauS: require log transformation
# - gf, tf: require probit transformation
?SSexG

# Simulate parameter vector 
p_vector <- sampled_pars(designSS, doMap = FALSE)
p_vector


p_vector <- c(
  # mu, mu_lMTRUE, sigma, tau
  log(0.6), log(0.8), log(0.06), log(0.3),
  # muS, sigmaS, tauS
  log(0.25), log(0.02), log(0.05),
  # gf, tf
  qnorm(0.15), qnorm(0.2)
)

# It is good practice to use mapped_pars to check if the specified design and parameters 
# are sensible given the model. For race models to produce higher than chance accuracy, 
# we need the matching runner to be faster than the mismatching runner, 
# which is the case in this design: mu is smaller for
# matching than for mismatching
mapped_pars(designSS,p_vector)

# Define stop-signal delay (SSD) staircase parameters

RNGkind("L'Ecuyer-CMRG")
set.seed(4293)

staircase <- list(SSD0 = 0.25, stairstep = 0.05, stairmin = 0, stairmax = Inf)


# Simulate data
#dat <- make_data(
#   p_vector, designSS, n_trials = 400, staircase = staircase,
#   functions = list(SSD = function(d) EMC2:::SSD_function(d, SSD = NA, p = 0.25)))


# Required columns in the data: S, R, SSD, rt, subjects
# - SSD must be numeric (Inf on go-trials) and in seconds
# - rt is NA for non-responses, and in seconds
head(dat)
str(dat)
# right-skewed RT distributions
plot_density(dat, factors = "S")

# Separate stop and go
plot_density(dat, factors = c("SS","S"), functions=list(SS=\(d) is.finite(d$SSD)))
# Here we can see go-RTs and stop-signal response RTs (SRRT)
# The model predicts that SRRTs are faster than Go-RTs, as is the case here in this
# simulated dataset.


# Plot staircase behavior across stop trials
layout(1)
stop_dat <- dat[is.finite(dat$SSD), ]
go_dat <- dat[!is.finite(dat$SSD), ]
plot(stop_dat$trials, stop_dat$SSD, ylim = c(0, 0.7),
     xlab = "Stop trials over time", ylab = "SSD")

# Staircase performance: proportion of successful inhibitions (NA = inhibition)
# at around 50% as intended
mean(is.na(stop_dat$rt))

# Go-trial accuracy
# around 70% correct responses, 30% errors
mean(go_dat$R == go_dat$S, na.rm = TRUE)

# Go-omission rate (should reflect gf)
mean(is.na(go_dat$rt))
# meaning there are around 11% of non-responses on go-trials

# ------------------------------------------------------------------------------
# Model Assumption Checks
# ------------------------------------------------------------------------------

# Context independence assumption: Go RT > Stop-Signal Response RT (SRRT)
mean(go_dat$rt, na.rm = TRUE) > mean(stop_dat$rt, na.rm = TRUE)

# Visual check
plot_density(dat, defective_factor = "SS", functions=list(SS=\(d) is.finite(d$SSD)))

# Inhibition function: response probability should increase with SSD
plot_ss_if(dat, use_global_quantiles = FALSE, probs = seq(0, 1, length.out = 5))
# Per default, the SSD-categories are defined in terms of the percentiles of the 
# SSD distribution for each participant, and then averaged over participants.
# If no posterior predictive data are supplied, the data is plotted with error bars 
# (plus/minus the standard error per SSD bin/category) 
?plot_ss_if

# Signal-response RT (SRRT) should increase with SSD
plot_ss_srrt(dat, use_global_quantiles = TRUE, probs = seq(0, 1, length.out = 5),
             ylim = c(0.6, 0.75))
# Per default, SSDs are pooled over participants before calculating percentiles, 
# so the same absolute SSD range is used to get mean SRRT for each participant, 
# and then these probabilities are averaged over participants (see use_global_quantiles).
?plot_ss_srrt

# ------------------------------------------------------------------------------
# Single-Level Model Fitting
# ------------------------------------------------------------------------------

# Define prior distributions for parameters
prior_SS <- prior(
  designSS, type = "single",
  pmean = c(log(0.6), 0, log(0.5), log(0.5), log(0.7), log(0.5), log(0.5), qnorm(0.5), qnorm(0.5)),
  psd = c(rep(1, 7), c(1, 1))
)
plot(prior_SS)

# Create EMC2 object for fitting (It still needs to be checked how compression affects recovery for this model, which
# is why we turn it off for estimation.)
emc_obj <- make_emc(data = dat, design = designSS, compress = FALSE,
                    rt_resolution = 1e-6, type = "single", prior_list = prior_SS)
# single_ss <- fit(emc_obj, fileName = "single_ss.RData")
# save(single_ss,file="single_ss.RData")

# Check convergence diagnostics (Rhat, ESS)
check(single_ss)
plot(single_ss)  # Trace plots

# Visualize posterior vs prior; check whether priors are updated
plot_pars(single_ss)
plot_pars(single_ss, use_prior_lim = FALSE)
# All priors have been updated. Note that the go-related parameters including gf are
# better updated than the stop-related parameters. The main reason for this is
# that there is less information (only 25% of the trials are stop trials),
# and the fact that the finishing time distribution of the stop runner must be 
# inferred indirectly via response rates only, in the absence of observable RTs.

# Posterior summaries (map=TRUE for seconds and probability scale)
summary(single_ss, map = TRUE)

# Recovery of true parameters from simulated data
p_vector
credint(single_ss)
recovery(single_ss, p_vector)

# ------------------------------------------------------------------------------
# Posterior Predictive Checks
# ------------------------------------------------------------------------------

# Simulate 100 datasets from posterior predictive distribution
ppsamps <- predict(single_ss, n_post = 100)
head(ppsamps)
plot_cdf(ppsamps, factors = "S")

# ------------------------------------------------------------------------------
# Compute Mean SSRT
# ------------------------------------------------------------------------------

# In practice, one typically wants to compute mean SSRT.
# In BEESTS, mean SSRT can be computed by adding the samples of mu and tau parameters (on the seconds scale)
# of the stop runner's finishing time distribution, since mu+tau is the mean of the ex-Gaussian distribution.
# SSRT = muS + tauS (mean of ex-Gaussian stop distribution)
pars <- get_pars(single_ss, map = TRUE, merge_chains = TRUE, return_mcmc = FALSE)
hist(pars["muS", , ] + pars["tauS", , ], xlab = "seconds", main = "Mean SSRT")
abline(v = mean(pars["muS", , ] + pars["tauS", , ]), lty = 3, lwd = 5, col = "darkblue")
mean_ssrt <- mean(pars["muS", , ] + pars["tauS", , ])

# Compare this model estimate to the integration method
# Integration method for SSRT estimation
presp <- mean(!is.na(stop_dat$rt))
ssd <- mean(stop_dat$SSD)
signal.resp.rt <- mean(stop_dat$rt, na.rm = TRUE)
# Handle missing RTs by assigning max observed RT
nosignal_resp <- go_dat[!is.na(go_dat$rt), ]
go_dat$rt.true <- ifelse(is.na(go_dat$rt), max(nosignal_resp$rt), go_dat$rt)
# Determine nth quantile of Go RT distribution
nthRT_all <- quantile(go_dat$rt.true, probs = presp, type = 6)
SSRT_integration <- nthRT_all - ssd
SSRT_integration
abline(v = as.numeric(SSRT_integration), lty = 3, lwd = 5, col = "darkgreen")
trueSSRT <- exp(p_vector[5]) + exp(p_vector[7]) 
abline(v = trueSSRT, lty = 3, lwd = 5, col = "darkred")

legend("topright",
       legend = c("True SSRT", "Posterior Mean", "Integration Method"),
       col = c("darkred", "darkblue", "darkgreen"),
       lwd = 2,
       bty = "n",
       lty = 3) 

# ------------------------------------------------------------------------------
# Compare Observed vs Model-Implied Go RTs
# ------------------------------------------------------------------------------

# Go RT is not directly given by model parameters.
# It results from a race between the competing go runners (here, matching vs. mismatching go runners).
# To derive mean Go RTs from the model estimates, one needs to simulate the race process
# or use posterior predictive samples.

# To see this, we first plot model-implied latencies from parameters
# Plot mu + tau for matching and mismatching go runners
# They represent the means of the go runner finishing time distributions, respectively

# Model-implied finishing times from parameter estimates and observed Go RTs
par(mfrow = c(2, 2))
hist(pars["mu_lMTRUE", , ] + pars["tau", , ], xlab = "Seconds",
     main = "Mean matching go runner vs. observed correct go RT", xlim = c(0.65, 0.9))
abline(v = mean(pars["mu_lMTRUE", , ] + pars["tau", , ]), lty = 3, lwd = 5, col = "darkblue")
abline(v = mean(go_dat$rt[go_dat$R == go_dat$S], na.rm = TRUE), lty = 3, lwd = 5, col = "darkgreen")
legend("topleft",
       legend = c("Observed", "Mean Go runner"),
       col = c("darkgreen", "darkblue"),
       lwd = 2,
       bty = "n",
       lty = 3) 


hist(pars["mu_lMFALSE", , ] + pars["tau", , ], xlab = "Seconds",
     main = "Mean mismatching go runner vs. observed incorrect go RT", xlim = c(0.7, 1.05))
abline(v = mean(go_dat$rt[go_dat$R != go_dat$S], na.rm = TRUE), lty = 3, lwd = 5, col = "darkgreen")
abline(v = mean(pars["mu_lMFALSE", , ] + pars["tau", , ]), lty = 3, lwd = 5, col = "darkblue")
legend("topleft",
       legend = c("Observed", "Mean Go runner"),
       col = c("darkgreen", "darkblue"),
       lwd = 2,
       bty = "n",
       lty = 3) 

# Clearly, these parameter-based estimates do not represent the observed RTs
# because the RTs that are observed represent the result of many race processes, 
# based on the joint posterior distribution of the model parameters. They are the winning
# RT subsets only, whereas the finishing time distributions show the entire distribution
# including the losers.

# Posterior predictive samples (ppsamps) are simulated datasets that essentially
# represent the result of many race processes.
# Let's compute mean RTs for correct Go responses for each posterior predictive sample, 
# resulting in the posterior predictive distribution of mean go RT.

# Posterior predictive Go RTs (from race outcomes)
correct_go <- ppsamps[!is.finite(ppsamps$SSD) & ppsamps$S == ppsamps$R, ]
mean_correct_go <- aggregate(rt ~ postn, data = correct_go, FUN = mean, na.rm = TRUE)
hist(mean_correct_go$rt, xlim = c(0.65, 0.9), xlab = "seconds", main = "Correct Go-RT")
abline(v = mean(go_dat$rt[go_dat$R == go_dat$S], na.rm = TRUE), lty = 3, lwd = 5, col = "darkgreen")
abline(v = mean(correct_go$rt, na.rm = TRUE), lty = 3, lwd = 5, col = "darkorange")
legend("topright",
       legend = c("Predicted", "Observed"),
       col = c("darkorange", "darkgreen"),
       lwd = 2,
       bty = "n") 

incorrect_go <- ppsamps[!is.finite(ppsamps$SSD) & ppsamps$S != ppsamps$R, ]
mean_incorrect_go <- aggregate(rt ~ postn, data = incorrect_go, FUN = mean, na.rm = TRUE)
hist(mean_incorrect_go$rt, xlim = c(0.7, 1.05), main = "Incorrect Go-RT")
abline(v = mean(go_dat$rt[go_dat$R != go_dat$S], na.rm = TRUE), lty = 3, lwd = 5, col = "darkgreen")
abline(v = mean(incorrect_go$rt, na.rm = TRUE), lty = 3, lwd = 5, col = "darkorange")
legend("topleft",
       legend = c("Predicted", "Observed"),
       col = c("darkorange", "darkgreen"),
       lwd = 2,
       bty = "n") 
# As you can see, the green line (posterior predictive mean) closely approximates 
# the observed correct and incorrect mean go RTs.

# ------------------------------------------------------------------------------
# Hierarchical Fitting: Real Dataset
# ------------------------------------------------------------------------------

# We'll now fit data from the following paper:

# Matzke, D., Curley, S., Gong, C. Q., & Heathcote, A. (2019). Inhibiting responses 
# to difficult choices. Journal of Experimental Psychology: General, 148(1), 
# 124–142. https://doi.org/10.1037/xge0000525

# In this study, participants completed a random dot motion task with a difficulty manipulation.
# The stop-signal, a gray square boarder around the go stimulus, appeared randomly on around 29% of the trials.
# Here we'll re-analyze the data and conduct a model comparison.


# data from Matzke et al. (2019); datgf
names(datgf)[c(1, 6)] <- c("subjects", "rt")
datgf <- datgf[, c("subjects", "S", "D", "R", "SSD", "rt")]
datgf$R[datgf$R == "NR"] <- NA
datgf$R <- factor(tolower(as.character(datgf$R)))
head(datgf)

# Sample size, trial numbers, and observed RT distributions for each stimulus and difficulty condition
length(unique(datgf$subjects))
table(datgf$subjects, datgf$D)
plot_density(datgf, factors = c("D", "S"))

# Define simple model with effects in mu.go
# Note that this is not (yet) the model reported in the paper.
ADmat <- cbind(d = c(0.5, -0.5))
designMu <- design(data = datgf, model = SSexG,
                   matchfun = function(d) d$S == d$lR,
                   contrasts = list(D = ADmat, lM = ADmat),
                   formula = list(mu ~ D * lM, sigma ~ 1, tau ~ 1,
                                  muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))

# mu: the grand mean/intercept parameter
# mu_Dd: the effect of difficulty on mu
# mu_lMd: the difference between the matching and mismatching go runners in mu
# mu_Dd:lMd: the interaction effect of math/mismatch and difficulty on mu

# Prior similar to Matzke et al. (2019)
prior_ss <- prior(designMu, mu_mean = c(log(.7), 0, 0, 0, log(0.5), log(0.5), log(0.7), log(0.5), log(0.5), qnorm(0.05), qnorm(0.05)),
                  mu_sd = c(rep(1, 9), 0.5, 0.5), type = "standard")

plot(prior_ss)
emc <- make_emc(datgf, designMu, prior = prior_ss, compress = FALSE, rt_resolution = 1e-20)
#hierarchical_mu <- fit(emc, cores_per_chain = 15, fileName = "hierarchical.RData")

# Check convergence and posterior summaries
check(hierarchical_mu)

# traceplots of the population-level parameters
plot(hierarchical_mu,selection = "mu", layout=c(2,3))
plot(hierarchical_mu,selection = "sigma2", layout=c(2,3))
plot(hierarchical_mu,selection = "correlation", layout=c(2,3))

# Posterior and prior distribution of the group-level mean parameters
plot_pars(hierarchical_mu, selection = "mu",layout=c(2,3))
# Posterior and prior distribution of the group-level variance parameters
plot_pars(hierarchical_mu, selection = "sigma2",layout=c(2,3))

# pp_samps_mu <- predict(hierarchical_mu, n_post = 100, n_cores = 15)
# Posterior predictive checks: there is considerable misfit!
plot_cdf(datgf, pp_samps_mu, factors = c("S", "D"))

# Inhibition functions with posterior predictives (slow!):
plot_ss_if(datgf, post_predict = pp_samps_mu, probs = seq(0, 1, length.out = 5))
##The misfit is worse in the easy condition
plot_ss_if(datgf, post_predict = pp_samps_mu, probs = seq(0, 1, length.out = 5),
           factors = "D")


# SRRT as a function of SSD bins with posterior predictives:
plot_ss_srrt(datgf, post_predict = pp_samps_mu,
             probs = seq(0, 1, length.out = 5))
plot_ss_srrt(datgf, post_predict = pp_samps_mu,
             probs = seq(0, 1, length.out = 5),
             factors = "D")


# ------------------------------------------------------------------------------
# Full Model (mu, sigma, tau modelled)
# ------------------------------------------------------------------------------

# We will try to reduce the misfit by also allowing sigma and tau to vary with match/mismatch
# This is also the model reported in the paper!
designFull <- design(data = datgf, model = SSexG,
                     matchfun = function(d) d$S == d$lR,
                     contrasts = list(D = ADmat, lM = ADmat),
                     formula = list(mu ~ D * lM, sigma ~ D * lM, tau ~ D * lM,
                                    muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1))

prior_ss <- prior(designFull, mu_mean = c(log(.7),
                                          0,
                                          0,
                                          0,
                                          log(0.5),
                                          0,
                                          0,
                                          0,
                                          log(0.5),
                                          0,
                                          0,
                                          0,
                                          log(0.7),
                                          log(0.5),
                                          log(0.5),
                                          qnorm(0.05),
                                          qnorm(0.05)),
                  mu_sd = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, .5, .5),
                  type ="standard")

emc <- make_emc(datgf, designFull, prior = prior_ss, compress = FALSE, rt_resolution = 1e-20)
#hierarchical_full <- fit(emc, cores_per_chain = 15, fileName = "hierarchical-full.RData")

check(hierarchical_full)

# Posterior predictive checks
# pp_samps_full <- predict(hierarchical_full,  n_post = 100, n_cores = 15)
plot_cdf(datgf, pp_samps_full, factors = c("S", "D"))
plot_ss_if(datgf, post_predict = pp_samps_full, factors = "D")
plot_ss_srrt(datgf, post_predict = pp_samps_full, factors = "D")

# Model comparison also suggests that the more complex model is preferred
compare(list(mu_model = hierarchical_mu, full_model = hierarchical_full), BayesFactor = FALSE)
