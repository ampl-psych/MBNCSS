# remotes::install_github("ampl-psych/EMC2",ref="SScpp1")
# remotes::install_github("ampl-psych/EMC2",ref="dev")
library(EMC2)

#### Setup ----

# Use the final selected model from the hierarchical lesson to use in as our
# recovery example.
matchfun=function(data)data$S==data$lR
ADmat <- cbind(d = c(-1/2,1/2))
E2 <- function(data) factor(data$E!="speed",labels=c("speed","nonspeed"))
E_incr <- contr.increasing(3)
colnames(E_incr) <- c("nonspeed", "acc")

# Make a design suitable for a single subject
design_LBABvE2_1 <- design(
  factors=list(subjects=1,E=levels(forstmann$E),S=levels(forstmann$S)),
  Rlevels=levels(forstmann$R),model=LBA,matchfun=matchfun,
  formula=list(v~E*lM,sv~lM,B~E2+lR,A~1,t0~1),
  contrasts=list(lM = ADmat, E = E_incr),
  constants=c(sv=log(1), `v_Eacc:lMd` = 0),
  functions=list(E2=E2))

# Load fits to data
print(load("samples/LBABvE2.RData"))
print(load("samples/LBABvE.RData"))

#### Single participant models ----

### Checking single participant parameter recovery ----

print(load("samples/E2_single.RData"))

# Lets see how well we can estimate parameters from a lot of simulated data.
# This test generally indicates if our model is working as expected, and so
# is used early in model development.

# Lets pick some typical parameter values from our fit to real data
mu <- credint(LBABvE2)$mu[,2]

# 30,000 trials should support accurate recovery
dat <- make_data(mu,design_LBABvE2_1,n_trials=5000)

# Is the conditional likelihood working as we expect? Profiles should peak 
# over true values.
profile_plot(dat,design_LBABvE2_1,mu,layout=c(2,6))

# # Now lets try fitting to the simulated data
# emc <- make_emc(dat,design_LBABvE2_1,type="single")
# E21 <- fit(emc,fileName="samples/E21.RData")

E21 <- get(load("samples/E21.RData"))

# As expected recovery looks excellent
recovery(E21,true_pars=mu)

# With a more realistic design still pretty good profiles 
dat1 <- make_data(mu,design_LBABvE2_1,n_trials=125)
profile_plot(dat1,design_LBABvE2_1,mu,layout=c(2,6))

# # Small sample recovery
# emc <- make_emc(dat1,design_LBABvE2_1,type="single")
# E211 <- fit(emc,fileName="samples/E211.RData")

E211 <- get(load("samples/E211.RData"))

# Ok, but with wider credible intervals
recovery(E211,mu)

### Checking many single-participant recoveries ----

# Now lets try this with lots (200) of single-fit replicates (as subjects)
design_LBABvE2 <- design(
  factors=list(subjects=1:200,E=levels(forstmann$E),S=levels(forstmann$S)),
  Rlevels=levels(forstmann$R),model=LBA,matchfun=matchfun,
  formula=list(v~E*lM,sv~lM,B~E2+lR,A~1,t0~1),
  contrasts=list(lM = ADmat, E = E_incr),
  constants=c(sv=log(1), `v_Eacc:lMd` = 0),
  functions=list(E2=E2))

# We will draw the 200 simulated participants from the full multivariate normal 
# population distribution estimated by our fit.

# First get the covariance matrix
var <- credint(LBABvE2,selection="sigma2")$sigma2[,2]
vcov <- diag(var)
dimnames(vcov) <- list(names(var),dimnames=names(var))
cov <- lapply(credint(LBABvE2,selection="covariance"),\(x)x[,2])
for (i in names(var)) vcov[i,colnames(vcov)!=i] <- cov[[i]] 


# # Sample 200 subjects and make data
# E2par <- make_random_effects(design_LBABvE2,group_means=mu,covariances = vcov)
# 
# # Now use the sampled parameters to make data
# datE2 <- make_data(E2par,design_LBABvE2,n_trials=125)

# # Run by fit_E2E2.R
# emc <- make_emc(datdesign_LBABvE2,type="single")
# save(emc,file="E2E2.RData")

# Load fits
E2E2 <- get(load("samples/E2E2.RData"))

# Generally good, some wild estimates but also with large credible intervals.
recovery(E2E2,E2par,layout=c(2,6))

# Lets look at the "coverage" statistic: on what proportion of replicates was
# the true value inside the CI?

# Default 95% CI coverage
CIs <- credint(E2E2,probs=c(.025,.5,.975))
coverage <- setNames(numeric(ncol(E2par)),names(CIs))
for (i in names(CIs)) {
  p <- E2par[,i]
  coverage[i] <- mean((CIs[[i]][,1] < p) & (CIs[[i]][,3]>p))
}

# Generally close to the nominal 95% value.
100*coverage

# save(mu,dat,emc,dat1,emc1,E2par,datE2,file="samples/E2_single.RData")

### A more general and fine-grained approach: simulation-based calibration (SBC) ----

## https://cran.r-project.org/web//packages/EMC2/vignettes/Simulation-based-Calibration.html


# Check if a model is sampling properly in the space defined by a prior

# We use the posterior samples from the fits to data to get a prior.
E2E2 <- get(load("samples/E2E2.RData"))

# Get broadly representative prior from all subject-level parameter
pars <- parameters(E2E2,selection="alpha")

# Get their mean and SD
pmean <- apply(pars[,-1],2,mean)
psd <- apply(pars[,-1],2,sd)

# Make a prior
priorE2 <- prior(design_LBABvE2,mu_mean=pmean,mu_sd=psd,type="single")
plot(priorE2,layout=c(2,7))

# # Run SBC (see also run_sbc.R, ~ 3.5hrs using 15 cores)
# SBC_E2 <- run_sbc(design_LBABvE2, priorE2, replicates = 200, trials = 125, plot_data = FALSE,
#                   iter = 1000, n_post = 1000, fileName = "SBC_E2.RData",
#                   cores_per_chain = 25)

# Load results
sbcE2 <- get(load("samples/SBC/SBC_E2.RData"))

# Distribution of proportional ranks of true values in posteriors. For well
# calibrated estimation this should be uniform. 
plot_sbc_hist(sbcE2,layout=c(2,6))

# Plot shows 95% credible intervals but hard to judge deviations as histogram 
# bin width affects appearance. 

# A better method is to use cumulative density functions, which should increase
# linearly. To make deviations easy to see the following plot subtracts the
# linear trend, so a flat line is expected, with the blue egg representing the
# 95% credible region.
plot_sbc_ecdf(sbcE2,layout=c(2,6))


#### Hierarchical models ----

### Checking hierarchical parameter recovery ----

# forstmann design
design_LBABvE2 <- design(data = forstmann,model=LBA,matchfun=matchfun,
                        formula=list(v~E*lM,sv~lM,B~E2+lR,A~1,t0~1),
                        contrasts=list(lM = ADmat, E = E_incr),
                        constants=c(sv=log(1), `v_Eacc:lMd` = 0),
                        functions=list(E2=E2))

# Make data from posterior medians of alpha, exact same design as the
# original data

# E2dat <- make_data(LBABvE2)
# save(E2dat,file="samples/hier/E2dat.RData")
load("samples/hier/E2dat.RData")

# emc <- make_emc(E2dat,design_LBABvE2)
# save(emc,file="E2.RData")
load("samples/hier/E2.RData")

# Fits converged
check(emc)

# Population means are well recovered
recovery(emc,true_pars=LBABvE2)  # default selection="mu"

# Clear shrinkage on variance estimates
recovery(emc,true_pars=LBABvE2,selection="sigma2")

# Quite a deal of uncertainty but little apparent bias
recovery(emc,true_pars=LBABvE2,selection="correlation",flatten=TRUE)

# Shrinkage apparent in most parameters, clearest in v
recovery(emc,true_pars=LBABvE2,selection="alpha",layout=c(2,3))


### Checking hierarchical model selection ----

# Lets test how well we can discriminate the "E2" model from the "E" model

# Simpler E model

design_LBABvE <- design(data=forstmann,model=LBA,matchfun=matchfun,
  formula=list(v~lM*E2,sv~lM,B~E2+lR,A~1,t0~1),
  contrasts=list(lM = ADmat),constants=c(sv=log(1)),
  functions=list(E2=E2))


# # Now fit all 2 x 2 combinations 
# 
# # Fit E model to E2 data
# emc <- make_emc(E2dat,design_LBABvE)
# save(emc,file="samples/E2_E.RData")
# 
# 
# # Fit true
# emc <- make_emc(Edat,design_LBABvE)
# save(emc,file="samples/E.RData")
# 
# # Fit E2 model to E data
# emc <- make_emc(Edat,design_LBABv2E)
# save(emc,file="samples/E_E2.RData")


E2 <- get(load("samples/hier/E2.RData"))
E2E <- get(load("samples/hier/E2_E.RData"))
E <- get(load("samples/hier/E.RData"))
EE2 <- get(load("samples/hier/E_E2.RData"))

# The more complex model is strongly selected over the simpler one (slowish)
compare(list(E2E2=E2,E2E=E2E),cores_for_props = 4)
#          MD wMD    DIC wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# E2E2 -16620   1 -17166    1 -16932     1        233 -17399 -17578 -17632
# E2E  -16453   0 -16915    0 -16686     0        229 -17144 -17306 -17373

# At the individual level results are more equivocal (12/19)
compare_subject(list(E2E2=E2,E2E=E2E))

# Reflects both greater uncertainty and that some participants have small 
# values of v_Eacc (e.g., ku4t   rt3t) 
sort(credint(LBABvE2,selection="alpha")$v_Eacc[,2])

# However, the simpler model is not perfered to the more complex one when it
# is true (slow). 
compare(list(EE=E,EE2=EE2),cores_for_props = 4)
#         MD  wMD    DIC  wDIC   BPIC wBPIC EffectiveN  meanD  Dmean   minD
# EE  -16475 0.29 -16943 0.213 -16721 0.213        222 -17165 -17321 -17387
# EE2 -16477 0.71 -16946 0.787 -16724 0.787        222 -17168 -17322 -17389

# Individual level selection also biased to the more complex model (e.g., DIC 7/19), 
# but much more clearly equivocal in model weights.
compare_subject(list(EE=E,EE2=EE2))


# These results suggest some caution about the support for the E2 model in 
# the empirical data.


#### Exercise ----

# 1. Repeat the SBC but allow for correlations in the prior with values taken
#    from the empirical fits  (see samples/run_sbcr.R) for a solution

