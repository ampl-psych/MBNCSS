#### GNG DDM ----
library(EMC2)

# Gomez, P., Ratcliff, R., & Perea, M. (2007). A Model of the Go/No-Go Task.
# Journal of Experimental Psychology: General, 136(3), 389–413.
# https://doi.org/10.1037/0096-3445.136.3.389

DDMGNG()

Rnogo <- function(d)factor(rep("nogo",nrow(d)),levels=c("go","nogo"))
Rgo <- function(d)factor(rep("go",nrow(d)),levels=c("go","nogo"))

# Timout
maxt <- 100  # no timeouts, 100 second response window
TIMEOUT <- function(d)rep(maxt,nrow(d))

dGNG <- design(factors=list(subjects=1:1,S=c("go","nogo")),
               Rlevels=c("go","nogo"),
    functions=list(
      TIMEOUT=TIMEOUT,
      Rnogo=Rnogo, # no go response level
      Rgo=Rgo),  # go response level
    contrasts=list(S=cbind(d=c(-1,1))),
    formula=list(v~S,a~1, Z~1, t0~1),
    model=DDMGNG)

# Lets look at some simulated data
p <- sampled_pars(dGNG,doMap=FALSE)
# Very high accuracy
p[] <- c(0,10,log(2),qnorm(.5),log(.4))  # ~100% accuracy

# Always correct because of high accuracy, correct for nogo
# => NA in RT
dat <- make_data(p,dGNG,10)
dat

# dadm and the functions passed to design
dadm <- EMC2:::design_model(dat,dGNG,compress=FALSE)
dadm

# Same number of rows as in data
nrow(dat)
nrow(dadm)

# Look at the functions
Rnogo(dadm)
Rgo(dadm)
TIMEOUT(dadm)

# Now lets use a stupidly short timeout
maxt <- .5
TIMEOUT <- function(d)rep(maxt,nrow(d))
dGNG <- design(model=DDMGNG,factors=list(
    subjects=1:1,S=c("go","nogo")),Rlevels=c("go","nogo"),
    functions=list(TIMEOUT=TIMEOUT,Rnogo=Rnogo,Rgo=Rgo),
    contrasts=list(S=cbind(d=c(-1,1))),formula=list(v~S,a~1, Z~1, t0~1))

# Many go trials have no response.
make_data(p,dGNG,10)


# No re-run design with a more realistic timeout
maxt <- 2
TIMEOUT <- function(d)rep(maxt,nrow(d))
dGNG <- design(model=DDMGNG,factors=list(
    subjects=1:1,S=c("go","nogo")),Rlevels=c("go","nogo"),
    functions=list(TIMEOUT=TIMEOUT,Rnogo=Rnogo,Rgo=Rgo),
    contrasts=list(S=cbind(d=c(-1,1))),formula=list(v~S,a~1, Z~1, t0~1))

# and more realistic error-prone responding, and also a bias to nogo
p[] <- c(0.5,1,log(2),qnorm(.4),log(.3))

# Make lots of data
dat <- make_data(p,dGNG,200)

# Now we have errors, more so for go
tapply(dat$S==dat$R,dat$S,mean)

# and fewer timeouts, but still quite a few
mean(is.na(dat$rt[dat$S=="go"]))

# Truncation clear at 2sec
plot_density(dat,factor="S")

# A little more obvious in a histogram
hist(dat$rt)


# Lets fit some data
emc <- make_emc(dat,dGNG,type="single")
# Processing data set 1
# Likelihood speedup factor: 125.8 (159 unique trials)


# NB: 8x or so slower than standard DDM
fit(emc,fileName="gng3.RData")


load("gng1.RData")
load("gng2.RData")
load("gng3.RData")
recovery(emc,p)

# Fit is excellent to RT distributions
pemc <- predict(emc,n_cores=9)
plot_cdf(emc,pemc,factors="S")
plot_density(emc,pemc,factors="S")

# and accuracy
tapply(dat$S==dat$R,dat$S,mean)
tapply(pemc$S==pemc$R,pemc$S,mean)

# and timeouts
mean(is.na(dat$rt[dat$S=="go"]))
mean(is.na(pemc$rt[pemc$S=="go"]))








