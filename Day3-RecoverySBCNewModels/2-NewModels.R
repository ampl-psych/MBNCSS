rm(list=ls())
# remotes::install_github("ampl-psych/EMC2",ref="dev")
library(EMC2)

# This lesson describes how to go beyond the standard models offered by EMC2, 
# It describes how to implement and test new models, model parameteriations and 
# ways of bounding parameters, the latter being sometimes necessary when there 
# are convergence problems. It also gives insight into the inner working of EMC2.

# We will use the Basic lesson data to illustrate (already filtered)
load("FlankerData.RData")

# We will use this contrast later
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))

### Profile plots and augmented data  ----

# When implementing a new model it is easy to make coding mistakes. A useful 
# check is to see if random and likelihood functions (both of which use the 
# transforms) are still consistent by simulating a large amount of data data and 
# making likelihood "profile" plots. These plots hold all parameters at the data 
# generating values except one, which is varied to check if the likelihood peaks 
# above the data generating value.

# To illustrate we use the RDM to simulate large data set so that 
# sampling noise is at a minimum. The result is 40,000 trials, 5000 in each of 
# the 8 design cells.
designRDMt0 <- design(data=dat,model=RDM,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,t0~lR+E))
sdat <- make_data(p_vectort0,designRDMt0,n_trials=5000)

# We see the profiles peak above the data-generating values. By default the 
# calculates the likelihood at 100 points and prints out a table containing the
# true values, the value for which the likelihood was maximum (a crude estimate
# of the best fitting parameter) and the difference between the two. We see
# this "miss" is very small in all cases, consistent with the model being  
# correctly specified. 
profile_plot(sdat, designRDMt0, p_vector = p_vectort0, layout = c(2,5),n_cores=12)

# Profile plot is initially quite slow because it must first create a "dadm"
# object using the EMC2's internal "design_model" function, which requires 
# several steps
# 1) data augmentation
# 2) data compression
# 3) binding the augmented and compressed data with a design and model

# For a race model data augmentation consists of stacking as many copies of the
# data together as there are accumulators in the model. For the binary RDM
# example here this doubles the size. To illustrate we first make a 
# dadm without compression and with only minimal rounding of rt's: 
dadm1 <- EMC2:::design_model(sdat,designRDMt0,compress=FALSE, rt_resolution=.001)
nrow(sdat)
nrow(dadm1)

# We see that each row of sdat is now repeated twice and three columns added:
# a) lR factor = latent response (with the same levels as R)
# b) lM logical = latent match, indicating which accumulator corresponds to the
#       the correct response given the stimulus (created by matchfun) and 
# c) winner logical = which accumulator corresponds to the response
head(dadm1)

# Next the data are compressed as determined by the rt_resolution parameter. 
# The default is 0.02 seconds, but here we use 0.05, so that only 1219 trials
# are kept that have the potential for having different likelihoods given any
# possible parameter vector. Calculating the compression can be time 
# consuming.
dadm2 <- EMC2:::design_model(sdat,designRDMt0,rt_resolution=.02)

# New now see that RTs are rounded to the nearest 0.05 and the dadm is much
# smaller (2 x 1219 rows)
nrow(dadm2)
head(dadm2)

# We see that if we now pass the dadm to profile plot it runs much more quickly
# as it does not have to do the time consuming compression operation. Although
# you mist still pass sdat this is no longer used.
profile_plot(sdat, designRDMt0, p_vector = p_vectort0, layout = c(2,5),
             dadm=dadm2,n_cores=12)

# Passing a dadm is the only way to examine the effect of compression on the
# profile. Generally we have found that compression of .05 or less has little
# impact on parameter recovery for race models.

#### Changing scales and interactions  ----

# Suppose that you have reason to believe that non-decision time (t0) is slower
# for right than left responses (e.g., due to the right response button being
# stiffer than the left response button) and that non-decision time is faster
# under speed than accuracy emphasis, and crucially for this illustration, that
# these two effects are additive on the natural (i.e., seconds) scale. For 
# example, suppose that the two effects have equal and opposite effects of 0.1s. 
# Additively implies the two effects exactly cancel for right (0.1s slower) 
# speed emphasis (0.1s faster) responses. 

# Unfortunately there is no easy way to make the effects additive using the 
# standard RDM model as it estimates t0 on a  log scale, and an additive model 
# (t0 ~ lR + E) on the log scale is not additive on the natural scale.
#
# To illustrate consider the standard model with an additive t0 effect.

designRDMt0 <- design(data=dat,model=RDM,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,t0~lR+E))

# Make a parameter vector where the baseline t0 = 0.2, right is 0.1s slower
# (so on a log scale we add log(1.5) which is the same as multiply by 1.5 on 
# the natural scale, an increase of 0.1 (.2*1.5=0.3) and speed is 0.1s faster, 
# (so on a log scale we add log(0.5) which is the same as multiply by 0.5 on 
# the natural scale, a decrease of 0.1 (.2*0.5=0.1).

p_vectort0 <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
                v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
                t0=log(.2),t0_lRright=log(1.5),t0_Espeed=log(0.5))

# Examining t0 for right speed we that rather than being 0.2 (i.e., the two
# effects cancel) it is 0.15 ( exp(log(.2) + log(1.5) + log(.5)) = .2*1.5*.5)
mapped_pars(designRDMt0,p_vectort0)

# We can address this problem by modifying the transforms argument
# of design.
designRDMt0N <- design(data=dat,model=RDM,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,t0~lR+E),
  transform=list(func=c(t0="identity")))

p_vectort0N <- c(B=log(2),B_Espeed=log(1),B_lRright=log(1),v=log(2),
                v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
                t0=.2,t0_lRright=.1,t0_Espeed=-.1)

# We get an additive effect on the real scale, .3/.2 for accuracy and .2/.1 for
# speed.
mapped_pars(designRDMt0N,p_vectort0N)

dat <- make_data(p_vectort0N,designRDMt0N,n_trials=1000)

# Can see steps in t0 likelihood as it falls below the .05 bound
profile_plot(dat,designRDMt0N,p_vectort0N,layout=c(2,5),n_cores=12)


# Estimating t0 on the natural scale raises the possibility of sampling negative
# times. EMC2 models usually have default bounds that prevent non-nonsensical 
# values from occurring. The bounds are given on the natural scale and have two 
# types, a matrix with rows giving ranges within which values are allowable 
# ($minmax) and a set of specific exception points ($exception) that are allowed
# outside this range.

# We see t0 is already protected from values < .05 seconds
RDM()$bound

# NB1:  When the data are not very constraining highly variable sampling can
#       occur near transformation bounds, making it hard to achieve convergence. 
#       In such cases it can be useful to tweak the lower bound (e.g., you might 
#       use pars[,"A"] > 1e-3) to archive convergence.


# Suppose we thought a higher bound of 0.15 seconds is more reasonable. That
# will cause p_vectort0N to break the bounds for t0 in the SPEED conditions 

designRDMt0Nb <- design(data=dat,model=RDM,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,t0~E+lR),
  transform=list(func=c(t0="identity")),
  bound=list(minmax = cbind(t0=c(.15,Inf))))


# Make a very small amount of data and a corresponding dadm
dat <- make_data(p_vectort0N,designRDMt0N,n_trials=1)
dadm <- EMC2:::design_model(dat,designRDMt0Nb)

# Here is how EMC2 calculates likelihoods in R
EMC2:::calc_ll_R(p_vectort0N, attr(dadm, "model")(), dadm)

# Lets take a look at the function that computes likelihoods and tweak it so 
# we get ll for each trial returned.
my_log_likelihood_race <- function (pars, dadm, model, min_ll = log(1e-10)) 
{
    if (any(names(dadm) == "RACE")) {
        pars[as.numeric(dadm$lR) > as.numeric(as.character(dadm$RACE)), 
            ] <- NA
    }
    if (is.null(attr(pars, "ok"))) {
        ok <- !logical(dim(pars)[1])
    }
    else ok <- attr(pars, "ok")
    lds <- numeric(dim(dadm)[1])
    lds[dadm$winner] <- log(model$dfun(rt = dadm$rt[dadm$winner], 
        pars = pars[dadm$winner, ]))
    n_acc <- length(levels(dadm$R))
    if (n_acc > 1) 
        lds[!dadm$winner] <- log(1 - model$pfun(rt = dadm$rt[!dadm$winner], 
            pars = pars[!dadm$winner, ]))
    lds[is.na(lds) | !ok] <- min_ll
    if (n_acc > 1) {
        ll <- lds[dadm$winner]
        if (n_acc == 2) {
            ll <- ll + lds[!dadm$winner]
        }
        else {
            ll <- ll + apply(matrix(lds[!dadm$winner], nrow = n_acc - 
                1), 2, sum)
        }
        ll[is.na(ll)] <- min_ll
        return(pmax(min_ll, ll[attr(dadm, "expand")]))
    }
    else return(pmax(min_ll, lds[attr(dadm, "expand")]))
}

# To use the underlying log_likelihood we need to expand the parameter vector
# to match the dadm
pars <- EMC2:::get_pars_matrix(p_vectort0N, dadm, attr(dadm, "model")())
View(cbind(dadm,pars))

# We see tat the SPEED trials are flagged as having parameters that break bounds.
pars

# We see min_ll = log(1e-10) in speed cells
cbind(dat,ll=my_log_likelihood_race(pars,dadm,attr(dadm, "model")()))

# Without bounds the speed cells have standard likelihoods
dadm <- EMC2:::design_model(dat,designRDMt0N)
pars <- EMC2:::get_pars_matrix(p_vectort0N, dadm, attr(dadm, "model")())
cbind(dat,my_log_likelihood_race(pars,dadm,attr(dadm, "model")()))


### A new EMC2 model ----

# Now lets return to making a new model. Models in EMC2 are functions returning 
# a list of values and functions. Lets look at the LBA model and go through its 
# components that we have to change to make an LBA with a different 
# parameterization.
LBA()

# Standard models like the RDM are also implemented in C++ through Rcpp (see the
# model_RDM.h in the package src directory), but this is only used when 
# sampling, where speed is critical. Otherwise (when simulating data etc.). 

# The new model will make will not have a C++ version. The presence of a C++ 
# version is signaled by the presence of the following list entry, so we will
# simply omit this entry from our new model
LBA()$c_name

# Note that even without this EMC2 still uses the C++ versions of functions
# to generate random samples from (rfun in the model list), and the density 
# (dfun) and cumulative density (pfun) of a single LBA accumulator defined at 
# the bottom of the R model function, it is only the likelihood function call 
# that is not in C.


# Our aim here is to parameter A as a fraction of B rather. To do that a
# couple of things need to change. Lets first make a copy of the model that
# we will alter.

# Here we will change the Ttransform function of the standard LBA model so in 
# the new model B is an estimate of the absolute  threshold (b) and A falls a 
# portion of the way between 0 and b

# First note the added "EMC2:::" to make internal EMC2 functions available.

LBApA <- function () { list(
    type="RACE",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "pnorm",t0 = "exp")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(0,1),B=c(0,Inf),
                            t0=c(0.05,Inf)),exception=c(A=0)),
    Ttransform = function(pars,dadm) {
      pars[, "A"] <- pars[,"A"]*pars[,"B"]
      pars <- cbind(pars, b = pars[,"B"])
      pars
    },
    # Random function for racing accumulator
    rfun=function(data,pars) EMC2:::rLBA(data$lR,pars,posdrift=TRUE,ok = attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) EMC2:::dLBA(rt,pars,posdrift = TRUE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) EMC2:::pLBA(rt,pars,posdrift = TRUE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      EMC2:::log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

# NB1: The Ttransform function is also used to combine sampled parameters once 
#     they are  mapped back to the natural scale in order to create parameters 
#     required by the likelihood. For example, in the LBA it has the line
#              pars <- cbind(pars,b=pars[,"B"] + pars[,"A"]),
#     which makes the b, the response threshold, by adding A (the range of the 
#     start-point distribution, to B, the gap from the top of the start-point 
#     distribution to the threshold (in the RDM this is done inside rfun,
#     dfun and pfun).

# NB2: Although we will not change the other elements here are some brief notes
#      1) type: RACE for race models, DDM for the DDM. 
#      2) p_types: vector of default values with the names of the model's 
#         parameter types.
#      3) log_likelihood: the likelihood function to use with the model (for
#         RACE type models usually EMC2:::log_likelihood_race, for DDM type
#         models usually EMC2:::log_likelihood_ddm)
#      4) The dfun and pfun elements just directly call the underlying RDM 
#         density and cumulative density functions, but rfun also calculates the
#         same "ok" test as Ttransform. This is used with data simulation (where
#         inadmissible values can also cause problems). lR is the latent 
#         response factor, which is used to pass the names of responses to the
#         data simulation functions. 



# Next we want to keep A on 0-1 by putting A on the probit scale 
# (transform$func["A"] <- "pnorm") 
LBA()$transform
LBApA()$transform

# and $bound$minmax[,"A"] to the unit interval
LBA()$bound
LBApA()$bound

# Finally, note that the new model has no c_name to indicate it does not 
# have a C version 
LBA()$c_name
LBApA()$c_name



# Now lets make a design 
designLBApA <- design(data=dat,model=LBApA,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,A~1,t0~1,sv~1),constants=c(sv=log(1)))

p_vector <- c(B=log(2),B_Espeed=log(.9),B_lRright=log(1),v=log(3),
  v_lMd=log(2),v_CIincongruent=log(1),'v_lMd:CIincongruent'=log(1),
  t0=log(.2),A=qnorm(.25))

# We see that A is a quarter of B, which makes it smaller for SPEED than 
# ACCURACY as caution is less for SPEED. This proportionality assumption
# can be useful / psychologically grounded so having the freedom to be able
# to alter parameterizations flexibly is important.
mapped_pars(designLBApA,p_vector)


# Lets now check the profiles
sdat <- make_data(p_vector,designLBApA,n_trials=5000)
plot_density(sdat,factors="S")

profile_plot(sdat, designLBApA, p_vector = p_vector, 
             layout = c(2,5),n_cores=12)
# We see that profiles and table are consistent with out implementation being
# correct. 

# When you make a new model it is a good idea to do parameter recovery studies 
# to make sure it is identifiable in your design of interest. Do this for the
# LBApA model

emc <-  make_emc(sdat,designLBApA,type="single")
save(emc,file="samples/newmodel/LBApA.RData")
# run by fit_LBApA.R

# Recovery is ok.
load("samples/newmodel/LBApA.RData")
recovery(emc,true_pars=matrix(p_vector,nrow=1))

# Speed comparison with a standard LBA
designLBA <- design(data=dat,model=LBA,
  matchfun=function(d)d$S==d$lR,contrasts=list(lM=ADmat),
  formula=list(B~E+lR,v~lM*CI,A~1,t0~1,sv~1),constants=c(sv=log(1)))
emc <-  make_emc(sdat,designLBA,type="single")
save(emc,file="samples/newmodel/LBA.RData")
# run by fit_LBA.R

# Looking at the fitLBA.Rout and fitLBApA.Rout files we see that the new model 
# is much slower because it lacks a C++ instantiation (18x).

# Altering the C would require changes in src/model_LBA.h in the following two
# functions

# NumericVector dlba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
#   //v = 0, sv = 1, B = 2, A = 3, t0 = 4
#   int n = sum(idx);
#   NumericVector out(n);
#   int k = 0;
#   for(int i = 0; i < rts.length(); i++){
#     if(idx[i] == TRUE){
#       if(NumericVector::is_na(pars(i,0))){
#         out[k] = 0;
#       } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
#         out[k] = dlba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
#       } else{
#         out[k] = min_ll;
#       }
#       k++;
#     }
#   }
#   return(out);
# }
# 
# NumericVector plba_c(NumericVector rts, NumericMatrix pars, LogicalVector idx, double min_ll, LogicalVector is_ok){
#   //v = 0, sv = 1, B = 2, A = 3, t0 = 4
#   int n = sum(idx);
#   NumericVector out(n);
#   int k = 0;
#   for(int i = 0; i < rts.length(); i++){
#     if(idx[i] == TRUE){
#       if(NumericVector::is_na(pars(i,0))){
#         out[k] = 0;
#       } else if((rts[i] - pars(i,4) > 0) && (is_ok[i] == TRUE)){
#         out[k] = plba_norm(rts[i] - pars(i,4), pars(i,3), pars(i,2) + pars(i,3), pars(i,0), pars(i,1), true);
#       } else{
#         out[k] = min_ll;
#       }
#       k++;
#     }
#   }
#   return(out);
# }

#### GNG DDM (if time) ----

# Gomez, P., Ratcliff, R., & Perea, M. (2007). A Model of the Go/No-Go Task. 
# Journal of Experimental Psychology: General, 136(3), 389–413. 
# https://doi.org/10.1037/0096-3445.136.3.389

DDMGNG()


# With TIMEOUT and the response that is nogo passed through to the likelihood

   # Ttransform = function(pars,dadm) {
   #    pars[,"SZ"] <- 2*pars[,"SZ"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min)
   #    pars <- cbind(pars,z=pars[,"Z"]*pars[,"a"], sz = pars[,"SZ"]*pars[,"a"],
   #                  TIMEOUT=dadm$TIMEOUT,Rnogo=as.numeric(dadm$Rnogo))
   #    pars
   #  }

# The standard random function sets rt to NA after timeout and R for "nogo"
# to NA
    # Random function
    rfun=function(data,pars) {
      out <- rDDM(data$R,pars, attr(pars, "ok"))
      out$rt[out$rt>pars[,"TIMEOUT"]] <- NA
      out$rt[as.numeric(out$R)==pars[,"Rnogo"]] <- NA
      out
    }

# And of course a different likelihood
        log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddmgng(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }

EMC2:::log_likelihood_ddmgng

# When we have an rt it is the standard DDM
        isna <- is.na(dadm$rt)
        ok <- attr(pars, "ok") & !isna
        like[ok] <- model$dfun(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
        
# When rt is NA you get 1 - probability of response before timeout 
# (so either a slow response that would have been made on the go boundary or 
# termination on the nogo boundary).

        ok <- attr(pars, "ok") & isna
        like[ok] <- pmax(0, pmin(1, (1 - model$pfun(dadm$TIMEOUT[ok], 
            dadm$Rgo[ok], pars[ok, , drop = FALSE]))))

# Here is the underlying DDM cdf
EMC2:::pDDM

# Which calls the WeinR package
library(WeinR)
WienR:::pWDM
        
# Running the model 
  
dGNG <- design(factors=list(subjects=1:1,S=1:2),Rlevels=1:2,
    functions=list(
      TIMEOUT=function(d)rep(2,nrow(d)),
      Rnogo=function(d)factor(rep(1,nrow(d)),levels=c(1:2)), # no go response level
      Rgo=function(d)factor(rep(2,nrow(d)),levels=c(1:2))),  # go response level
    contrasts=list(S=cbind(d=c(-1,1))),
    formula=list(v~S,a~1, Z~1, t0~1),
    model=DDMGNG)

# Lets look at some simulated data
p <- sampled_pars(dGNG)
p[] <- c(0,1,log(1),qnorm(.5),log(.4))  # ~70% accuracy

# We see non-responses and reponses
make_data(p,dGNG,10)

# If we slow it down (a=4) we get more timeouts 
p[] <- c(0,1,log(4),qnorm(.5),log(.4))  # ~70% accuracy
make_data(p,dGNG,10)

# Now lets compare fitting speed
gng <- make_data(p,dGNG,100)
emc <- make_emc(gng,dGNG,type="single")
system.time(tmp <- init_chains(emc))


dDDM <- design(factors=list(subjects=1:1,S=1:2),Rlevels=1:2,
    contrasts=list(S=cbind(d=c(-1,1))),
    formula=list(v~S,a~1, Z~1, t0~1),
    model=DDM)

ddmd <- make_data(p,dDDM,100)
emc1 <- make_emc(ddmd,dDDM,type="single")

system.time(tmp <- init_chains(emc1)) # ~8x

# As an exercise compare fitting time.
emc <- fit(emc) 
# Time difference of 58.08074 secs

emc <- fit(emc1)  # ~8x
# Time difference of 6.796554 secs

