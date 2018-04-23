###### Merging SDM and NEVA: 03-30-2018
################
## Key components: No longer sampling using MCMC update -- exhaustive evaulation of all potential GAM models. Adding in an option to allow for either using a fixed baseline OR using a new baseline and future for each of the potential GAM models. 

## Preliminaries
library(mgcv)
library(tidyverse)
library(mvtnorm)
library(Smisc)
library(rgeos)
library(akima)
library(sf)
library(viridis)
library(cowplot)
library(corrr)
library(maptools)

## Core functions
# Beta parameterization from mean and variance
est_beta_params_func<- function(mu, var) {
  alpha<- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta<- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model -- combined
loglike_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params){
  
  if(FALSE){
    gam.mod = gam.mod0
    season = season.use
    cand.params = cand.params 
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = c(10, 2, 0)
    nevaV = NULL
    base.update = FALSE
    dist = "beta"
    fix.params = "All"
  }
  
  # NEVA conversion components
  # New baseline and future for specific candidate parameters
  # Get predictor matrix for each of the observations in base.preds and fut.preds
  base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
  fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
  lpmat.base<- predict.gam(gam.mod0, newdata = base.preds.use, type = "lpmatrix")
  lpmat.fut<- predict.gam(gam.mod0, newdata = fut.preds.use, type = "lpmatrix")
  
  # Predictions, combining lpmats and candidate parameters
  # Multiply candidate parameters by predictor matrix
  if(base.update){
    base.prob<- lpmat.base %*% t(cand.params)
    fut.prob<- lpmat.fut %*% t(cand.params)
  } else {
    base.prob<- lpmat.base %*% coef(gam.mod0)
    fut.prob<- lpmat.fut %*% t(cand.params)
  }
  
  # Bin baseline distribution for directional effects and then bin vulnerability
  # Assuming normal distribution...
  if(dist == "norm"){
    # For proposed GAM, get new baseline and future statistics
    bmn<- mean(base.prob, na.rm = T)
    bsd<- sd(base.prob, na.rm = T)
    fmn<- mean(fut.prob, na.rm = T)
    fsd<- sd(fut.prob, na.rm = T)
    
    Dbrks<- dir.brks*bsd+bmn # Directional effect cuts (+/- 1SD), sets bins to integrate future pdf over
    Vbrks<- vuln.brks*bsd/2+bmn # Vulnerability cuts, sets bins to integrate future pdf over
    md<- length(Vbrks)/2+1 # Numeric value (2)
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pnorm(Dbrks, fmn, fsd) # Cumulative distribution function, p(x) <= -1 SD and p(x) <= 1SD from N(fm2, fs2) dist
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pnorm(Vbrks, fmn, fsd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  } 
  
  # Beta distribution
  if(dist == "beta"){
    ilink <- family(gam.mod0)$linkinv
    base.prob<- ilink(base.prob)
    fut.prob<- ilink(fut.prob)
    
    # For proposed GAM, get new baseline and future statistics
    bmn<- mean(base.prob, na.rm = T)
    bsd<- sd(base.prob, na.rm = T)
    fmn<- mean(fut.prob, na.rm = T)
    fsd<- sd(fut.prob, na.rm = T)
    
    # Parameterize beta distribution from mn and sd
    base.params<- est_beta_params_func(mu = bmn, var = bsd^2)
    fut.params<- est_beta_params_func(mu = fmn, var = fsd^2)
    
    Dbrks<- quantile(na.omit(base.prob), prob = c(0.33, 0.67))
    Vbrks<- quantile(na.omit(base.prob), prob = seq(from = 0, to = 1, length.out = 6))
    md<- length(Vbrks)/2+1
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pbeta(Dbrks, fut.params$alpha, fut.params$beta)
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pbeta(Vbrks, fut.params$alpha, fut.params$beta) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  }
  
  if(dist == "norm.diff"){
    # For proposed GAM, get new baseline and future statistics
    diff<- base.prob - fut.prob
    diff.mn<- mean(diff, na.rm = T)
    diff.sd<- sd(diff, na.rm = T)
    
    # Quantiles
    Dbrks<- quantile(na.omit(diff), prob = c(0.25, 0.75))
    Vbrks<- quantile(na.omit(diff), prob = seq(from = 0, to = 1, length.out = 6))
    md<- length(Vbrks)/2+1
    
    ## Directional effect piece: likelihood of NEVA votes given our bmn, bsd, fmn, fsd
    if(!is.null(nevaD)){
      dtmp<- pnorm(Dbrks, diff.mn, diff.sd)
      pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2]) # Per bin probabilities for multinomial distribution
      like.dir<- dmultinom(nevaD, prob = pd, log = TRUE) # Likelihood of observed directional effect NEVA votes given multinomial described by pd per bin probabilities
    }
    
    ## Vulnerability piece
    if(!is.null(nevaV)){
      vtmp<- pnorm(Vbrks, diff.mn, diff.sd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
      vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])
      
      # Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
      vd<- rep(NA, md)
      vd[1]<- vtmp[md]
      
      for(k in 1:(md-1)){
        vd[k+1]<- vtmp[md+k]+vtmp[md-k]
      }
      
      like.vuln<- dmultinom(nevaV, prob = vd, log = TRUE) # Likelihood of observed vulnerability NEVA votes given multinomial described by pd per bin probabilities
    }
    
    if(is.null(nevaD)){
      likelihood<- like.vuln
    } else if(is.null(nevaV)){
      likelihood<- like.dir
    } else {
      likelihood<- like.dir + like.vuln # Likelihood of getting the directional effect votes and the vulnerability votes
    }
  }
  return(likelihood) # Loglikelihood
}

# Probability of our model 
logprior_func<- function(gam.mod0, cand.params, fix.params){
  # Likelihood of the GAM baseline future (fitted curves)
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod<- coef(gam.mod0) # Coefficients for every term/basis function in the model
  vcov.mod<- vcov(gam.mod0) 
  
  # Likelihood of these candidate parameter values, given coef (mean) and sigma (se)
  cand.params.logprior<- dmvnorm(cand.params, mean = coefs.mod, sigma = vcov.mod, log = TRUE)
  return(cand.params.logprior)
}

# Posterior (loglike + logprior)
logposterior_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params) == "numeric"){
    cand.params<- matrix(cand.params, nrow = 1, ncol = length(cand.params), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
  }
  
  loglike.out<- loglike_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params)
  
  # Probability of prior
  logprior.out<- logprior_func(gam.mod0, cand.params)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Prior" = logprior.out, "Posterior" = loglike.out + logprior.out))
}

# Bayes Algorithm
sdm_neva_bayes<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params){
  
  # Calculate likelihood, prior and posterior given candidate
  post<- logposterior_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks, dist, fix.params)
  likes<- cbind(post$Likelihood, post$Prior, post$Posterior)
  return(likes)
}

temp.scale<- function(new.temp, base.temp.mean, base.temp.sd){
  if(is.na(new.temp)){
    temp.scaled<- NA
    return(temp.scaled)
  } else {
    temp.scaled<- (new.temp - base.temp.mean)/base.temp.sd
    return(temp.scaled)
  }
}

#### A few examples
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- c("ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID", "AMERICAN LOBSTER")

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") 

# Training vs. testing
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither"))

# Bottom trawl strata
bstrat<- st_read("./Data/BottomTrawlStrata/BTS_Strata.shp")

# Get names of strata
bstrat.names<- unique(bstrat$STRATA)

# Reduce dataset
dat<- dat[dat$STRATUM %in% bstrat.names,]

# Training data
dat.train.f<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) %>%
  filter(., COMNAME %in% spp.example)
dat.train.f$PRESENCE.ABUNDANCE<- ifelse(dat.train.f$ABUNDANCE > 0, 1, 0)
dat.train.f$WT.ABUNDANCE<- ifelse(dat.train.f$PRESENCE.ABUNDANCE == 0, 1, dat.train.f$ABUNDANCE)

# Temps to rescale other time periods
base.depth.mean.f<- mean(abs(dat.train.f$DEPTH))
base.depth.sd.f<- sd(abs(dat.train.f$DEPTH))
base.shelf.mean.f<- mean(dat.train.f$SHELF_POS)
base.shelf.sd.f<- sd(dat.train.f$SHELF_POS)
base.temp.mean.f<- mean(dat.train.f$SEASONALMU.OISST, na.rm = T)
base.temp.sd.f<- sd(dat.train.f$SEASONALMU.OISST, na.rm = T)
fall.rescale.df<- data.frame(SEASON = "FALL", mean.t = base.temp.mean.f, sd.t = base.temp.sd.f, mean.depth = base.depth.mean.f, sd.depth = base.depth.sd.f, mean.shelf = base.shelf.mean.f, sd.shelf = base.shelf.sd.f)

dat.train.s<- dat %>%
  filter(., TRAIN.TEST == "TRAIN" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST))) %>%
  filter(., COMNAME %in% spp.example)
dat.train.s$PRESENCE.ABUNDANCE<- ifelse(dat.train.s$ABUNDANCE > 0, 1, 0)
dat.train.s$WT.ABUNDANCE<- ifelse(dat.train.s$PRESENCE.ABUNDANCE == 0, 1, dat.train.s$ABUNDANCE)

# Temps to rescale other variables
base.depth.mean.sp<- mean(abs(dat.train.s$DEPTH))
base.depth.sd.sp<- sd(abs(dat.train.s$DEPTH))
base.shelf.mean.sp<- mean(dat.train.s$SHELF_POS)
base.shelf.sd.sp<- sd(dat.train.s$SHELF_POS)
base.temp.mean.sp<- mean(dat.train.s$SEASONALMU.OISST, na.rm = T)
base.temp.sd.sp<- sd(dat.train.s$SEASONALMU.OISST, na.rm = T)
spring.rescale.df<- data.frame(SEASON = "SPRING", mean.t = base.temp.mean.sp, sd.t = base.temp.sd.sp, mean.depth = base.depth.mean.sp, sd.depth = base.depth.sd.sp, mean.shelf = base.shelf.mean.sp, sd.shelf = base.shelf.sd.sp)

## Testing dataframes
dat.test.f<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "FALL") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
dat.test.f$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
dat.test.f$SHELF_POS.Scale<- mapply(temp.scale, dat.test.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
dat.test.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.f$SEASONALMU.OISST, fall.rescale.df$mean.t, fall.rescale.df$sd.t)

dat.test.f<- dat.test.f %>%
  filter(., COMNAME %in% spp.example)
dat.test.f$PRESENCE.ABUNDANCE<- ifelse(dat.test.f$ABUNDANCE > 0, 1, 0)
dat.test.f$WT.ABUNDANCE<- ifelse(dat.test.f$PRESENCE.ABUNDANCE == 0, 1, dat.test.f$ABUNDANCE)

## Testing dataframes
dat.test.s<- dat %>%
  filter(., TRAIN.TEST == "TEST" & SEASON == "SPRING") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
dat.test.s$DEPTH.Scale<- mapply(temp.scale, abs(dat.test.s$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
dat.test.s$SHELF_POS.Scale<- mapply(temp.scale, dat.test.s$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
dat.test.s$SEASONALMU.OISST.Scale<- mapply(temp.scale, dat.test.s$SEASONALMU.OISST, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

dat.test.s<- dat.test.s %>%
  filter(., COMNAME %in% spp.example)
dat.test.s$PRESENCE.ABUNDANCE<- ifelse(dat.test.s$ABUNDANCE > 0, 1, 0)
dat.test.s$WT.ABUNDANCE<- ifelse(dat.test.s$PRESENCE.ABUNDANCE == 0, 1, dat.test.s$ABUNDANCE)

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat.train.f %>%
  bind_rows(., dat.train.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat.test.f %>%
  bind_rows(., dat.test.s) %>%
  group_by(., COMNAME, SEASON) %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

# Running the sdm_neva_bayes function. Easiest to map this function??
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds03232018.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds03232018.rds"
fall.preds<- readRDS(fall.preds)

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.)))

base.preds.sp<- base.preds.sp %>%
  left_join(., spring.rescale.df, by = "SEASON")
base.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
base.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, base.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
base.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.sp$Baseline, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

base.preds.sp<- base.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) 
base.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(base.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
base.preds.f$SHELF_POS.Scale<- mapply(temp.scale, base.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
base.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, base.preds.f$Baseline, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
base.preds.f<- base.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

## Future
fut.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Spring.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("SPRING", nrow(.))) %>%
  left_join(., spring.rescale.df, by = "SEASON")
fut.preds.sp$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.sp$DEPTH), spring.rescale.df$mean.depth, spring.rescale.df$sd.depth)
fut.preds.sp$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.sp$SHELF_POS, spring.rescale.df$mean.shelf, spring.rescale.df$sd.shelf)
fut.preds.sp$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.sp$Spring.2055, spring.rescale.df$mean.t, spring.rescale.df$sd.t)

fut.preds.sp<- fut.preds.sp %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SEASON" = rep("FALL", nrow(.))) %>%
  left_join(., fall.rescale.df, by = "SEASON")
fut.preds.f$DEPTH.Scale<- mapply(temp.scale, abs(fut.preds.f$DEPTH), fall.rescale.df$mean.depth, fall.rescale.df$sd.depth)
fut.preds.f$SHELF_POS.Scale<- mapply(temp.scale, fut.preds.f$SHELF_POS, fall.rescale.df$mean.shelf, fall.rescale.df$sd.shelf)
fut.preds.f$SEASONALMU.OISST.Scale<- mapply(temp.scale, fut.preds.f$Fall.2055, fall.rescale.df$mean.t, fall.rescale.df$sd.t)
fut.preds.f<- fut.preds.f %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

## Also going to want these as prediction dataframes...
# Prediction dataframe
pred.dat.f<- with(fut.preds$Data[[1]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("FALL", 1500)))
rescaled.dat.f<- with(fut.preds$Data[[1]],
                      data.frame("SST" = seq(min(Fall.2055, na.rm = T), max(Fall.2055, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(DEPTH), na.rm = T), max(abs(DEPTH), na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("FALL", 500)))

pred.dat.s<- with(fut.preds$Data[[2]],
                  data.frame(SEASONALMU.OISST.Scale = c(seq(min(SEASONALMU.OISST.Scale, na.rm = T), max(SEASONALMU.OISST.Scale, na.rm = T), length = 500), rep(mean(SEASONALMU.OISST.Scale, na.rm = T), 1000)),
                             DEPTH.Scale = c(rep(mean(DEPTH.Scale, na.rm = T), 500), seq(min(DEPTH.Scale, na.rm = T), max(DEPTH.Scale, na.rm = T), length = 500), rep(mean(DEPTH.Scale, na.rm = T), 500)),
                             SHELF_POS.Scale = c(rep(mean(SHELF_POS.Scale, na.rm = T), 1000), seq(min(SHELF_POS.Scale, na.rm = T), max(SHELF_POS.Scale, na.rm = T), length = 500)), "SEASON" = rep("SPRING", 1500)))
rescaled.dat.s<- with(fut.preds$Data[[2]],
                      data.frame("SST" = seq(min(Spring.2055, na.rm = T), max(Spring.2055, na.rm = T), length = 500),
                                 "Depth" = seq(min(abs(DEPTH), na.rm = T), max(abs(DEPTH), na.rm = T), length = 500),
                                 "Shelf_Pos" = seq(min(SHELF_POS, na.rm = T), max(SHELF_POS, na.rm = T), length = 500),
                                 "Season" = rep("SPRING", 500)))

pred.dat<- bind_rows(pred.dat.f, pred.dat.s)
rescaled.dat<- bind_rows(rescaled.dat.f, rescaled.dat.s)

#### A few examples
## Loop over a few different scenarios...
vote.fac<- 1
dir.scenarios<- data.frame("Dir.Scenario" = c("Negative", "Neutral", "Positive"), "Negative" = c(10, 1, 0)*vote.fac, "Neutral" = c(2, 10, 2)*vote.fac, "Positive" = c(0, 1, 10)*vote.fac)

vuln.scenarios<- data.frame("Vuln.Scenario" = c("Low", "Mod", "High", "VeryHigh"), "Low" = c(10, 1, 0, 0)*vote.fac, "Mod" = c(2, 10, 1, 0)*vote.fac, "High" = c(0, 1, 10, 2)*vote.fac, "VeryHigh" = c(0, 0, 1, 10)*vote.fac)

dir.scenarios.c<- data.frame("Dir.Scenario" = c("Negative", "Neutral", "Positive"), "Negative" = c(10, 1, 0)*vote.fac, "Neutral" = c(2, 10, 2)*vote.fac, "Positive" = c(0, 1, 10)*vote.fac)
vuln.scenarios.c<- data.frame("Vuln.Scenario" = c("Low", "Mod", "High", "VeryHigh"), "Low" = c(10, 1, 1, 0)*vote.fac, "Mod" = c(2, 10, 1, 0)*vote.fac, "High" = c(0, 1, 10, 2)*vote.fac, "VeryHigh" = c(0, 0, 1, 10)*vote.fac)
scenarios.c<- expand.grid(dir.scenarios.c$Dir.Scenario, vuln.scenarios.c$Vuln.Scenario)
names(scenarios.c)<- c("Dir.Scenario", "Vuln.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., dir.scenarios.c, by = "Dir.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., vuln.scenarios.c, by = "Vuln.Scenario")

scenarios.c$Scenario<- paste(scenarios.c$Dir.Scenario, scenarios.c$Vuln.Scenario, sep = ".")

names(dir.scenarios)[1]<- "Scenario"
names(vuln.scenarios)[1]<- "Scenario"

# Set up the loop
scenarios.type<- "Both"
scenarios.run<- scenarios.c
mod.sims<- 10000
base.update.use<- TRUE
dist.use<- "beta"
dir.brks.use<- c(-1, 1)
vuln.brks.use<- c(-3, -2, -1, 1, 2, 3)
out.dir.stem<- "~/Desktop/NormalVoting_041818/"
out.dir.stem<- "~/Desktop/Testing_042018/"

fix.params<- NULL

if(dir.exists(out.dir.stem)){
  print("Directory Exists")
} else {
  dir.create(out.dir.stem)
  print("Directory Created")
}

## Loop
for(i in 1:nrow(dat.train)){
  
  # Get the data and fit the model
  dat.use<- dat.train[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  gam.mod0<- gam(PRESENCE.ABUNDANCE ~ s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  coef.mod0<- coef(gam.mod0)
  vc.mod0<- vcov(gam.mod0)
  saveRDS(gam.mod0, file = paste(out.dir.stem, "gamfit", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), ".rds", sep = ""))
  
  # Draws of candidate models from the posterior
  cand.params.all<- rmvnorm(mod.sims, mean = coef.mod0, sigma = vc.mod0) 
  
  # Fix any of these?
  if(!is.null(fix.params)){
    if(fix.params == "All"){
      # Modify to match all column names of cand.params except intercept
      fix.params.use<- names(coef.mod0)[-which(names(coef.mod0) %in% "(Intercept)")]
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
      
      # Randomly generate intercept values...
      cand.params.all[,1]<- runif(mod.sims, -3, 3)
    } else {
      # Modify to match column names of cand.params
      fix.params.use<- paste(paste("s(", fix.params, ")", sep = ""), ".", seq(from = 1, to = 9, by = 1), sep = "")
      
      # Cand.params.all index...
      cand.params.all.ind<- which(colnames(cand.params.all) %in% fix.params.use, arr.ind = T)
      
      # Coef index
      coef.mod0.fix<- coef.mod0[names(coef.mod0) %in% fix.params.use]
      coef.mod0.fix.mat<- matrix(coef.mod0.fix, nrow = 1, ncol = length(coef.mod0.fix), byrow = T)
      
      # Move over values
      cand.params.all[,cand.params.all.ind]<- coef.mod0.fix.mat[rep(1:nrow(coef.mod0.fix.mat), times = nrow(cand.params.all)),]
    }
  }
  
  for(j in 1:nrow(scenarios.run)){
    scenario.use<- scenarios.run[j,]
    
    out.likes<- array(dim = c(mod.sims, 3))
    out.cands<- array(dim = c(mod.sims, ncol(cand.params.all)))
    
    for(k in 1:nrow(cand.params.all)){
      cand.params<- cand.params.all[k,]
      out.cands[k,]<- cand.params
      
      if(scenarios.type == "Dir") {
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = NULL, base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use, dist = dist.use, fix.params = fix.params)
      } 
      
      if(scenarios.type == "Vuln"){
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = NULL, nevaV = c(scenario.use$Low, scenario.use$Mod, scenario.use$High, scenario.use$VeryHigh), base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use, dist = dist.use, fix.params = fix.params)
      }
      
      if(scenarios.type == "Both"){
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = c(scenario.use$Low, scenario.use$Mod, scenario.use$High, scenario.use$VeryHigh), base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use, dist = dist.use, fix.params = fix.params)
      }
    }
    
    # Bring together results in a list
    mcmc.results<- list(out.likes, out.cands)
    
    # Processing and results plots
    likes.df<- data.frame(mcmc.results[[1]])
    likes.nsamps<- nrow(likes.df)
    
    posts.df<- data.frame(mcmc.results[[2]])
    posts.nsamps<- nrow(posts.df)
    
    # Likes dataframe and plots of mixing
    names(likes.df)<- c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior")
    likes.df$Iteration<- rep(seq(from = 1, to = likes.nsamps, by = 1))
    likes.df<- likes.df %>%
      gather(., Sample, Value, -Iteration)
    likes.df$Sample<- factor(likes.df$Sample, levels = c("Likelihood (NEVA|SDM)", "Prior P(SDM)", "Posterior"))
    
    out.plot<- ggplot(likes.df) +
      geom_line(aes(x = Iteration, y = Value), alpha = 0.25) +
      #scale_color_manual(values = c('#4daf4a', '#377eb8', '#e41a1c')) + 
      ylab("Log Likelihood") + 
      facet_wrap(~Sample, scales = "free") + 
      theme_bw()
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
    
    # Plots
    plot.dat<- posts.df
    colnames(plot.dat)<- names(coef(gam.mod0))
    plot.dat$Iteration<- rep(seq(from = 1, to = posts.nsamps, by = 1))
    plot.dat<- plot.dat %>%
      gather(., Parameter, Value, -Iteration)
    plot.dat$Parameter<- factor(plot.dat$Parameter, levels = names(coef(gam.mod0)))
    
    out.plot<- ggplot(plot.dat) + 
      geom_line(aes(x = Value), stat = "density") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") + 
      theme_bw()
    
    # Add original means
    orig<- data.frame("Parameter" = names(coef(gam.mod0)), "Value" = coef(gam.mod0))
    
    out.plot2<- out.plot +
      geom_vline(data = orig, aes(xintercept = Value), color = "red") +
      facet_wrap(~Parameter, nrow = 4, scales = "free") 
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), "_PostDist.jpg", sep = ""), out.plot2, width = 11, height = 8, dpi = 125, units = "in")
    
    # Maps
    source("~/GitHub/COCA/Code/plot_func.R")
    maps<- plot_func(gam.mod0 = gam.mod0, season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = likes.df[likes.df$Sample == "Posterior",], posts.samps = posts.df)
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), "_Maps.jpg", sep = ""), maps, width = 11, height = 8, dpi = 125, units = "in")
    rm(maps)
    
    # All possible prediction curves, the original smooth, and then the one we have selected...
    # Some prep
    ilink<- family(gam.mod0)$linkinv
    pred.dat.use<- pred.dat[pred.dat$SEASON == season.use, ]
    rescaled.dat.use<- rescaled.dat[rescaled.dat$Season == season.use, ]
    
    # Predictor matrix
    pred.mat<- predict(gam.mod0, newdata = pred.dat.use, type = "lpmatrix")
    
    # Original predictions
    pred0<- ilink(pred.mat %*% coef(gam.mod0))
    
    # Best prediction
    like.posts = likes.df[likes.df$Sample == "Posterior",]
    mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
    pred.best.vec<- as.numeric(posts.df[mod.ind,])
    names(pred.best.vec)<- names(coef(gam.mod0))
    pred.best<- ilink(pred.mat %*% pred.best.vec)
    
    # All predictions
    pred.all<- apply(cand.params.all, 1, function(x) ilink(pred.mat %*% x))
    
    # Plotting
    # SST
    want<- 1:500
    sst.all<- data.frame(pred.all[want,])
    colnames(sst.all)<- paste("Mod.", seq(from = 1, to = ncol(sst.all)), sep = "")
    sst.all<- sst.all %>%
      gather(., "Model", "Pred")
    sst.all$Value<- rep(rescaled.dat.use$SST, length(unique(sst.all$Model)))
    sst.all$Parameter<- rep("SST", nrow(sst.all))
    
    sst.base<- ggplot(sst.all, aes(x = Value, y = Pred, group = Model)) + 
      geom_line(show.legend = F, color = "#bdbdbd") +
      ylim(c(0,1)) 
    
    sst.dat<- data.frame("Parameter" = rep("SST", 1000), "Value" = rep(rescaled.dat.use$SST, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    
    sst.out<- sst.base +
      geom_line(data = sst.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
      scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
      xlab("SST") +
      theme_bw() 
    
    # Depth
    want<- 501:1000
    depth.all<- data.frame(pred.all[want,])
    colnames(depth.all)<- paste("Mod.", seq(from = 1, to = ncol(depth.all)), sep = "")
    depth.all<- depth.all %>%
      gather(., "Model", "Pred")
    depth.all$Value<- rep(rescaled.dat.use$Depth, length(unique(depth.all$Model)))
    depth.all$Parameter<- rep("Depth", nrow(depth.all))
    
    depth.base<- ggplot(depth.all, aes(x = Value, y = Pred, group = Model)) + 
      geom_line(show.legend = F, color = "#bdbdbd") +
      ylim(c(0,1)) 
    
    depth.dat<- data.frame("Parameter" = rep("Depth", 1000), "Value" = rep(rescaled.dat.use$Depth, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    
    depth.out<- depth.base +
      geom_line(data = depth.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
      scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
      xlab("Depth") +
      theme_bw() 
    
    # Along 
    # want<- 1001:1500
    # shelf.all<- data.frame(pred.all[want,])
    # colnames(shelf.all)<- paste("Mod.", seq(from = 1, to = ncol(shelf.all)), sep = "")
    # shelf.all<- shelf.all %>%
    #   gather(., "Model", "Pred")
    # shelf.all$Value<- rep(rescaled.dat.use$Shelf_Pos, length(unique(shelf.all$Model)))
    # shelf.all$Parameter<- rep("Shelf_Pos", nrow(depth.all))
    # 
    # shelf.base<- ggplot(shelf.all, aes(x = Value, y = Pred, group = Model)) + 
    #   geom_line(show.legend = F, color = "#bdbdbd") +
    #   ylim(c(0,1)) 
    # 
    # shelf.dat<- data.frame("Parameter" = rep("Shelf_Pos", 1000), "Value" = rep(rescaled.dat.use$Shelf_Pos, 2), "Pred" = c(pred0[want], pred.best[want]), "Model" = c(rep("SDM", 500), rep("SDM + NEVA", 500)))
    # 
    # shelf.out<- shelf.base +
    #   geom_line(data = shelf.dat, aes(x = Value, y = Pred, group = Model, color = Model)) +
    #   scale_color_manual(name = "Model", values = c('#e41a1c','#377eb8'), labels = c("SDM", "SDM + NEVA")) +
    #   theme_bw()
    # 
    # Arrange and save
    # out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), shelf.out + theme(legend.position="none"), nrow = 1, ncol = 3, align = "hv", scale = 1)
    out<- plot_grid(sst.out + theme(legend.position="none"), depth.out + theme(legend.position="none"), nrow = 1, ncol = 2, align = "hv", scale = 1)
    legend<- get_legend(sst.out)
    out<- plot_grid(out, legend, rel_widths = c(3, .3))
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), "_PredCurv.jpg", sep = ""), out, width = 11, height = 8, dpi = 125, units = "in")
    
    # Update
    print(paste(tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), " is done!", sep = ""))
    
    # Save MCMC results
    saveRDS(mcmc.results, file = paste(out.dir.stem, "mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use$Scenario), ".rds", sep = ""))
  }
}


