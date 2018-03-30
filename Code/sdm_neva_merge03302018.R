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

## Core functions
# NegativeLL function -- the likelihood of NEVA votes given the base and future, which comes out of a potential gam model -- combined
loglike_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  
  if(FALSE){
    gam.mod = mod
    season = "Fall"
    cand.params = rmvnorm(1, mean = coef.mod, sigma = vc.mod) 
    base.preds = base.preds
    fut.preds = fut.preds
    nevaD = c(0, 3, 9)
    nevaV = c(0, 0, 3, 9)
    base.update = FALSE
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
  
  # For proposed GAM, get new baseline and future statistics
  bmn<- mean(base.prob, na.rm = T)
  bsd<- sd(base.prob, na.rm = T)
  fmn<- mean(fut.prob, na.rm = T)
  fsd<- sd(fut.prob, na.rm = T)
  
  # Bin baseline distribution for directional effects and then bin vulnerability
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
  return(likelihood) # Loglikelihood
}

# Probability of our model 
logprior_func<- function(gam.mod0, cand.params){
  # Likelihood of the GAM baseline future (fitted curves)
  # Gather coefficiences and varcov matrix, needed to define the normal dist of each candidate parameter value
  coefs.mod<- coef(gam.mod0) # Coefficients for every term/basis function in the model
  vcov.mod<- vcov(gam.mod0) 
  
  # Likelihood of these candidate parameter values, given coef (mean) and sigma (se)
  cand.params.logprior<- dmvnorm(cand.params, mean = coefs.mod, sigma = vcov.mod, log = TRUE)
  return(cand.params.logprior)
}

# Posterior (loglike + logprior)
logposterior_func<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  # Likelihood of NEVA dir and vuln votes, given GAM baseline future (fitted curves)
  if(class(cand.params) == "numeric"){
    cand.params<- matrix(cand.params, nrow = 1, ncol = length(cand.params), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
  }
  
  loglike.out<- loglike_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks)
  
  # Probability of prior
  logprior.out<- logprior_func(gam.mod0, cand.params)
  
  # Add em all together
  return(data.frame("Likelihood" = loglike.out, "Prior" = logprior.out, "Posterior" = loglike.out + logprior.out))
}

# Bayes Algorithm
sdm_neva_bayes<- function(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks){
  
  # Calculate likelihood, prior and posterior given candidate
  post<- logposterior_func(gam.mod0, season, cand.params, base.preds, fut.preds, nevaD, nevaV, base.update, dir.brks, vuln.brks)
  likes<- cbind(post$Likelihood, post$Prior, post$Posterior)
  return(likes)
}

# Running the sdm_neva_bayes function. Easiest to map this function??
# Need base.preds, fut.preds, nevaD and nevaV
spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds03232018.rds"
spring.preds<- readRDS(spring.preds)

fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds03232018.rds"
fall.preds<- readRDS(fall.preds)

base.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Baseline)),
         "SEASON" = rep("SPRING", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Baseline, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Baseline)),
         "SEASON" = rep("FALL", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

base.preds<- base.preds.f %>%
  bind_rows(., base.preds.sp)

fut.preds.sp<- spring.preds %>%
  dplyr::select(., x, y, Spring.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Spring.2055)),
         "SEASON" = rep("SPRING", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds.f<- fall.preds %>%
  dplyr::select(., x, y, Fall.2055, DEPTH, SHELF_POS) %>%
  mutate(., "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(Fall.2055)),
         "SEASON" = rep("FALL", nrow(.))) %>%
  group_by(., SEASON) %>%
  nest(.key = "Data")

fut.preds<- fut.preds.f %>%
  bind_rows(., fut.preds.sp)

#### A few examples
# Data path
dat.path<- "~/GitHub/COCA/Data/model.dat.rds"

# Fish assessment species
fish.spp<- read.csv("~/GitHub/COCA/Data/Assesmentfishspecies.csv")
spp.example<- c("ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID", "AMERICAN LOBSTER")

# Read it in, filter to one species, do some quick formatting to fit the GAM
dat<- readRDS(dat.path) %>% 
  filter(., SVSPP %in% fish.spp$SVSPP) %>%
  left_join(., fish.spp, by = "SVSPP") %>%
  mutate(., "YEAR" = factor(EST_YEAR, levels = seq(from = min(EST_YEAR), to = max(EST_YEAR), by = 1)),
         "STRATUM.FACTOR" = factor(STRATUM, levels = unique(STRATUM)),
         "SED.SIZE" = factor(SED.SIZE, levels = unique(SED.SIZE)),
         "BIOMASS.WMEAN.ANOM" = as.numeric(BIOMASS.WMEAN.ANOM),
         "BIOMASS.WMEAN.Scale" = as.numeric(scale(BIOMASS.WMEAN)),
         "SHELF_POS.Scale" = as.numeric(scale(SHELF_POS)),
         "DEPTH.Scale" = as.numeric(scale(abs(DEPTH))),
         "SEASONALMU.OISST.Scale" = as.numeric(scale(SEASONALMU.OISST)),
         "BOTTOMTEMP.Scale" = as.numeric(scale(BOTTEMP))) %>%
  filter(., COMNAME %in% spp.example)
dat$PRESENCE.ABUNDANCE<- ifelse(dat$ABUNDANCE > 0, 1, 0)
dat$WT.ABUNDANCE<- ifelse(dat$PRESENCE.ABUNDANCE == 0, 1, dat$ABUNDANCE)

# Look at the catch per tow distributions -- a lot of variability?
dat.mean<- dat %>%
  group_by(COMNAME, SEASON) %>%
  summarize_at("ABUNDANCE", c(mean, sd), na.rm = T)
names(dat.mean)<- c("COMNAME", "SEASON", "MEAN", "SD")
dat.mean <- dat.mean %>%
  mutate(., "ABUNDANCE.UP" = MEAN + 2*SD,
         "ABUNDANCE.LOW" = max(c(MEAN - 2*SD), 0))

dat.summ<- dat %>%
  left_join(., dat.mean, by = c("COMNAME", "SEASON")) %>%
  mutate(., "OUTSIDE" = ifelse(ABUNDANCE > ABUNDANCE.UP, 1, 0))

dat.summ.out<- dat.summ %>%
  group_by(COMNAME, SEASON) %>%
  summarise(SUM.OUT = sum(OUTSIDE))

dat.summ.pres<- dat.summ %>%
  group_by(COMNAME, SEASON) %>%
  summarize(SUM.PRES = sum(PRESENCE.ABUNDANCE))

out<- dat.summ.out %>%
  left_join(., dat.summ.pres, by = c("COMNAME", "SEASON")) %>%
  mutate(., "Prop" = SUM.OUT/SUM.PRES)

### Fitting models
train.start<- "1982-01-01"
train.end<- "2010-12-31"
test.start<- "2011-01-01"
test.end<- "2016-01-01"

dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= train.start & as.Date(dat$DATE) <= train.end, "TRAIN", 
                        ifelse(as.Date(dat$DATE) >= test.start & as.Date(dat$DATE) <= test.end, "TEST", "Neither")) 

# Create nested dataframes, one for testing, one for training
# Training
dat.train<- dat %>%
  group_by(., COMNAME, SEASON, TRAIN.TEST) %>%
  dplyr::filter(., TRAIN.TEST == "TRAIN") %>%
  nest(.key = "TRAIN.DATA") %>%
  arrange(COMNAME)

# Testing
dat.test<- dat %>%
  group_by(., COMNAME, SEASON, TRAIN.TEST) %>%
  dplyr::filter(., TRAIN.TEST == "TEST") %>%
  nest(.key = "TEST.DATA") %>%
  arrange(COMNAME) 

## Loop over a few different scenarios...
vote.fac<- 1
dir.scenarios<- data.frame("Dir.Scenario" = c("Negative", "Neutral", "Positive"), "Negative" = c(10, 1, 0)*vote.fac, "Neutral" = c(2, 10, 2)*vote.fac, "Positive" = c(0, 1, 10)*vote.fac)
nevaV<- NULL

vuln.scenarios<- data.frame("Vuln.Scenario" = c("Low", "Mod", "High", "VeryHigh"), "Low" = c(10, 1, 0, 0)*vote.fac, "Mod" = c(2, 10, 1, 0)*vote.fac, "High" = c(0, 1, 10, 2)*vote.fac, "VeryHigh" = c(0, 0, 1, 10)*vote.fac)
nevaD<- NULL

dir.scenarios.c<- data.frame("Dir.Scenario" = c("Negative", "Positive"), "Negative" = c(10, 0)*vote.fac, "Neutral" = c(2, 2)*vote.fac, "Positive" = c(0, 10)*vote.fac)
vuln.scenarios.c<- data.frame("Vuln.Scenario" = c("Low", "VeryHigh"), "Low" = c(10, 0)*vote.fac, "Mod" = c(2, 0)*vote.fac, "High" = c(0, 2)*vote.fac, "VeryHigh" = c(0, 10)*vote.fac)
scenarios.c<- expand.grid(dir.scenarios.c$Dir.Scenario, vuln.scenarios.c$Vuln.Scenario)
names(scenarios.c)<- c("Dir.Scenario", "Vuln.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., dir.scenarios.c, by = "Dir.Scenario")

scenarios.c<- scenarios.c %>%
  left_join(., vuln.scenarios.c, by = "Vuln.Scenario")

scenarios.c$Scenario<- paste(scenarios.c$Dir.Scenario, scenarios.c$Vuln.Scenario, sep = ".")

# Set up the loop
scenarios.type<- "Dir"
scenarios.run<- dir.scenarios
mod.sims<- 10000
base.update.use<- TRUE
dir.brks.use<- c(-1, 1)
vuln.brks.use<- c(-3, -2, -1, 1, 2, 3)
out.dir.stem<- "~/Desktop/MCMC_Results5/"

## Loop
for(i in 1:nrow(dat.train)){
  
  # Get the data and fit the model
  dat.use<- dat.train[i,]
  season.use<- as.character(dat.use$SEASON[[1]])
  gam.mod0<- gam(PRESENCE.ABUNDANCE ~ s(SHELF_POS.Scale, fx = FALSE, k = 5, bs = 'cs') + s(DEPTH.Scale, k = 5, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, k = 5, fx = FALSE, bs = 'cs'), drop.unused.levels = T, data = dat.use$TRAIN.DATA[[1]], family = binomial(link = logit), select = TRUE)
  coef.mod0<- coef(gam.mod0)
  vc.mod0<- vcov(gam.mod0)
  
  # Draws of candidate models from the posterior
  cand.params.all<- rmvnorm(mod.sims, mean = coef.mod0, sigma = vc.mod0) 
  
  for(j in 1:nrow(scenarios.run)){
    scenario.use<- scenarios.run[j,]
    
    out.likes<- array(dim = c(mod.sims, 3))
    out.cands<- array(dim = c(mod.sims, ncol(cand.params.all)))
    
    for(k in 1:nrow(cand.params.all)){
      cand.params<- cand.params.all[k,]
      out.cands[k,]<- cand.params
      
      if(scenarios.type == "Dir") {
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = NULL, base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use)
      } 
      
      if(scenarios.type == "Vuln"){
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = NULL, nevaV = c(scenario.use$Low, scenario.use$Mod, scenario.use$High, scenario.use$VeryHigh), base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use)
      }
      
      if(scenarios.type == "Both"){
        out.likes[k,]<- sdm_neva_bayes(gam.mod0, season = season.use, cand.params, base.preds, fut.preds, nevaD = c(scenario.use$Negative, scenario.use$Neutral, scenario.use$Positive), nevaV = c(scenario.use$Low, scenario.use$Mod, scenario.use$High, scenario.use$VeryHigh), base.update = base.update.use, dir.brks = dir.brks.use, vuln.brks = vuln.brks.use)
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
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_LikePriorPost.jpg", sep = ""), out.plot, width = 11, height = 8, dpi = 125, units = "in")
    
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
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_PostDist.jpg", sep = ""), out.plot2, width = 11, height = 8, dpi = 125, units = "in")
    
    # Maps
    source("~/GitHub/COCA/Code/plot_func.R")
    maps<- plot_func(gam.mod = gam.mod0, season = season.use, base.preds = base.preds, fut.preds = fut.preds, like.posts = likes.df[likes.df$Sample == "Posterior",], posts.samps = posts.df)
    ggsave(paste(out.dir.stem, tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), "_Maps.jpg", sep = ""), maps, width = 11, height = 8, dpi = 125, units = "in")
    rm(maps)
    
    # Save MCMC results
    saveRDS(mcmc.results, file = paste(out.dir.stem, "mcmc_", tolower(dat.use$COMNAME), "_", tolower(dat.use$SEASON), "_", as.character(scenario.use[1,1]), ".rds", sep = ""))
  }
}
