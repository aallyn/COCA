neva_projection_func<- function(qual.spp.path, qual.avg.path = NULL, sdm.projections, type, dir.eff.path = NULL, mod.out) {
  
  # Details -----------------------------------------------------------
  # The function plots our quantitative projections against Jon Hare's classification
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "gridExtra", "truncnorm")
  package_check(packages.needed)
  
  if(FALSE) {
    qual.spp.path = "./Data/Species.DirEff.VulnRank.Diffs.csv"
    qual.spp.path = "./Data/JHareQualitativeDataResults.csv" # Use for distribution option
    qual.avg.path =  "./Data/AverageAdjustment.DirEff.VulnRank.csv"
    type = "Distribution"
    sdm.projections = "./Data/sdm.projections.SST_01172018.rds"
    dir.eff.path = "./Data/JHareDirectionalEffect.csv"
    mod.out<- "SST_01172018"
  }
  
  if(type == "Average") {
    # Read in qualitative datasets
    # qual.spp.dat has species, directional effect, vulnerability rank and then the average mean p(presence) difference for each time period from baseline 
    qual.spp.dat<- read_csv(qual.spp.path)
    qual.spp.dat$Match.Col<- paste(qual.spp.dat$HareDirEff, qual.spp.dat$HareRankVulnerability, qual.spp.dat$SEASON, qual.spp.dat$Season.Time, sep = ".")
    
    # qual.avg.dat has directional effect, vulnerability rank, season.time and average average mean p(presence) difference ACROSS all species within that directional effect/vulnerability rank bin
    qual.avg.dat<- read_csv(qual.avg.path)
    
    # First step, need to create a range of values for each directional effect that spans across the vulnerability ranks in the qual.avg.dat. This will be based on the Adjustment.Mean values and will basically be the min, 25%, 75% and max from the range. 
    qual.avg.range<- qual.avg.dat %>%
      dplyr::group_by(., HareDirEff, SEASON, Season.Time) %>%
      dplyr::summarize(., "Min" = min(Adjustment.Mean, na.rm = T),
                       "Max" = max(Adjustment.Mean, na.rm = T))
    
    cuts.res<- data.frame(matrix(nrow = 4*3*length(unique(qual.avg.range$Season.Time))*length(unique(qual.avg.range$SEASON)), ncol = 5))
    names(cuts.res)<- c("Cut", "Dir.Eff.Match", "Rank.Match", "Season.Time.Match", "Season")
    start<- 1
    end<- 4
    
    for(i in 1:nrow(qual.avg.range)) {
      temp<- data.frame(qual.avg.range[i,])
      cuts.res[start:end,1]<- switch(as.character(temp$HareDirEff),
                                     "Negative" = seq(from = temp$Max, to = temp$Min, length = 4),
                                     "Neutral" = seq(from = temp$Min, to = temp$Max, length = 4),
                                     "Positive" = seq(from = temp$Min, to = temp$Max, length = 4))
      cuts.res[start:end,2]<- rep(as.character(temp$HareDirEff), 4)
      cuts.res[start:end,3]<- c("Low", "Moderate", "High", "Very High")
      cuts.res[start:end,4]<- rep(as.character(temp$Season.Time), 4)
      cuts.res[start:end, 5]<- rep(as.character(temp$SEASON), 4)
      
      start<- start+4
      end<- end+4
    }
    
    cuts.res$Match.Col<- paste(cuts.res$Dir.Eff.Match, cuts.res$Rank.Match, cuts.res$Season, cuts.res$Season.Time.Match, sep = ".")
    qual.spp.dat$Adjustment.Use<- cuts.res$Cut[match(qual.spp.dat$Match.Col, cuts.res$Match.Col)]
    
    # Now that we have the adjustment value for each species, we need to adjust the current seasonal predicted maps using this adjustment.
    # Get projection data
    sdm.proj<- readRDS(sdm.projections)
    
    # ilink as all predictions for presence are on the link scale and we will need to backtransform them
    ilink<- family(sdm.proj$Model.Fitted[[1]])$linkinv

    # To make this easier, we are only going to be adjusting the values from the presence model component. 
    sdm.proj.red<- sdm.proj %>%
      dplyr::filter(., Component == "Presence") %>%
      dplyr::select(., -Model.Type, -Model.Fitted)
    
    # Similarly, we will only be using adjustments calculated for presence (not biomass). So, we can reduce the qual.spp.dat as well to make things a bit easier to use. 
    qual.spp.dat.sub<- qual.spp.dat %>%
      dplyr::select(., COMNAME, SEASON, Component, HareRankVulnerability, HareDirEff, Season.Time, Adjustment.Use) %>%
      filter(., Component == "Presence" & Season.Time == "diff.2055") 
    
    # Okay, now to make use of the "map"ing functions, we need all of this information in one dataframe. So, here adding the projections list of dataframes to our qual.spp.dat.sub tibble.
    qual.proj<- sdm.proj.red %>%
      dplyr::left_join(., qual.spp.dat.sub, by = c("COMNAME", "SEASON", "Component"))
    
    # Map changes to baseline p(presence) projections
    qual_projections_func<- function(df.use, Adjustment.use, Dir.Eff.use) {
      proj.use<- df.use
      
      # Add in switch for directional effect, positive (add), neutral (don't do anything), negative (subtract)
      qual.projection<- switch(as.character(Dir.Eff.use),
                               "Negative" = ilink(proj.use$Baseline) + Adjustment.use,
                               "Positive" = ilink(proj.use$Baseline) + Adjustment.use,
                               "Neutral" = ilink(proj.use$Baseline))
      qual.projection[qual.projection < 0]<- 0
      qual.projection[qual.projection > 1]<- 1 
      out<- data.frame("x" = proj.use$x, "y" = proj.use$y, "qual.projection" = qual.projection)
      return(out)
    }
    
    qual.proj<- qual.proj %>%
      mutate("NEVA.Proj" = purrr::pmap(list(df.use = Projections, Adjustment.use = Adjustment.Use, Dir.Eff.use = as.character(HareDirEff)), possibly(qual_projections_func, NA)))
    
    qual.proj<- dplyr::select(qual.proj, c("COMNAME", "SEASON", "Component", "NEVA.Proj"))
    
    all.proj<- sdm.proj %>%
      left_join(., qual.proj, by = c("COMNAME", "SEASON", "Component"))
    
    # Save it
    saveRDS(qual.proj, file = paste("./Data/neva.projections.", mod.out, ".rds", sep = ""))
  }
  
  if(type == "Distribution"){
    # Read in SDM data
    sdm.proj<- readRDS(sdm.projections) %>%
      dplyr::filter(., Component == "Presence")
    
    # Get model link for use later on
    ilink<- family(sdm.proj$Model.Fitted[[1]])$linkinv
    
    # Read in qualitative data
    # Directional effect data
    dir.eff.dat<- read.csv(dir.eff.path) 
    
    # Sensitivity and exposure
    qual.dat<- read_csv(qual.spp.path)
    
    # Let's do a bit of processing here. First, we need to translate the votes into probabilities...
    # Directional effect first
    dir.eff.dat$Total.Votes<- rowSums(dir.eff.dat[,3:5])
    
    prob_func<- function(count, total) {
      prob.out<- count/total
      return(prob.out)
    }
    
    dir.eff.dat<- dir.eff.dat %>%
      mutate(., "Neg.Prob" = map2(Negative, Total.Votes, prob_func),
             "Neut.Prob" = map2(Neutral, Total.Votes, prob_func),
             "Pos.Prob" = map2(Positive, Total.Votes, prob_func))
    dir.eff.dat$Species<- toupper(dir.eff.dat$Species)
    names(dir.eff.dat)[1]<- "COMNAME"
    
    # Reduce
    dir.eff.dat<- dir.eff.dat %>%
      dplyr::select(., COMNAME, Neg.Prob, Neut.Prob, Pos.Prob)
    
    # Vulnerability
    # Sum up within attribute category...
    qual.dat.sums<- qual.dat %>%
      group_by(., Species, Functional.Group, Attribute.Category) %>%
      summarise_at(.vars = vars(Low, Moderate, High, Very.High), .funs = sum)
    
    qual.dat.sums<- qual.dat.sums %>%
      mutate(., "Total.Votes" = Low + Moderate + High + Very.High)
    
    # Now, need to translate these to probabilities
    qual.dat.sums<- qual.dat.sums %>%
      mutate(., "Low.Prob" = map2(Low, Total.Votes, prob_func),
             "Moderate.Prob" = map2(Moderate, Total.Votes, prob_func),
             "High.Prob" = map2(High, Total.Votes, prob_func),
             "Very.High.Prob" = map2(Very.High, Total.Votes, prob_func)) %>%
      data.frame()
    qual.dat.sums$Species<- toupper(qual.dat.sums$Species)
    names(qual.dat.sums)[1]<- "COMNAME"
    
    # Reduce
    qual.dat.sums<- qual.dat.sums %>%
      dplyr::select(., COMNAME, Attribute.Category, Low.Prob, Moderate.Prob, High.Prob, Very.High.Prob)
    
    # Separate sensitivity and exposure to make merging easier
    qual.dat.sens<- qual.dat.sums[qual.dat.sums$Attribute.Category == "Sensitivity.Attribute",]
    names(qual.dat.sens)[3:6]<- c("Low.Prob.Sens", "Moderate.Prob.Sens", "High.Prob.Sens", "Very.High.Prob.Sens")
    qual.dat.sens<- qual.dat.sens %>%
      dplyr::select(., -Attribute.Category)
    
    qual.dat.exp<- qual.dat.sums[qual.dat.sums$Attribute.Category == "Exposure.Factor",]
    names(qual.dat.exp)[3:6]<- c("Low.Prob.Exp", "Moderate.Prob.Exp", "High.Prob.Exp", "Very.High.Prob.Exp")
    qual.dat.exp<- qual.dat.exp %>%
      dplyr::select(., -Attribute.Category)
    
    # Now, can we add these pieces of information to the sdm.proj dataframe to then write functions to work on each species (row)?
    sdm.proj<- sdm.proj %>%
      left_join(., dir.eff.dat, by = "COMNAME") %>%
      left_join(., qual.dat.sens, by = "COMNAME", -Attribute.Category) %>%
      left_join(., qual.dat.exp, by = "COMNAME", -Attribute.Category)
    
    # Alright, I think we are ready to get going...
    # A few functions first...
    # SDM prediction normal distribution function
    rnorm_func<- function(samples, mu, sd){
      out.link<- rnorm(n = samples, mean = mu, sd = sd)
      return(out.link)
    }
    
    # SDM prediction distribution on response scale
    sdm_resp_func<- function(sdm.vec.link, ilink.use = ilink){
      out.resp<- ilink.use(sdm.vec.link)
      return(out.resp)
    }
    
    # Inner NEVA sample function
    neva_samp_func_multinorm<- function(samples = 1000, sdm.vec.link, neg.prob, neut.prob, pos.prob, sdm.se, low.prob.sens, moderate.prob.sens, high.prob.sens, very.high.prob.sens, low.prob.exp, moderate.prob.exp, high.prob.exp, very.high.prob.exp, ilink.use = ilink){
      
      if(FALSE){
        samples<- 1000
        sdm.vec.link<- proj.dat.temp$sdm.vec.link[[1]]
        ilink.use<- ilink
        neg.prob<- proj.dat.temp$neg.prob[[1]]
        neut.prob<- proj.dat.temp$neut.prob[[1]]
        pos.prob<- proj.dat.temp$pos.prob[[1]]
        sdm.se<- proj.dat.temp$sdm.se[[1]]
        low.prob.sens<- proj.dat.temp$low.prob.sens[[1]]
        moderate.prob.sens<- proj.dat.temp$moderate.prob.sens[[1]]
        high.prob.sens<- proj.dat.temp$high.prob.sens[[1]]
        very.high.prob.sens<- proj.dat.temp$very.high.prob.sens[[1]]
        low.prob.exp<- proj.dat.temp$low.prob.exp[[1]]
        moderate.prob.exp<- proj.dat.temp$moderate.prob.exp[[1]]
        high.prob.exp<- proj.dat.temp$high.prob.exp[[1]]
        very.high.prob.exp<- proj.dat.temp$very.high.prob.exp[[1]]
      }
      
      ## NEVA means from directional effect
      # Marginal probabilities
      dir.probs<- c(neg.prob, neut.prob, pos.prob)
      dir.samps<- data.frame(t(rmultinom(n = samples, size = 1, prob = dir.probs)))
      
      # Get percentiles of SDM vector
      sdm.quants<- quantile(sdm.vec.link, probs = c(0.15, 0.5, 0.85))
      
      # Marginal distribution sampling index
      dir.samps[,]<- ifelse(dir.samps == 1, TRUE, FALSE)
      
      # Extract percentile value based on marginal distribution sampling index (negative, neutral, positive)
      dir.keep<- apply(dir.samps, 1, function(x) sdm.quants[x])
      
      ## NEVA vulnerability
      # Marginal probabilities for exposure and sensitivity
      probs.exp<- c(low.prob.exp, moderate.prob.exp, high.prob.exp, very.high.prob.exp)
      probs.sens<- c(low.prob.sens, moderate.prob.sens, high.prob.sens, very.high.prob.sens)
      samps.exp<- data.frame(t(rmultinom(n = samples, size = 1, prob = probs.exp)))
      samps.sens<- data.frame(t(rmultinom(n = samples, size = 1, prob = probs.sens)))
      
      # Extract values...
      # Marginal distribution sampling index
      samps.exp[,]<- ifelse(samps.exp == 1, TRUE, FALSE)
      samps.sens[,]<- ifelse(samps.sens == 1, TRUE, FALSE)
      
      # Extract values
      exp.quant<- c(1, 2, 3, 4)
      sens.quant<- c(1, 2, 3, 4)
      samps.exp.keep<- apply(samps.exp, 1, function(x) exp.quant[x])
      samps.sens.keep<- apply(samps.sens, 1, function(x) sens.quant[x])
      
      # Vulnerability = exposure*sensitivity
      samps.vuln<- samps.exp.keep*samps.sens.keep
      
      # Rank based on JH rule
      samps.rank<- cut(samps.vuln, breaks = c(0, 3, 6, 9, 16), labels = c("Low", "Moderate", "High", "Very.High"))
      
      # Now, map each of these ranks to...an SD inflation factor?
      sd.inflation<- data.frame("Rank" = c("Low", "Moderate", "High", "Very.High"), "Adjustment" = seq(from = 1, to = 1.5, length.out = 4))
      neva.vuln.keep<- sd.inflation[match(samps.rank, sd.inflation$Rank),]
      
      neva.vec<- ilink(rnorm(samples, mean = dir.keep, sd = sdm.se*neva.vuln.keep$Adjustment))
      return(neva.vec)
    }
    
    # Outer NEVA function
    neva_dist_func<- function(samples = 1000, Projections, Projections.p.se, Neg.Prob, Neut.Prob, Pos.Prob, Low.Prob.Sens, Moderate.Prob.Sens, High.Prob.Sens, Very.High.Prob.Sens, Low.Prob.Exp, Moderate.Prob.Exp, High.Prob.Exp, Very.High.Prob.Exp){
      
      if(FALSE){
        row.use<- 1
        Projections = sdm.proj$Projections[[row.use]]
        Projections.p.se = sdm.proj$Projections.p.se[[row.use]]
        Neg.Prob = sdm.proj$Neg.Prob[[row.use]]
        Neut.Prob = sdm.proj$Neut.Prob[[row.use]]
        Pos.Prob = sdm.proj$Pos.Prob[[row.use]]
        Low.Prob.Sens = sdm.proj$Low.Prob.Sens[[row.use]]
        Moderate.Prob.Sens = sdm.proj$Moderate.Prob.Sens[[row.use]]
        High.Prob.Sens = sdm.proj$High.Prob.Sens[[row.use]]
        Very.High.Prob.Sens = sdm.proj$Very.High.Prob.Sens[[row.use]]
        Low.Prob.Exp = sdm.proj$Low.Prob.Exp[[row.use]]
        Moderate.Prob.Exp = sdm.proj$Moderate.Prob.Exp[[row.use]]
        High.Prob.Exp = sdm.proj$High.Prob.Exp[[row.use]]
        Very.High.Prob.Exp = sdm.proj$Very.High.Prob.Exp[[row.use]]
        samples = 1000
      }
      
      proj.dat.temp<- as_tibble(data.frame("x" = Projections$x, "y" = Projections$y, "sdm.mu" = Projections$Baseline, "sdm.se" = Projections.p.se$Baseline, "neg.prob" = Neg.Prob, "neut.prob" = Neut.Prob, "pos.prob" = Pos.Prob, "low.prob.sens" = Low.Prob.Sens, "moderate.prob.sens" = Moderate.Prob.Sens, "high.prob.sens" = High.Prob.Sens, "very.high.prob.sens" = Very.High.Prob.Sens, "low.prob.exp" = Low.Prob.Exp, "moderate.prob.exp" = Moderate.Prob.Exp, "high.prob.exp" = High.Prob.Exp, "very.high.prob.exp" = Very.High.Prob.Exp))
      proj.dat.temp<- proj.dat.temp %>%
        mutate(., "sdm.vec.link" = pmap(list(samples = samples, mu = sdm.mu, sd = sdm.se), possibly(rnorm_func, NA)),
               "sdm.vec.resp" = purrr::map(sdm.vec.link, possibly(sdm_resp_func, NA)),
               "NEVA.vec" = pmap(list(samples = samples, sdm.vec.link = sdm.vec.link, neg.prob = neg.prob, neut.prob = neut.prob, pos.prob = pos.prob, sdm.se = sdm.se, low.prob.sens = low.prob.sens, moderate.prob.sens = moderate.prob.sens, high.prob.sens = high.prob.sens, very.high.prob.sens = very.high.prob.sens, low.prob.exp = low.prob.exp, moderate.prob.exp = moderate.prob.exp, high.prob.exp = high.prob.exp, very.high.prob.exp = very.high.prob.exp, ilink.use = list(ilink)), possibly(neva_samp_func_multinorm, NA)),
               "NEVA.mu" = map_dbl(NEVA.vec, median, na.rm = T),
               "NEVA.sd" = map_dbl(NEVA.vec, sd, na.rm = T))
      
      # Reduce kept columns
      neva.proj<- dplyr::select(proj.dat.temp, c("x", "y", "NEVA.mu", "NEVA.sd", "sdm.vec.resp", "NEVA.vec"))
      return(neva.proj)
    }
    
    # Run it
    sdm.proj<- sdm.proj %>%
      dplyr::mutate(., "NEVA.Proj" = pmap(list(Projections = Projections, Projections.p.se = Projections.p.se, Neg.Prob = Neg.Prob, Neut.Prob = Neut.Prob, Pos.Prob = Pos.Prob, Low.Prob.Sens = Low.Prob.Sens, Moderate.Prob.Sens = Moderate.Prob.Sens, High.Prob.Sens = High.Prob.Sens, Very.High.Prob.Sens = Very.High.Prob.Sens, Low.Prob.Exp = Low.Prob.Exp, Moderate.Prob.Exp = Moderate.Prob.Exp, High.Prob.Exp = High.Prob.Exp, Very.High.Prob.Exp = Very.High.Prob.Exp), possibly(neva_dist_func, NA)))
    
    saveRDS(sdm.proj, file = paste("./Data/neva.projections.", mod.out, ".rds", sep = ""))
  }
}