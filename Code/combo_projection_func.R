combo_projection_func<- function(neva.projs) {
  
  # Details -----------------------------------------------------------
  # This function will bring in the NEVA projections to make a combined (one future) SDM/NEVA projection
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "gridExtra")
  package_check(packages.needed)
  
  if(FALSE) {
    neva.projs =  "./Data/neva.projections.rds"
    sdm.projs = "./Data/sdm.projections.rds"
  }
  
  # Read in both datasets
  neva.proj<- readRDS(neva.projs) 
  sdm.proj<- readRDS(sdm.projs)
  
  # Combo method function
  qual_quant_combo<- function(sdm, neva, Season.Time) {
    sdm.col<- which(grepl(Season.Time, paste("diff", colnames(sdm), sep = ".")))
    sdm.use<- sdm[,sdm.col]
    total.p.sdm.future<- sum(sdm.use, na.rm = T)
    total.p.neva.future<- sum(neva$qual.projection, na.rm = T)
    total.p.avg<- mean(c(total.p.sdm.future, total.p.neva.future))
    combo<- round(((total.p.avg/total.p.sdm.future)*sdm.use), 3)
    data.frame("x" = sdm$x, "y" = sdm$y, "Projected" = combo)
  }
  
  combo<- neva.proj %>%
    mutate("Combo.Projections" = purrr::pmap(list(sdm = Projections, neva = NEVA.Projections, Season.Time = Season.Time), possibly(qual_quant_combo, NA)))
  
  # Alright, lets get all of these projections together
  sdm.temp<- sdm.proj %>%
    filter(., Component == "Presence") %>%
    dplyr::select(., COMNAME, SEASON, Projections) %>%
    unnest()
  names(sdm.temp)[5:8]<- c("Base", "SDM.2025.Proj", "SDM.2040.Proj", "SDM.2055.Proj")
  
  sdm.base<- sdm.temp %>%
    dplyr::select(., COMNAME, SEASON, Base, SDM.2025.Proj, SDM.2040.Proj, SDM.2055.Proj) %>%
    group_by(., COMNAME, SEASON) %>%
    nest(., Base)
  names(sdm.base)[3]<- "Base"
  
  sdm.2025<- sdm.temp %>%
    dplyr::select(., COMNAME, SEASON, Base, SDM.2025.Proj, SDM.2040.Proj, SDM.2055.Proj) %>%
    group_by(., COMNAME, SEASON) %>%
    nest(., SDM.2025.Proj)
  names(sdm.2025)[3]<- "SDM.2025.Proj"
  
  sdm.2040<- sdm.temp %>%
    dplyr::select(., COMNAME, SEASON, Base, SDM.2025.Proj, SDM.2040.Proj, SDM.2055.Proj) %>%
    group_by(., COMNAME, SEASON) %>%
    nest(., SDM.2040.Proj)
  names(sdm.2040)[3]<- "SDM.2040.Proj"
  
  sdm.2055<- sdm.temp %>%
    dplyr::select(., COMNAME, SEASON, Base, SDM.2025.Proj, SDM.2040.Proj, SDM.2055.Proj) %>%
    group_by(., COMNAME, SEASON) %>%
    nest(., SDM.2055.Proj)
  names(sdm.2055)[3]<- "SDM.2055.Proj"
  
  sdm.all<- sdm.base %>%
    left_join(., sdm.2025, by = c("COMNAME", "SEASON")) %>%
    left_join(., sdm.2040, by = c("COMNAME", "SEASON")) %>%
    left_join(., sdm.2055, by = c("COMNAME", "SEASON"))
  
  neva.temp<- combo %>%
    dplyr::select(., COMNAME, SEASON, Season.Time, NEVA.Projections) %>%
    unnest() %>%
    group_by(., COMNAME, SEASON, Season.Time) %>%
    nest() %>%
    spread(., Season.Time, data)
  names(neva.temp)[3:5]<- c("NEVA.2025.Proj", "NEVA.2040.Proj", "NEVA.2055.Proj")
  
  combo.temp<- combo %>%
    dplyr::select(., COMNAME, SEASON, Season.Time, Combo.Projections) %>%
    unnest() %>%
    group_by(., COMNAME, SEASON, Season.Time) %>%
    nest() %>%
    spread(., Season.Time, data)
  names(combo.temp)[3:5]<- c("Combo.2025.Proj", "Combo.2040.Proj", "Combo.2055.Proj")
  
  all.proj<- sdm.all %>%
    left_join(., neva.temp, by = c("COMNAME", "SEASON")) %>%
    left_join(., combo.temp, by = c("COMNAME", "SEASON"))
  
  # Let's also make a percent difference table
  perc_diff_func<- function(base, new){
    base.use<- data.frame(base)
    names(base.use)[1]<- "Prob"
    new.use<- data.frame(new)
    names(new.use)[1]<- "Prob"
    base.use$Prob[base.use$Prob == 0]<- 0.001
    diff<- round(((new.use$Prob - base.use$Prob)/base.use$Prob)*100,3)
    return(diff)
  }
  
  all.proj.diffs<- all.proj %>%
    gather(., "Scenario", "Proj", -COMNAME, -SEASON, -Base) %>%
    mutate(., "Perc.Diff" = purrr::map2(Base, Proj, possibly(perc_diff_func, NA))) %>% 
    dplyr::select(., -Base, -Proj) %>%
    spread(., Scenario, Perc.Diff)
  
  names(all.proj.diffs)<- c("COMNAME", "SEASON", "Combo.2025.Diff", "Combo.2040.Diff", "Combo.2055.Diff", "NEVA.2025.Diff", "NEVA.2040.Diff", "NEVA.2055.Diff", "SDM.2025.Diff", "SDM.2040.Diff", "SDM.2055.Diff")
  
  # Grab projection coordinates
  coords<- data.frame("x" = neva.proj$Projections[[1]]$x, "y" = neva.proj$Projections[[1]]$y)
  coords.list<- vector("list", 1)
  coords.list[]<- replicate(1, coords, simplify = FALSE)
  coords.merge<- tibble("Coords" = coords.list)
  
  projections.full<- all.proj %>%
    left_join(., all.proj.diffs, by = c("COMNAME", "SEASON")) 
  projections.full$Coords<- rep(coords.merge, nrow(projections.full))

  # Save it all
  saveRDS(projections.full, file = "./Data/AllProjections.rds")
}