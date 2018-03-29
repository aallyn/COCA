#######################
# Data prep, model fitting, projections and port overlap
#######################
# Source all of the R functions
files <- paste("./Code/", list.files(path = "./Code/", pattern = "func"), sep = "")
funcs <- lapply(files, source)

# Prep the modeling dataset -----------------------------------------------
source("./Code/trawl_dat_prep_func.R")
dat.out<- trawl_dat_prep(survdat.file = "Survdat_Nye2016.RData", sst.nc = "EC_sst_1981_2015_OISST-V2-AVHRR_agg_combined.nc")

# Fit Delta Log Normal GAMs -----------------------------------------------
mods<- sdm_fit_func("./Data/model.dat.rds", model.framework = "GAM.DLN", focal.spp = NULL, model.formula = "s(SHELF_POS.Scale, fx = FALSE, bs = 'cs') + s(DEPTH.Scale, fx = FALSE, bs = 'cs') + s(SEASONALMU.OISST.Scale, fx = FALSE, bs = 'cs')", split.type = "transfer", train.start<- "1982-01-01", train.end<- "2010-12-31", test.start<- "2011-01-01", test.end<- "2016-01-01", mod.out = "SST_01172018")

# Make quantitative SDM projections ---------------------------------------
sdm.preds<- sdm_projection_func(mod.dat =  "./Data/dat.fitted.SST_01172018.rds", fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds11202017.rds", spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds11202017.rds", mod.out = "SST_01172018")

# # Bootstrap NEVA data to get qualitative expectations and certainty -------
# neva.out<- neva_boot_func(qual.sensitivity.exposure.path = "~/Dropbox/Andrew/Work/GMRI/COCA/Data/JHareQualitativeDataResults.csv", qual.directional.effect.path = "~/Dropbox/Andrew/Work/GMRI/COCA/Data/JHareDirectionalEffect.csv", boot.run = 10000)
# 
# # Compare quantitative SDM with qualitative NEVA --------------------------
# library(tidyverse)
# qual.dat<- read_csv("./Data/neva.boot.results.csv")
# qual.dat$HareRankVulnerability<- factor(qual.dat$HareRankVulnerability, levels = c("Low", "Moderate", "High", "Very High"))
# qual.dat$Vulnerability.Rank<- factor(qual.dat$Vulnerability.Rank, levels = c("Low", "Moderate", "High", "Very High"))
# qual.dat$HareCertaintyVulnerability<- factor(qual.dat$HareCertaintyVulnerability, levels = c("Low", "Moderate", "High", "Very High"))
# qual.dat$Vulnerability.Certainty<- factor(qual.dat$Vulnerability.Certainty, levels = c("Low", "Moderate", "High", "Very High"))
# 
# # Only want the ones that don't match, vulnerability rank first -- total of nine species
# vuln.mismatch<- qual.dat %>%
#   filter(., HareRankVulnerability != Vulnerability.Rank)
# 
# vuln.mismatch.plot<- ggplot(data = vuln.mismatch, aes(x = factor(HareRankVulnerability), y = factor(Vulnerability.Rank), color = HareCertaintyVulnerability, label = COMNAME)) +
#   geom_text(position=position_jitter(w=.3, h=.3)) +
#   xlab("Hare Rank") +
#   ylab("Our Rank") +
#   theme_bw()
# png(file = "./ExploratoryResults/VulnearbilityMismatch.png", width = 14, height = 8, units = "in", res = 100)
# plot(vuln.mismatch.plot)
# dev.off()
# 
# # Now directional effect
# qual.dat$HareDir.Eff<- factor(qual.dat$HareDir.Eff, levels = c("Negative", "Neutral", "Positive"))
# qual.dat$Dir.Eff<- factor(qual.dat$Dir.Eff, levels = c("Negative", "Neutral", "Positive"))
# qual.dat$HareDir.EffCertainty<- factor(qual.dat$HareDir.EffCertainty, levels = c("Negative", "Neutral", "Positive"))
# qual.dat$Dir.Eff.Certainty<- factor(qual.dat$Dir.Eff.Certainty, levels = c("Negative", "Neutral", "Positive"))
# 
# # Only want the ones that don't match, vulnerability rank first -- None
# direff.mismatch<- qual.dat %>%
#   filter(., HareDir.Eff != Dir.Eff)
# 
# direff.mismatch.plot<- ggplot(data = vuln.mismatch, aes(x = factor(HareRankVulnerability), y = factor(Vulnerability.Rank), color = HareCertaintyVulnerability, label = COMNAME)) +
#   geom_text(position=position_jitter(w=.3, h=.3)) +
#   xlab("Hare Rank") +
#   ylab("Our Rank") +
#   theme_bw()
# png(file = "./ExploratoryResults/VulnearbilityMismatch.png", width = 14, height = 8, units = "in", res = 100)
# plot(vuln.mismatch.plot)
# dev.off()
# 
# # Not too sure what is causing this or why we would have so many correct and only a few wrong. To avoid any issues, just going to use Jon's results from the PLosONE paper.
# 
# # Now, looking at our results versus Jon's. What we are really hoping for is to translate the vulnerability ranks to a numeric value (which has a -1 to 1 range) based on the directional effect bin a species falls into. This function outputs a few different plots to visualize the relationships between our observed changes in model components from the SDM and NEVA directional effect, and vulnerability.
# sdm.vs.neva<- quant_qual_plots_func(quant.dat = "./Data/sdm.projections.rds", qual.dat = "./Data/neva.boot.results.csv")

# Create NEVA future spatially explicit projections ------------------------------------------------------
# Old NEVA 
neva.proj<- neva_projection_func(qual.spp.path = "./Data/Species.DirEff.VulnRank.Diffs.csv", qual.avg.path =  "./Data/AverageAdjustment.DirEff.VulnRank.csv", dir.eff.path = NULL, sdm.projections = "./Data/sdm.projections.SST_01172018.rds", type = "Average", mod.out = "SST_01172018_oldNEVA")

# New NEVA SST 
neva.proj<- neva_projection_func(qual.spp.path = "./Data/JHareQualitativeDataResults.csv", qual.avg.path =  NULL, dir.eff.path = "./Data/JHareDirectionalEffect.csv", sdm.projections = "./Data/sdm.projections.SST_01172018.rds", type = "Distribution", mod.out = "SST_01172018_newNEVA")

# Combined SDM and NEVA projections ---------------------------------------
# Combo OLD NEVA
combo.proj<- combo_projection_func(neva.projs = "./Data/neva.projections.SST_01172018_oldNEVA.rds", sdm.projs = "./Data/sdm.projections.SST_01172018.rds", mod.dat =  "./Data/dat.fitted.SST_01172018.rds", fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds11202017.rds", spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds11202017.rds", type = "Average", mod.out = "SST_01172018_oldNEVA")

# Combo new NEVA SST and BT
combo.proj<- combo_projection_func(neva.projs = "./Data/neva.projections.SST_01172018_newNEVA.rds", sdm.projs = "./Data/sdm.projections.SST_01172018.rds", mod.dat =  "./Data/dat.fitted.SST_01172018.rds", fall.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/fall.rast.preds11202017.rds", spring.preds = "~/Dropbox/Andrew/Work/GMRI/AllData/spring.rast.preds11202017.rds", type = "Distribution", mod.out = "SST_01172018_newNEVA")

# Projection Maps ---------------------------------------------------------
proj.maps<- proj_maps_func(projections.dat = "./Data/AllProjections.rds") # Not sure why this works if you manually run the code inside the function, but it doesn't work as a function.

# Port Differences ---------------------------------------------------------
port.diffs<- port_diffs_func(projections.dat = "./Data/AllProjections.oldNEVA.biomass.rds", port.foots = "./Data/VTR fishing footprints by community and gear type 2011-2015.rds", plot = FALSE, mod.out = "oldNEVA.biomass")

port.diffs<- port_diffs_func(projections.dat = "./Data/AllProjections.newNEVA.biomass.rds", port.foots = "./Data/VTR fishing footprints by community and gear type 2011-2015.rds", plot = FALSE, mod.out = "newNEVA.biomass")

# Port Differences Weighted by Landings ---------------------------------------------------------
port.weighted.diffs<- port_weightedchange_func(port.dat = "./Data/PortAverageChangesinFishAvailability.Biomass.rds", focal.ports = c("PORTLAND.ME", "STONINGTON.ME", "NEW_BEDFORD.MA", "POINT_JUDITH.RI"), landings.dat = "./Data/prop_gear-spp_combos.csv") 


#######################
# Results
#######################
library(stringr)

# Model Fit Diagnostics ---------------------------------------------------
library(tidyverse)

mod.diag<- read_csv("~/GitHub/COCA/Results/GAM.DLN.mod.table.csv")
mod.diag$Dev.Exp.p[mod.diag$Dev.Exp.p == -Inf]<- NA
mod.diag$temp.devexp.prop.p[mod.diag$temp.devexp.prop.p > 1 | mod.diag$temp.devexp.prop.p < 0]<- NA
mod.diag$Dev.Exp.temp.p<- mod.diag$temp.devexp.prop.p * mod.diag$Dev.Exp.p
mod.diag$temp.devexp.prop.b[mod.diag$temp.devexp.prop.b == -Inf]<- NA
mod.diag$Dev.Exp.temp.b<- mod.diag$temp.devexp.prop.b * mod.diag$Dev.Exp.b

# Need to move from wide to long for plotting
mod.diag.l<- mod.diag %>%
  gather(., "Statistic", "Value", -X1, -COMNAME, -SEASON)

# Join with species functional group info
func.groups<- read.csv("~/Dropbox/Andrew/Work/GMRI/COCA/Data/JHareSppFunctionalGroup.csv")
func.groups$COMNAME<- toupper(func.groups$COMNAME)

mod.diag.l<- mod.diag.l %>%
  left_join(., func.groups, by = "COMNAME")

# Deviance Explained Presence and AUC
stats<- c("Dev.Exp.p", "Dev.Exp.temp.p", "AUC")
mod.diag.l$Functional.Group<- factor(mod.diag.l$Functional.Group, levels = c("Diadromous", "Coastal", "Pelagic", "Groundfish", "Invertebrates", "Elasmobranch"))

mod.sub<- mod.diag.l %>%
  filter(., Statistic %in% stats) %>%
  arrange(., SEASON, Functional.Group, Statistic, -Value, COMNAME) %>%
  na.omit()
mod.sub$COMNAME<- factor(mod.sub$COMNAME, levels = unique(mod.sub$COMNAME))

# Plots -- one for each of the functional groups, seasonal panels...
for(i in 1:length(unique(mod.sub$Functional.Group))) {
  group.use<- unique(mod.sub$Functional.Group)[i]
  temp.devexpl.dat<- mod.sub %>%
    filter(., Functional.Group == group.use) 
  plot.out<- ggplot(data = subset(temp.devexpl.dat, Statistic == "Dev.Exp.p"  | Statistic == "Dev.Exp.temp.p"), aes(x = str_to_title(COMNAME), y = Value, fill = Statistic)) +
    geom_col(position = "stack", width = 0.125) +
    xlab("Species") +
    scale_fill_manual(name = "Statistic", values = c("Black", "#377eb8", "#ff7f00")) +
    scale_y_continuous(name = "Deviance Explained", expand = c(0,0.01), limits = c(0, 110)) +
    geom_text(data = subset(temp.devexpl.dat, Statistic == "AUC"), aes(x = str_to_title(COMNAME), y = 104, label = format(Value, nsmall = 0, digits = 2, scientific = FALSE))) +
    facet_wrap(~SEASON, ncol = 6) +
    theme_bw(base_size = 16) +
    coord_flip() +
    ggtitle(paste(group.use)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  png(file = paste("~/GitHub/COCA/Results/", group.use, "DevExpAUC.png", sep = ""), width = 18, height = 14, units = "in", res = 175)
  plot(plot.out)
  dev.off()
}

# Deviance Explained Biomass
stats<- c("Dev.Exp.b", "Dev.Exp.temp.b")
mod.diag.l$Functional.Group<- factor(mod.diag.l$Functional.Group, levels = c("Diadromous", "Coastal", "Pelagic", "Groundfish", "Invertebrates", "Elasmobranch"))

mod.sub<- mod.diag.l %>%
  filter(., Statistic %in% stats) %>%
  arrange(., SEASON, Functional.Group, Statistic, -Value, COMNAME) %>%
  na.omit()
mod.sub$COMNAME<- factor(mod.sub$COMNAME, levels = unique(mod.sub$COMNAME))

# Plots...
for(i in 1:length(unique(mod.sub$Functional.Group))) {
  group.use<- unique(mod.sub$Functional.Group)[i]
  temp.devexpl.dat<- mod.sub %>%
    filter(., Functional.Group == group.use) 
  plot.out<- ggplot(data = temp.devexpl.dat, aes(x = str_to_title(COMNAME), y = Value, fill = Statistic)) +
    geom_col(position = "stack", width = 0.125) +
    xlab("Species") +
    scale_fill_manual(name = "Statistic", values = c("#377eb8", "#ff7f00")) +
    scale_y_continuous(name = "Deviance Explained", expand = c(0,0.01), limits = c(0, 110)) +
    facet_wrap(~SEASON, ncol = 6) +
    theme_bw(base_size = 16) +
    coord_flip() +
    ggtitle(paste(group.use)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  png(file = paste("~/GitHub/COCA/Results/", group.use, "DevExpBiomass.png", sep = ""), width = 18, height = 14, units = "in", res = 175)
  plot(plot.out)
  dev.off()
}

# Projection Maps ---------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(cowplot)
library(gridExtra)

mods<- c("SST_01172018_newNEVA", "SST_01172018_oldNEVA")
focal.spp<- c("AMERICAN LOBSTER", "ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID")

# Presence/Absence
for(h in seq_along(mods)){
  
  spp<- tolower(unique(readRDS("~/GitHub/COCA/Data/AllProjections.biomass.rds")$COMNAME))
  if(!is.null(focal.spp)) {
    spp<- spp[spp %in% tolower(focal.spp)]
  }
  files<- paste("~/GitHub/COCA/Data/", spp, "ProjectionsandMaps.", mods[h], ".rds", sep = "")
  
  for(i in seq_along(files)) {
    
    print(files[i])
    plot.dat<- readRDS(files[i])
    
    # Typical panel plot: Baseline, SDM.Proj, NEVA.Proj, Combo.Proj and then differences PRESENCE
    time.scenario.keep<- c("Base", "SDM.2055.Proj", "NEVA.2055.Proj", "Combo.2055.Proj", "SDM.2055.Diff", "NEVA.2055.Diff", "Combo.2055.Diff")
    plot.dat.sub<- plot.dat %>%
      filter(., Time.Scenario %in% time.scenario.keep)
    
    # Okay, loop over two seasons, one plot for each....
    seasons<- c("FALL", "SPRING")
    
    for(j in seq_along(seasons)) {
      dat.temp<- plot.dat.sub %>%
        filter(., SEASON == seasons[j])
      
      if(any(sapply(dat.temp$Plot, class) == "logical")){
        next
      }
      
      pred.legend<- get_legend(dat.temp$Plot[dat.temp$Time.Scenario == "SDM.2055.Proj"][[1]]) 
      diff.legend<- get_legend(dat.temp$Plot[dat.temp$Time.Scenario == "SDM.2055.Diff"][[1]])
      
      # Bottom-right
      base.plot<- dat.temp$Plot[dat.temp$Time.Scenario == "Base"][[1]] + theme(legend.title = element_blank())
      proj.plots<- grid.arrange(plot_grid(dat.temp$Plot[dat.temp$Time.Scenario == "SDM.2055.Proj"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                          dat.temp$Plot[dat.temp$Time.Scenario == "NEVA.2055.Proj"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                          dat.temp$Plot[dat.temp$Time.Scenario == "Combo.2055.Proj"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()),
                                          pred.legend,
                                          nrow = 1, ncol = 4, scale = 1, labels=c("SDM.2055.Proj", "NEVA.2055.Proj", "Combo.2055.Proj", ""), hjust = -0.85, vjust = 28))
      dev.off()
      diff.plots<- grid.arrange(plot_grid(dat.temp$Plot[dat.temp$Time.Scenario == "SDM.2055.Diff"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                          dat.temp$Plot[dat.temp$Time.Scenario == "NEVA.2055.Diff"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                          dat.temp$Plot[dat.temp$Time.Scenario == "Combo.2055.Diff"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()),
                                          diff.legend,
                                          nrow = 1, ncol = 4, scale = 1, labels=c("SDM.2055.Diff", "NEVA.2055.Diff", "Combo.2055.Diff", ""), hjust = -0.85, vjust = 28))
      dev.off()
      tile.plots<- plot_grid(proj.plots, diff.plots, nrow = 2, scale = 1, align = "hv")
      jpeg(filename = paste("~/GitHub/COCA/Results/Maps/", spp[i], seasons[j], mods[h], ".PresenceMap.jpg", sep = ""), width = 14, height = 10, units = "in", res = 175)
      out<- plot_grid(base.plot, tile.plots, rel_widths = c(0.75,2), labels = c("Baseline", "Projections and Differences"), hjust = c(-2, -1.5), vjust = c(22,6))
      plot(out)
      dev.off()
      rm(base.plot, tile.plots, proj.plots, diff.plots)
    }
    print(paste(spp[i], " is done!", sep = ""))
    rm(plot.dat, plot.dat.sub)
  }
}

# Biomass
for(h in seq_along(mods)){
  
  spp<- tolower(unique(readRDS("~/GitHub/COCA/Data/AllProjections.biomass.rds")$COMNAME))
  if(!is.null(focal.spp)) {
    spp<- spp[spp %in% tolower(focal.spp)]
  }
  files<- paste("~/GitHub/COCA/Data/", spp, "ProjectionsandMaps.", mods[h], ".rds", sep = "")
  
  for(i in seq_along(files)) {

    print(files[i])
    plot.dat<- readRDS(files[i])
    
    # Okay, loop over two seasons, one plot for each....
    seasons<- c("FALL", "SPRING")
    
    # Typical panel plot: Baseline, SDM.Proj, NEVA.Proj, Combo.Proj and then differences BIOMASS
    time.scenario.keep<- c("Base.Biomass", "SDM.2055.Biomass", "Combo.2055.Biomass", "SDM.2055.Biomass.Diff", "Combo.2055.Biomass.Diff")
    plot.dat.sub<- plot.dat %>%
      filter(., Time.Scenario %in% time.scenario.keep)
    
    if(any(sapply(plot.dat.sub$Plot, class) == "logical")){
      next
    }
    
    for(k in seq_along(seasons)) {
      dat.temp.b<- plot.dat.sub %>%
        filter(., SEASON == seasons[k])
      
      if(any(sapply(dat.temp.b$Plot, class) == "logical")){
        next
      }
      
      pred.legend.b<- try(get_legend(dat.temp.b$Plot[dat.temp.b$Time.Scenario == "SDM.2055.Biomass"][[1]]), silent = T) 
      diff.legend.b<- try(get_legend(dat.temp.b$Plot[dat.temp.b$Time.Scenario == "SDM.2055.Biomass.Diff"][[1]]), silent = T)
      
      if(class(pred.legend.b) == "try-error" | class(diff.legend.b) == "try-error") {
        next
      }
      
      # Bottom-right
      base.plot.bio<- dat.temp.b$Plot[dat.temp.b$Time.Scenario == "Base.Biomass"][[1]] + theme(legend.title = element_blank())
      proj.plots.bio<- grid.arrange(plot_grid(dat.temp.b$Plot[dat.temp.b$Time.Scenario == "SDM.2055.Biomass"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                              dat.temp.b$Plot[dat.temp.b$Time.Scenario == "Combo.2055.Biomass"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()),
                                              pred.legend.b,
                                              nrow = 1, ncol = 4, scale = 1, labels=c("SDM.2055.Biomass", "Combo.2055.Biomass", "", ""), hjust = -0.85, vjust = 28))
      dev.off()
      diff.plots.bio<- grid.arrange(plot_grid(dat.temp.b$Plot[dat.temp.b$Time.Scenario == "SDM.2055.Biomass.Diff"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()), 
                                              dat.temp.b$Plot[dat.temp.b$Time.Scenario == "Combo.2055.Biomass.Diff"][[1]] + theme(legend.position = "none", axis.text.x = element_blank(), axis.text.y = element_blank(),  axis.ticks = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()),
                                              diff.legend.b,
                                              nrow = 1, ncol = 4, scale = 1, labels=c("SDM.2055.Biomass.Diff", "Combo.2055.Biomass.Diff", "", ""), hjust = -0.85, vjust = 28))
      dev.off()
      tile.plots.bio<- plot_grid(proj.plots.bio, diff.plots.bio, nrow = 2, scale = 1, align = "hv")
      jpeg(filename = paste("~/GitHub/COCA/Results/Maps/", spp[i], seasons[k], mods[h], ".BiomassMap.jpg", sep = ""), width = 14, height = 10, units = "in", res = 175)
      out<- plot_grid(base.plot.bio, tile.plots.bio, rel_widths = c(0.75,2), labels = c("Baseline.Bio", "Projections and Differences"), hjust = c(-2, -1.5), vjust = c(22,6))
      plot(out)
      dev.off()
      rm(base.plot.bio, tile.plots.bio, proj.plots.bio, diff.plots.bio)
    }
    print(paste(spp[i], " is done!", sep = ""))
    rm(plot.dat, plot.dat.sub)
  }
}

# Synthesis plots ---------------------------------------------------------
###############
### Study area ###
###############
library(sp)
library(sf)
library(maptools)
library(raster)

# Spatial projections
proj.wgs84<- CRS("+init=epsg:4326") #WGS84
proj.utm<- CRS("+init=epsg:2960") #UTM 19

# NELME
nelme<- readShapePoly("~/Dropbox/Andrew/Work/GMRI/AllGIS/nelme.shp")
proj4string(nelme)<- proj.wgs84
nelme.df<- fortify(nelme, LME_NAME)

#Bounds
xlim.use<- c(-77, -64.5)
ylim.use<- c(35.05, 46)

states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")

us <- raster::getData("GADM",country="USA",level=1)
us.states <- us[us$NAME_1 %in% states,]
us.states <- gSimplify(us.states, tol = 0.025, topologyPreserve = TRUE)
canada <- raster::getData("GADM",country="CAN",level=1)
ca.provinces <- canada[canada$NAME_1 %in% provinces,]
ca.provinces <- gSimplify(ca.provinces, tol = 0.025, topologyPreserve = TRUE)

us.states.f<- fortify(us.states, NAME_1)
ca.provinces.f<- fortify(ca.provinces, NAME_1)

ggplot() + 
  #geom_polygon(data = nelme.df, aes(x = long, y = lat, group = group), fill = NA, color = "black") +
  geom_map(data = us.states.f, map = us.states.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  geom_map(data = ca.provinces.f, map = ca.provinces.f,
           aes(map_id = id, group = group),
           fill = "gray65", color = "gray45", size = 0.15) +
  ylim(ylim.use) + ylab("Lat") +
  scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
  coord_fixed(1.3) + 
  theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black")) 

# GoM
gom<- readShapePoly("~/Dropbox/Andrew/Work/GMRI/AllGIS/PhysioRegions_Maine/PhysioRegions_wgs84.shp")
proj4string(gom)<- proj.wgs84
regs.keep<- as.character(gom@data$Region[grepl("Seamount", as.character(gom@data$Region))])
gom<- gom[!as.character(gom@data$Region) %in% regs.keep,]
gom.buff<- raster::buffer(gom, width = 0.1, dissolve = T)

# Southern regions
south<- erase(nelme.buff, gom.buff)

###############
### Shelfwide assessments ###
###############
shelfwide_func<- function(projections = "./Data/AllProjections.rds") {
  library(tidyverse)
  library(maptools)
  library(rgeos)
  library(stringr)
  library(raster)
  source("~/Dropbox/Andrew/Work/GMRI/AllRFunctions/Char.String.Func.R")
  suppressWarnings(sapply(list.files(pattern = "[.]R$", path = "~/Dropbox/Andrew/Work/GMRI/AllRFunctions/", full.names = TRUE), source))
  
  if(FALSE) {
    type = "separate.2055"
    projections<- "./Data/AllProjections.newNEVA.no.biomass.rds"
    mod.out<- "newNEVA.no.biomass"
  }
  
  # QuantQual Projections data
  dat<- readRDS(projections)
  
  # Spatial projections
  proj.wgs84<- CRS("+init=epsg:4326") #WGS84
  proj.utm<- CRS("+init=epsg:2960") #UTM 19
  
  # Only need difference data
  preds.df.sub<- dat %>%
    gather(., "Data.Desc", "Data", -COMNAME, -SEASON, -Coords) %>%
    filter(., grepl("2055.Diff", Data.Desc)) %>%
    arrange(., COMNAME, Data.Desc) 
  
  # NELME
  nelme<- readShapePoly("~/Dropbox/Andrew/Work/GMRI/AllGIS/nelme.shp")
  proj4string(nelme)<- proj.wgs84
  nelme.buff<- raster::buffer(nelme, width = 0.1, dissolve = T)
  
  # GoM
  gom<- readShapePoly("~/Dropbox/Andrew/Work/GMRI/AllGIS/PhysioRegions_Maine/PhysioRegions_wgs84.shp")
  proj4string(gom)<- proj.wgs84
  regs.keep<- as.character(gom@data$Region[grepl("Seamount", as.character(gom@data$Region))])
  gom<- gom[!as.character(gom@data$Region) %in% regs.keep,]
  gom.buff<- raster::buffer(gom, width = 0.1, dissolve = T)
  
  # Southern regions
  south<- erase(nelme.buff, gom.buff)
  
  # Overlay func
  overlay_func<- function(df, region){
    pts.temp<- df
    coordinates(pts.temp)<- ~x+y
    proj4string(pts.temp)<- proj.wgs84
    
    switch(region,
           NELME = mean(df[,3], na.rm = T),
           GOM = mean(data.frame(pts.temp[!is.na(over(pts.temp, gom.buff)),])[,3], na.rm = T),
           South = mean(data.frame(pts.temp[!is.na(over(pts.temp, south)),])[,3], na.rm = T))
  }
  
  preds.df.sub<- preds.df.sub %>%
    mutate(., "NELME.Mean.Change" = purrr::map2(Data, list("NELME"), possibly(overlay_func, NA)),
           "GOM.Mean.Change" = purrr::map2(Data, list("GOM"), possibly(overlay_func, NA)),
           "South.Mean.Change" = purrr::map2(Data, list("South"), possibly(overlay_func, NA)))
  
  
  ## Alright, we are now after a species - scenario - season - region - mean change dataframe...
  res<- preds.df.sub %>%
    dplyr::select(., -Data) %>%
    mutate(., "Scenario" = gsub("[..]", "", gsub('[0-9]', '', Data.Desc)),
           "Time.Period" = as.numeric(regmatches(Data.Desc, gregexpr("[[:digit:]]+", Data.Desc)))) %>%
    dplyr::select(., - Coords, -Data.Desc) %>%
    gather(., "Region", "Mean.Change", -COMNAME, -Scenario, -Time.Period, -SEASON) %>%
    arrange(., COMNAME, SEASON, Scenario, Region) %>%
    unnest() %>%
    data.frame()
  
  res$Scenario<- factor(res$Scenario)
  res$Time.Period<- factor(res$Time.Period)
  res$SEASON<- factor(res$SEASON)
  res$Region<- factor(res$Region, levels = c("South.Mean.Change", "GOM.Mean.Change", "NELME.Mean.Change"))
  
  # Lets add a functional group column...
  # Merge with functional groups....
  func.groups<- read.csv("~/Dropbox/Andrew/Work/GMRI/COCA/Data/JHareSppFunctionalGroup.csv")
  func.groups$COMNAME<- toupper(func.groups$COMNAME)
  
  df.all<- res %>%
    left_join(., func.groups, by = "COMNAME")
  df.all<- df.all[!is.na(df.all$Functional.Group),]
  
  df.all$Functional.Group<- factor(df.all$Functional.Group, levels = c("Diadromous", "Coastal", "Pelagic", "Groundfish", "Invertebrates", "Elasmobranch"))
  df.all<- df.all %>%
    arrange(., Functional.Group, Mean.Change, COMNAME)
  df.all$COMNAME<- factor(df.all$COMNAME, levels = unique(df.all$COMNAME))
  
  # Save it
  write.csv(df.all, file = paste("./Results/", mod.out, ".species.2055.meanchange.csv", sep = ""))
  
  # Plots
  if(type == "separate.2055"){
    plot.type<- c("Negative", "Positive")
    
    for(i in seq_along(plot.type)) {
      plot.use<- plot.type[i]
      
      if(plot.use == "Negative") {
        df<- df.all %>%
          filter(., Mean.Change <= 0)
        
        df<- rbind(df[,], cbind(expand.grid(COMNAME = levels(df.all$COMNAME), SEASON = levels(df.all$Season), Scenario = levels(df.all$Scenario), Time.Period = levels(df.all$Time.Period), Region = levels(df.all$Region), Functional.Group = levels(df.all$Functional.Group), Mean.Change = NA)))
        
      } else {
        df<- df.all %>%
          filter(., Mean.Change > 0)
        
        df<- rbind(df[,], cbind(expand.grid(COMNAME = levels(df.all$COMNAME), Season = levels(df.all$Season), Scenario = levels(df.all$Scenario), Time.Period = levels(df.all$Time.Period), Region = levels(df.all$Region), Functional.Group = levels(df.all$Functional.Group), Mean.Change = NA)))
      }
      
      # NELME
      # Combined P/A difference
      dodge <- position_dodge(width = 0.75)
      
      nelme.combo.df<- df %>%
        filter(., as.character(Region) == "NELME.Mean.Change") %>%
        filter(., grepl("ComboDiff", as.character(Scenario))) %>%
        dplyr::select(., -Scenario)
      
      nelme.combo.df<- nelme.combo.df %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      nelme.combo<- ggplot(data = nelme.combo.df, aes(COMNAME, Mean.Change, fill = SEASON)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        scale_fill_manual("Season", values = c("#fd8d3c", "#2171b5")) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        coord_flip()
      
      png(file = paste("./Results/", plot.use, ".NELME_2055ComboAvgMeanChange.", mod.out, ".png", sep = ""), width = 10, height = 14, units = "in", res = 250)
      plot(nelme.combo)
      dev.off()
      
      # # Biomass
      # nelme.bio.df<- df %>%
      #   filter(., as.character(Region) == "NELME.Mean.Change") %>%
      #   filter(., grepl("Diff.B", as.character(Scenario))) 
      # 
      # nelme.bio.df<- nelme.bio.df %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # nelme.bio<- ggplot(data = nelme.bio.df, aes(COMNAME, Mean.Change, fill = Season)) + 
      #   geom_bar(stat = "identity", width = 0.5, position = dodge) +
      #   scale_fill_manual("Season", values = c("#fd8d3c", "#2171b5")) +
      #   geom_hline(yintercept = 0) +
      #   theme_bw() +
      #   theme(text = element_text(size = 20)) +
      #   coord_flip()
      # 
      # png(file = paste(out.dir, plot.use, "NELME_2055BioAvgMeanChange.png", sep = ""), width = 10, height = 14, units = "in", res = 250)
      # plot(nelme.bio)
      # dev.off() 
      
      # Region specific changes: Presence/Absence
      dat<- df %>%
        filter(., grepl("Combo", as.character(Scenario))) %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      gom.combo<- ggplot(data = subset(dat, as.character(Region) == "GOM.Mean.Change"), aes(COMNAME, Mean.Change, fill = SEASON)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        scale_fill_manual("Season", values = c("#fd8d3c", "#2171b5")) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        coord_flip()
      
      png(file = paste("./Results/", plot.use, ".GoM_2055CombinedAvgMeanChange.", mod.out, ".png", sep = ""), width = 10, height = 14, units = "in", res = 250)
      plot(gom.combo)
      dev.off()
      
      south.combo<- ggplot(data = subset(dat, as.character(Region) == "South.Mean.Change"), aes(COMNAME, Mean.Change, fill = SEASON)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        scale_fill_manual("Season", values = c("#fd8d3c", "#2171b5")) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        theme(text = element_text(size = 20)) +
        coord_flip()
      
      png(file = paste("./Results/", plot.use, ".South_2055CombinedAvgMeanChange.", mod.out, ".png", sep = ""), width = 10, height = 14, units = "in", res = 250)
      plot(south.combo)
      dev.off()
      
      # dat<- df %>%
      #   filter(., Scenario == "Spring.Combo.Diff.") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # spring.reg.combo<- ggplot(data = dat, aes(COMNAME, Mean.Change)) + 
      #   geom_bar(stat = "identity", width = 0.5, position = dodge) +
      #   geom_hline(yintercept = 0) +
      #   theme_bw() +
      #   facet_wrap(~Region) +
      #   coord_flip()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055CombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(spring.reg.combo)
      # dev.off()
      # 
      # # Region specific changes: Biomass
      # dat<- df %>%
      #   filter(., Scenario == "Fall.Diff.B.y") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # fall.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change)) + 
      #   geom_bar(stat = "identity", width = 0.5, position = dodge) +
      #   geom_hline(yintercept = 0) +
      #   theme_bw() +
      #   facet_wrap(~Region) +
      #   coord_flip()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055BioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(fall.reg.bio)
      # dev.off()
      # 
      # dat<- df %>%
      #   filter(., Scenario == "Spring.Diff.B.y") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # spring.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change)) + 
      #   geom_bar(stat = "identity", width = 0.5, position = dodge) +
      #   geom_hline(yintercept = 0) +
      #   theme_bw() +
      #   facet_wrap(~Region) +
      #   coord_flip()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055BioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(spring.reg.bio)
      # dev.off()
      # 
      # ## Functional group means
      # # Presence/absence
      # fall.means<- df %>%
      #   filter(., Scenario == "Fall.Combo.Diff.") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # # compare different sample populations across various temperatures
      # fall.reg.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change)) +
      #   geom_hline(yintercept = 0) +
      #   geom_boxplot() +
      #   facet_wrap(~ Region, nrow = 3) +
      #   coord_flip() +
      #   theme_bw()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055FunGroupCombinedAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(fall.reg.fungroup)
      # dev.off()
      # 
      # spring.means<- df %>%
      #   filter(., Scenario == "Spring.Combo.Diff.") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # 
      # # compare different sample populations across various temperatures
      # spring.reg.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change)) +
      #   geom_hline(yintercept = 0) +
      #   geom_boxplot() +
      #   facet_wrap(~ Region, nrow = 3) +
      #   coord_flip() +
      #   theme_bw()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055FunGroupCombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(spring.reg.fungroup)
      # dev.off()
      # 
      # # Biomass
      # fall.means<- df %>%
      #   filter(., Scenario == "Fall.Diff.B.y") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # 
      # # compare different sample populations across various temperatures
      # fall.reg.bio.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change)) +
      #   geom_hline(yintercept = 0) +
      #   geom_boxplot() +
      #   facet_wrap(~ Region, nrow = 3) +
      #   coord_flip() +
      #   theme_bw()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055FunGroupBioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(fall.reg.fungroup)
      # dev.off()
      # 
      # spring.means<- df %>%
      #   filter(., Scenario == "Spring.Diff.B.y") %>%
      #   group_by(COMNAME) %>%
      #   mutate(., All.NA = all(is.na(Mean.Change))) %>%
      #   filter(., All.NA == FALSE)
      # 
      # # compare different sample populations across various temperatures
      # spring.reg.bio.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change)) +
      #   geom_hline(yintercept = 0) +
      #   geom_boxplot() +
      #   facet_wrap(~ Region, nrow = 3) +
      #   coord_flip() +
      #   theme_bw()
      # 
      # png(file = paste(out.dir, plot.use, "Regions_2055FunGroupBioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      # plot(spring.reg.fungroup)
      # dev.off()
    }
  }
  
  if(type == "Separate") {
    plot.type<- c("Negative", "Positive")
    
    for(i in seq_along(plot.type)) {
      plot.use<- plot.type[i]
      
      if(plot.use == "Negative") {
        df<- df.all %>%
          filter(., Mean.Change <= 0)
        
        df<- rbind(df[,], cbind(expand.grid(COMNAME = levels(df.all$COMNAME), Season = levels(df.all$Season), Scenario = levels(df.all$Scenario), Time.Period = levels(df.all$Time.Period), Region = levels(df.all$Region), Functional.Group = levels(df.all$Functional.Group), Mean.Change = NA)))
        
      } else {
        df<- df.all %>%
          filter(., Mean.Change > 0)
        
        df<- rbind(df[,], cbind(expand.grid(COMNAME = levels(df.all$COMNAME), Season = levels(df.all$Season), Scenario = levels(df.all$Scenario), Time.Period = levels(df.all$Time.Period), Region = levels(df.all$Region), Functional.Group = levels(df.all$Functional.Group), Mean.Change = NA)))
      }
      
      # NELME
      # Combined P/A difference
      dodge <- position_dodge(width = 0.75)
      
      nelme.combo.df<- df %>%
        filter(., as.character(Region) == "NELME.Mean.Change") %>%
        filter(., grepl("Combo.Diff", as.character(Scenario))) %>%
        dplyr::select(., -Scenario)
      
      nelme.combo.df<- nelme.combo.df %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      nelme.combo<- ggplot(data = nelme.combo.df, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Season) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "NELME_ComboAvgMeanChange.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(nelme.combo)
      dev.off()
      
      # Biomass
      nelme.bio.df<- df %>%
        filter(., as.character(Region) == "NELME.Mean.Change") %>%
        filter(., grepl("Diff.B", as.character(Scenario))) 
      
      nelme.bio.df<- nelme.bio.df %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      nelme.bio<- ggplot(data = nelme.bio.df, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Season) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "NELME_BioAvgMeanChange.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(nelme.bio)
      dev.off() 
      
      # Region specific changes: Presence/Absence
      dat<- df %>%
        filter(., Scenario == "Fall.Combo.Diff.") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      fall.reg.combo<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Region) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "Regions_CombinedAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(fall.reg.combo)
      dev.off()
      
      dat<- df %>%
        filter(., Scenario == "Spring.Combo.Diff.") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      spring.reg.combo<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Region) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "Regions_CombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(spring.reg.combo)
      dev.off()
      
      # Region specific changes: Biomass
      dat<- df %>%
        filter(., Scenario == "Fall.Diff.B.y") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      fall.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Region) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "Regions_BioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(fall.reg.bio)
      dev.off()
      
      dat<- df %>%
        filter(., Scenario == "Spring.Diff.B.y") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      spring.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
        geom_bar(stat = "identity", width = 0.5, position = dodge) +
        geom_hline(yintercept = 0) +
        theme_bw() +
        facet_wrap(~Region) +
        coord_flip()
      
      png(file = paste(out.dir, plot.use, "Regions_BioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(spring.reg.bio)
      dev.off()
      
      ## Functional group means
      # Presence/absence
      fall.means<- df %>%
        filter(., Scenario == "Fall.Combo.Diff.") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      # compare different sample populations across various temperatures
      fall.reg.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
        geom_hline(yintercept = 0) +
        geom_boxplot() +
        facet_wrap(~ Region, nrow = 3) +
        coord_flip() +
        theme_bw()
      
      png(file = paste(out.dir, plot.use, "Regions_FunGroupCombinedAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(fall.reg.fungroup)
      dev.off()
      
      spring.means<- df %>%
        filter(., Scenario == "Spring.Combo.Diff.") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      
      # compare different sample populations across various temperatures
      spring.reg.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
        geom_hline(yintercept = 0) +
        geom_boxplot() +
        facet_wrap(~ Region, nrow = 3) +
        coord_flip() +
        theme_bw()
      
      png(file = paste(out.dir, plot.use, "Regions_FunGroupCombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(spring.reg.fungroup)
      dev.off()
      
      # Biomass
      fall.means<- df %>%
        filter(., Scenario == "Fall.Diff.B.y") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      
      # compare different sample populations across various temperatures
      fall.reg.bio.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
        geom_hline(yintercept = 0) +
        geom_boxplot() +
        facet_wrap(~ Region, nrow = 3) +
        coord_flip() +
        theme_bw()
      
      png(file = paste(out.dir, plot.use, "Regions_FunGroupBioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(fall.reg.fungroup)
      dev.off()
      
      spring.means<- df %>%
        filter(., Scenario == "Spring.Diff.B.y") %>%
        group_by(COMNAME) %>%
        mutate(., All.NA = all(is.na(Mean.Change))) %>%
        filter(., All.NA == FALSE)
      
      # compare different sample populations across various temperatures
      spring.reg.bio.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
        geom_hline(yintercept = 0) +
        geom_boxplot() +
        facet_wrap(~ Region, nrow = 3) +
        coord_flip() +
        theme_bw()
      
      png(file = paste(out.dir, plot.use, "Regions_FunGroupBioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
      plot(spring.reg.fungroup)
      dev.off()
    }
  }
  
  if(type == "Both") {
    plot.use<- "Both"
    
    df<- rbind(df.all[,], cbind(expand.grid(COMNAME = levels(df.all$COMNAME), Season = levels(df.all$Season), Scenario = levels(df.all$Scenario), Time.Period = levels(df.all$Time.Period), Region = levels(df.all$Region), Functional.Group = levels(df.all$Functional.Group), Mean.Change = NA)))
    
    # NELME
    # Combined P/A difference
    dodge <- position_dodge(width = 0.75)
    
    nelme.combo.df<- df %>%
      filter(., as.character(Region) == "NELME.Mean.Change") %>%
      filter(., grepl("Combo.Diff", as.character(Scenario))) %>%
      dplyr::select(., -Scenario)
    
    nelme.combo.df<- nelme.combo.df %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    nelme.combo<- ggplot(data = nelme.combo.df, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Season) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "NELME_ComboAvgMeanChange.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(nelme.combo)
    dev.off()
    
    # Biomass
    nelme.bio.df<- df %>%
      filter(., as.character(Region) == "NELME.Mean.Change") %>%
      filter(., grepl("Diff.B", as.character(Scenario))) 
    
    nelme.bio.df<- nelme.bio.df %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    nelme.bio<- ggplot(data = nelme.bio.df, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Season) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "NELME_BioAvgMeanChange.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(nelme.bio)
    dev.off() 
    
    # Region specific changes: Presence/Absence
    dat<- df %>%
      filter(., Scenario == "Fall.Combo.Diff.") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    fall.reg.combo<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Region) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "Regions_CombinedAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(fall.reg.combo)
    dev.off()
    
    dat<- df %>%
      filter(., Scenario == "Spring.Combo.Diff.") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    spring.reg.combo<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Region) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "Regions_CombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(spring.reg.combo)
    dev.off()
    
    # Region specific changes: Biomass
    dat<- df %>%
      filter(., Scenario == "Fall.Diff.B.y") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    fall.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Region) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "Regions_BioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(fall.reg.bio)
    dev.off()
    
    dat<- df %>%
      filter(., Scenario == "Spring.Diff.B.y") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    spring.reg.bio<- ggplot(data = dat, aes(COMNAME, Mean.Change, fill = Time.Period)) + 
      geom_bar(stat = "identity", width = 0.5, position = dodge) +
      geom_hline(yintercept = 0) +
      theme_bw() +
      facet_wrap(~Region) +
      coord_flip()
    
    png(file = paste(out.dir, plot.use, "Regions_BioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(spring.reg.bio)
    dev.off()
    
    ## Functional group means
    # Presence/absence
    fall.means<- df %>%
      filter(., Scenario == "Fall.Combo.Diff.") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    # compare different sample populations across various temperatures
    fall.reg.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      facet_wrap(~ Region, nrow = 3) +
      coord_flip() +
      theme_bw()
    
    png(file = paste(out.dir, plot.use, "Regions_FunGroupCombinedAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(fall.reg.fungroup)
    dev.off()
    
    spring.means<- df %>%
      filter(., Scenario == "Spring.Combo.Diff.") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    
    # compare different sample populations across various temperatures
    spring.reg.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      facet_wrap(~ Region, nrow = 3) +
      coord_flip() +
      theme_bw()
    
    png(file = paste(out.dir, plot.use, "Regions_FunGroupCombinedAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(spring.reg.fungroup)
    dev.off()
    
    # Biomass
    fall.means<- df %>%
      filter(., Scenario == "Fall.Diff.B.y") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    
    # compare different sample populations across various temperatures
    fall.reg.bio.fungroup<- ggplot(fall.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      facet_wrap(~ Region, nrow = 3) +
      coord_flip() +
      theme_bw()
    
    png(file = paste(out.dir, plot.use, "Regions_FunGroupBioAvgMeanChange_Fall.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(fall.reg.fungroup)
    dev.off()
    
    spring.means<- df %>%
      filter(., Scenario == "Spring.Diff.B.y") %>%
      group_by(COMNAME) %>%
      mutate(., All.NA = all(is.na(Mean.Change))) %>%
      filter(., All.NA == FALSE)
    
    # compare different sample populations across various temperatures
    spring.reg.bio.fungroup<- ggplot(spring.means, aes(x = Functional.Group, y = Mean.Change, fill = Time.Period)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      facet_wrap(~ Region, nrow = 3) +
      coord_flip() +
      theme_bw()
    
    png(file = paste(out.dir, plot.use, "Regions_FunGroupBioAvgMeanChange_Spring.png", sep = ""), width = 12, height = 14, units = "in", res = 250)
    plot(spring.reg.fungroup)
    dev.off()
  }
}
###############
### End shelfwide assessments ###
###############

###############
### Combining results from shelfwide assessments ###
###############
oldNEVA.bio<- read.csv(paste("./Results/", "oldNEVA.biomass.species.2055.meanchange.csv", sep = ""))
oldNEVA.bio$Model.Form<- rep("oldNEVA.bio", nrow(oldNEVA.bio))
oldNEVA.no.bio<- read.csv(paste("./Results/", "oldNEVA.no.biomass.species.2055.meanchange.csv", sep = ""))
oldNEVA.no.bio$Model.Form<- rep("oldNEVA.no.bio", nrow(oldNEVA.no.bio))
newNEVA.bio<- read.csv(paste("./Results/", "newNEVA.biomass.species.2055.meanchange.csv", sep = ""))
newNEVA.bio$Model.Form<- rep("newNEVA.bio", nrow(newNEVA.bio))
newNEVA.no.bio<- read.csv(paste("./Results/", "newNEVA.no.biomass.species.2055.meanchange.csv", sep = ""))
newNEVA.no.bio$Model.Form<- rep("newNEVA.no.bio", nrow(newNEVA.no.bio))

mean.change<- bind_rows(oldNEVA.bio, oldNEVA.no.bio, newNEVA.bio, newNEVA.no.bio) 

###############
### End shelfwide assessments ###
###############



###############
### Port change weighted by landings ###
###############
## Port-specific changes
port_weightedchange_func<- function(focal.ports = c("PORTLAND", "STONINGTON", "NEW BEDFORD", "POINT JUDITH"), groups = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12"), out.dir = out.dir) {
  
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  
  if(FALSE) {
    groups<- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12")
    focal.ports<- c("STONINGTON", "PORTLAND", "NEW BEDFORD", "POINT JUDITH")
    out.dir<- "~/Dropbox/Andrew/Work/GMRI/COCA/Spring2017/"
  }
  
  # Read in data from different groups, then bind the rows together into one dataframe for plotting
  for(i in 1:length(groups)) {
    group.use<- groups[i]
    temp<- readRDS(paste(out.dir, group.use, ".PortChangesinAvailability.rds", sep = ""))
    
    if(i == 1) {
      port.dat<- temp
    } else {
      port.dat<- dplyr::bind_rows(port.dat, temp)
    }
  }
  
  # Want a new dataset now that has: Species, Scenario, Port, Difference, Importance and Prop Importance
  port.dat.l<- port.dat %>%
    gather(., "Scenario", "Data", -COMNAME) %>%
    mutate(., "Remove" = ifelse(Data == "", NA, Data)) %>%
    filter(., !is.na(Remove)) %>%
    unnest(Data)
  
  # Cleaning up "Port" description for matching to landings data
  split1<- strsplit(as.character(port.dat.l$Port), "_")
  port.dat.l$Port.Name<- unlist(lapply(split1, "[", 1))
  names(port.dat.l)[3]<- "Port.Gear.Long"
  
  # Need to add importance and x/y...
  landings.dat<- read.csv("~/Dropbox/Andrew/Work/GMRI/COCA/Data/landport1115.csv")
  ports.geo<- read.csv("~/Dropbox/Andrew/Work/GMRI/COCA/Data/geocodedports.csv")
  ports.geo$Port.Name<- gsub(" ", ".", ports.geo$PORT_NAME) 
  ports.geo.sub<- ports.geo %>%
    filter(., Port.Name %in% unique(as.character(port.dat.l$Port.Name))) %>%
    dplyr::select(., Port.Name, lon, lat)
  
  port.dat.l<- port.dat.l %>%
    left_join(., ports.geo.sub, by = "Port.Name")
  
  # From this wide datafile, we want to add the value where species name matches column names of landings.dat and port name matches port row.
  column.indices<- match(gsub(" ", ".", port.dat.l$COMNAME), colnames(landings.dat))
  row.indices<- match(port.dat.l$Port.Name, gsub(" ", ".", landings.dat$PORT))
  landings.dat$filler<- rep(NA, nrow(landings.dat))
  
  # Adjust NA columns (no match in species name)
  column.indices[is.na(column.indices)]<- ncol(landings.dat)
  
  # Add the landings and calculate proportion of landings...by port?
  port.dat.l$landings<- as.numeric(landings.dat[cbind(row.indices, column.indices)])
  
  port.dat.l<- port.dat.l %>%
    dplyr::group_by(., Port.Name, Scenario) %>%
    dplyr::mutate(weighted.importance = round(landings/sum(landings, na.rm = T), 3),
                  weighted.difference = Difference*weighted.importance)
  
  # Write this out
  write.csv(data.frame(port.dat.l), file = paste(out.dir, "SpeciesScenarioPortGearTypeChanges.csv"))
  
  # Calculate port specific sum
  port.aggregated.weighted.difference<- port.dat.l %>%
    dplyr::group_by(Port.Name, lon, lat, Scenario) %>%
    dplyr::summarize(., TotalChange = sum(weighted.difference, na.rm = T)) %>%
    data.frame()
  
  # Write this out and return it
  write.csv(data.frame(port.aggregated.weighted.difference), file = paste(out.dir, "PortAggregatedImportanceWeightedChange.csv"))
  
  # Plots for focal ports
  port.imp.plot<- port.aggregated.weighted.difference %>%
    filter(., grepl("2025", Scenario) | grepl("2040", Scenario) | grepl("2055", Scenario)) %>%
    filter(., Port.Name %in% focal.ports)
  port.imp.plot$Port.Name<- factor(port.imp.plot$Port.Name, levels = focal.ports)
  
  port.imp.plot.out<- ggplot(data = port.imp.plot, aes(x = Port.Name, y = TotalChange, fill = Scenario)) +
    geom_bar(stat = "identity", position = "dodge")  +
    ylab("Proportional Landings Weighted Change in P(presence)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(text = element_text(size = 18))
  
  png(file = paste(out.dir, "FocalPorts_AverageChangeWeightedByImp.png", sep = ""), width = 12, height = 8, units = "in", res = 250)
  plot(port.imp.plot.out)
  dev.off()
}






























