port_weightedchange_func<- function(port.dat, focal.ports = c("PORTLAND.ME", "STONINGTON.ME", "NEW_BEDFORD.MA", "POINT_JUDITH.RI"), landings.dat = "./Data/prop_gear-spp_combos.csv") {
  
  # Details -----------------------------------------------------------
  # This function calculates 
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "viridis", "stringr")
  package_check(packages.needed)
  
  if(FALSE) {
   port.dat = "./Data/PortAverageChangesinFishAvailability.Biomass.rds"
   focal.ports = c("PORTLAND.ME", "STONINGTON.ME", "NEW_BEDFORD.MA", "POINT_JUDITH.RI")
   landings.dat = "./Data/prop_gear-spp_combos.csv"
  }

  # Want a new dataset now that has: Species, Scenario, Port, Difference, Importance and Prop Importance
  port.dat.l<- readRDS(port.dat) %>%
    gather(., "Scenario", "Data", -COMNAME) %>%
    unnest(Data) 
  
  # Need to add importance and x/y...
  landings.dat<- read_csv(landings.dat)
  
  # Geocode the landings
  library(ggmap)
  geo.name<- paste(landings.dat$Port_Name, "_", landings.dat$State, sep = "")
  xy.locals<- geocode(tolower(as.character(unique(geo.name))))
  xy.locals.bind<- bind_cols(data.frame(as.character(unique(landings.dat$Port_Name_State))), xy.locals)
  names(xy.locals.bind)[1]<- "Port_Name_State"
  landings.dat<- landings.dat %>%
    dplyr::mutate(., "lon" = xy.locals.bind$lon[match(landings.dat$Port_Name_State, xy.locals.bind$Port_Name_State)],
           "lat" = xy.locals.bind$lat[match(landings.dat$Port_Name_State, xy.locals.bind$Port_Name_State)])
  
  # Need to make some modifications so that we can move the landings.dat info over to our port.dat.l. In port.dat.l, remove JGS.PROPORTION
  port.dat.l$Port<- gsub(".JGS.SAFE.PROPORTION", "", port.dat.l$Port)
  
  # Add gear type to landings.dat Port_Name_State
  landings.dat$port.dat.match<- paste(landings.dat$Port_Name_State, ".", landings.dat$`Gear Type`, sep = "")
  
  # Reduce landings.dat a bit to make joining easier...
  keep.cols<- c("Jon_Species", "port.dat.match", "Prop_Volume", "lon", "lat")
  landings.sub<- landings.dat %>%
    dplyr::select(., one_of(keep.cols))
  names(landings.sub)[c(1, 2)]<- c("COMNAME", "Port")
  
  port.dat.l<- port.dat.l %>%
    left_join(., landings.sub, by = c("COMNAME", "Port"))
  
  # Calculate Difference weighted by prop_volume
  port.dat.l<- port.dat.l %>%
    dplyr::mutate(weighted.difference = Difference*Prop_Volume)
  
  # Write this out
  write.csv(data.frame(port.dat.l), file = "./Results/SpeciesScenarioPortGearTypeChanges.Biomass.csv")
  
  # Calculate port specific sum
  port.dat.l<- port.dat.l %>%
    dplyr::mutate(., "Port.Group" = as.vector(sapply(Port, FUN = function(x) paste(strsplit(x,"[.]")[[1]][1:2],collapse="."))))
  port.aggregated.weighted.difference<- port.dat.l %>%
    dplyr::group_by(Scenario, Port.Group) %>%
    dplyr::summarize(., TotalChange = sum(weighted.difference, na.rm = T)) %>%
    data.frame()
  
  # Write this out and return it
  write.csv(data.frame(port.aggregated.weighted.difference), file = "./Results/PortAggregatedImportanceWeightedChange.Biomass.csv")
  
  # Plots for focal ports
  port.imp.plot<- port.aggregated.weighted.difference %>%
    filter(., grepl("2025", Scenario) | grepl("2040", Scenario) | grepl("2055", Scenario)) %>%
    filter(., Port.Group %in% focal.ports) %>%
    mutate(., Port.Name = factor(Port.Group, levels = focal.ports))

  port.imp.plot.out<- ggplot(data = port.imp.plot, aes(x = Port.Name, y = TotalChange, fill = Scenario)) +
    geom_bar(stat = "identity", position = "dodge")  +
    ylab("Proportional Landings \n Weighted Change in P(presence)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    theme(text = element_text(size = 18))
  
  png(file = "./Results/FocalPorts_AverageChangeWeightedByImp.Biomass.png", width = 12, height = 8, units = "in", res = 250)
  plot(port.imp.plot.out)
  dev.off()
}