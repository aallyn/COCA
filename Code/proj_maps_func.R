proj_maps_func<- function(projections.dat, stats.table) {
  
  # Details -----------------------------------------------------------
  # The function plots projection maps
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "gridExtra", "sf", "sp", "rgeos", "akima", "viridis")
  package_check(packages.needed)
  
  if(FALSE) {
    projections.dat = "./Data/AllProjections.SST_01172018_newNEVA.rds"
    stats.table = "./Results/GAM.DLN.mod.table.SST_01172018.csv"
  }
  
  focal.spp<- c("ACADIAN REDFISH", "AMERICAN LOBSTER", "ATLANTIC COD", "SUMMER FLOUNDER", "LONGFIN SQUID")
  focal.spp<- c("AMERICAN LOBSTER")
  
  projections.dat = "./Data/AllProjections.SST_01172018_newNEVA.rds"
  mod.out = "SST_01172018_newNEVA"
  stats.table = "./Results/GAM.DLN.mod.table.SST_01172018.csv"

  # Read in the data
  proj.dat<- readRDS(projections.dat) 
  
  if(!is.null(focal.spp)){
    proj.dat<- proj.dat %>%
      filter(., COMNAME %in% focal.spp)
  }
  
  # We also will want to make note of the NEVA categories in the plots I think...
  qual.dat<- read_csv("./Data/neva.boot.results.csv") %>%
    dplyr::select(., COMNAME, HareRankVulnerability, HareRankVulnerability.Certainty, HareDirEff, HareDirEff.Certainty) 
  qual.dat<- bind_rows(qual.dat, qual.dat)
  qual.dat$COMNAME<- toupper(qual.dat$COMNAME)
  qual.dat$SEASON<- rep(c("SPRING", "FALL"), each = nrow(qual.dat)/2)

  stats<- read_csv(stats.table) %>%
    dplyr::select(., COMNAME, SEASON, AUC, Corr.b)
  proj.dat<- proj.dat %>%
    left_join(., stats, by = c("COMNAME", "SEASON")) %>%
    left_join(., qual.dat, by = c("COMNAME", "SEASON"))
  
  # Let's move to a long/tidy format, each row will get a map
  proj.dat.l<- proj.dat %>%
    gather(., Time.Scenario, Data, -COMNAME, -SEASON, -Coords, -AUC, -Corr.b, -HareRankVulnerability, -HareRankVulnerability.Certainty, -HareDirEff, -HareDirEff.Certainty)
  
  # Spatial stuff
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # NELME domaine
  nelme<- st_read("./Data/nelme.shp")
  st_crs(nelme)<- proj.wgs84
  nelme.sp<- as(nelme, "Spatial")
  
  #Bounds
  xlim.use<- c(-77, -65)
  ylim.use<- c(35, 45)
  
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
  
  # Ideally, want a bit more fine scale of a map for visual purposes. To do that, we will be doing some interpolation and then cropping the map to the NELME border. This can take a bit of time, so instead of doing it for EVERY dataset, lets do it once here so we can get a vector of the points to keep. Then just pass that into the plotting function.
  coords.df<- data.frame(do.call("cbind", proj.dat.l$Coords[[1]]))
  pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
  pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                          xo=seq(-87.99457, -57.4307, length = 115),
                          yo=seq(22.27352, 48.11657, length = 133))
  pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
  pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
  
  # Okay, now we can run the interpolation and mapping function for each of these datasets, we want a map. So, let's write a mapping function
  map_func<- function(COMNAME, SEASON, Coords, Time.Scenario, Data, nelme, AUC, Corr.b, HareRankVulnerability, HareRankVulnerability.Certainty, HareDirEff, HareDirEff.Certainty) {
    
    if(FALSE){
      row.use<- 14
      COMNAME = proj.dat.l$COMNAME[[row.use]]
      SEASON = proj.dat.l$SEASON[[row.use]]
      Coords = proj.dat.l$Coords[[row.use]]
      Time.Scenario = proj.dat.l$Time.Scenario[[row.use]]
      Data = proj.dat.l$Data[[row.use]]
      nelme = nelme
      AUC = proj.dat.l$AUC[[row.use]]
      Corr.b = proj.dat.l$Corr.b[[row.use]]
      HareRankVulnerability = proj.dat.l$HareRankVulnerability[[row.use]]
      HareRankVulnerability.Certainty = proj.dat.l$HareRankVulnerability.Certainty[[row.use]]
      HareDirEff = proj.dat.l$HareDirEff[[row.use]]
      HareDirEff.Certainty = proj.dat.l$HareDirEff.Certainty[[row.use]]
    }
    
    auc.plot<- round(AUC, 3)
    corr.b.plot<- round(Corr.b, 3)
    harevuln.plot<- substr(HareRankVulnerability, 1, 5)
    harevuln.cert.plot<- substr(HareRankVulnerability.Certainty, 1, 5)
    haredireff.plot<- substr(HareDirEff, 1, 5)
    haredireff.cert.plot<- substr(HareDirEff.Certainty, 1, 5)
    # limits.use<- if(grepl("Biomass", Time.Scenario)){
    #   c(min(Data[,3]), max(Data[,3]))
    #   } else if(min(Data[,3]) >= 0 & max(Data[,3]) <= 1){
    #   c(0, 1)
    #     } else if(min(Data[,3]) >= -1 & max(Data[,3]) <= 1) {
    #   c(-1, 1)
    #       } else if(min(Data[,3]) <= -1 & max(Data[,3]) >= 1) {
    #         c(min(Data[,3]), max(Data[,3]))
    #       }
    
    
    pred.df<- na.omit(data.frame("x" = as.numeric(unlist(Data[,1])), "y" = as.numeric(unlist(Data[,2])), "layer" = as.numeric(unlist(Data[,3]))))
    pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                            xo=seq(-87.99457, -57.4307, length = 115),
                            yo=seq(22.27352, 48.11657, length = 133))
    pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(round(pred.df.interp$z, 2)))
    pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
    
    # Clip to nelme
    pred.df.temp<- pred.sp[which(st_intersects(pred.sp, nelme, sparse = FALSE) == TRUE),]
    coords.keep<- as.data.frame(st_coordinates(pred.df.temp))
    row.names(coords.keep)<- NULL
    pred.df.use<- data.frame(cbind(coords.keep, "z" = as.numeric(pred.df.temp$z)))
    names(pred.df.use)<- c("X", "Y", "z")
    
    # Plot
    if(grepl("Diff", Time.Scenario)) {
      plot.out<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_gradient2(name = paste(tolower(COMNAME), tolower(SEASON), sep = " "), low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white") +
        geom_label(aes(x = -69.5, y = 35.05), label = paste(auc.plot, "/", corr.b.plot, sep = " ")) + 
        geom_label(aes(x = -69.5, y = 36.05), label = paste("Vuln: ", harevuln.plot, "/", harevuln.cert.plot, sep = "")) +
        geom_label(aes(x = -69.5, y = 37.05), label = paste("DirEff: ", haredireff.plot, "/", haredireff.cert.plot, sep = "")) +
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
    } else {
      plot.out<- ggplot() + 
        geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
        scale_fill_viridis(option = "viridis", name = paste(tolower(COMNAME), tolower(SEASON), sep = " "), na.value = "white") +
        geom_label(aes(x = -69.5, y = 35.05), label = paste(auc.plot, "/", corr.b.plot, sep = " ")) + 
        geom_label(aes(x = -69.5, y = 36.05), label = paste("Vuln: ", harevuln.plot, "/", harevuln.cert.plot, sep = "")) +
        geom_label(aes(x = -69.5, y = 37.05), label = paste("DirEff: ", haredireff.plot, "/", haredireff.cert.plot, sep = "")) +
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
    }
    return(plot.out)
  }
  
  proj.maps<-  proj.dat.l %>%
    mutate(., "Plot" = purrr::pmap(list(COMNAME = COMNAME, SEASON = SEASON, Coords = Coords, Time.Scenario = Time.Scenario, Data = Data, nelme = list(nelme), AUC = AUC, Corr.b = Corr.b, HareRankVulnerability = HareRankVulnerability, HareRankVulnerability.Certainty = HareRankVulnerability.Certainty, HareDirEff = HareDirEff, HareDirEff.Certainty = HareDirEff.Certainty), possibly(map_func, NA)))
  
  # Save each species
  spp.vec<- unique(proj.maps$COMNAME)
  
  for(i in 1:length(spp.vec)) {
    spp.use<- spp.vec[i]
    proj.maps.out<- proj.maps %>%
      filter(., COMNAME == spp.use)
    saveRDS(proj.maps.out, file = paste("~/GitHub/COCA/Data/", tolower(spp.use), "ProjectionsandMaps.", mod.out, ".rds", sep = ""))
    print(paste(spp.use, "is saved", sep = ""))
  }
}