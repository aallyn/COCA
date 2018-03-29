port_diffs_func<- function(projections.dat, port.foots, plot = FALSE, mod.out) {
  
  # Details -----------------------------------------------------------
  # This function calculates changes in probability of presence within a fishing footprint for 7 different gear types.
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "raster", "Hmisc")
  package_check(packages.needed)
  
  if(FALSE) {
    projections.dat = "./Data/AllProjections.oldNEVA.SSTandBT.rds"
    port.foots = "./Data/VTR fishing footprints by community and gear type 2011-2015.rds"
    plot = FALSE
    mod.out<- "oldNEVA.SSTandBT"
  }
  
  # Data
  proj.dat<- readRDS(projections.dat) %>%
    gather(., Scenario, Data, -COMNAME, -SEASON, -Coords)
  all.foot.dat<- readRDS(port.foots)
  ports.names<- all.foot.dat$JGS.COMMUNITY
  
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # Some work on the names
  ports.names<- gsub("\\_(?=[^_]*\\_)", " ", ports.names, perl = TRUE)
  ports.names<- gsub(' +', ' ', ports.names)
  ports.names<- gsub("_", ".", ports.names)
  ports.names<- gsub(" ", "_", ports.names)
  
  # Fishing port footprints -- gear type specific
  gear.types<- all.foot.dat$COST_ID
  gear.types<- ifelse(gear.types == 1, "Dredge",
                      ifelse(gear.types == 2, "Gillnet",
                             ifelse(gear.types == 3, "Longline",
                                    ifelse(gear.types == 4, "Pot/Trap",
                                           ifelse(gear.types == 5, "Purse/Seine",
                                                  ifelse(gear.types == 6, "Trawl", "Other"))))))
  port.foot.names<- paste(ports.names, gear.types, sep = "-")
  ports.all.foots<-all.foot.dat$JGS.NOAA.SAFE.COMMUNITY.GEAR.FOOTPRINTS
  names(ports.all.foots)<- port.foot.names
  
  # Get proportion layer we want for each port-gear type
  port.foots<- unlist(lapply(ports.all.foots, "[", 3))
  
  # Plotting port footprints
  if(plot) {
    
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
    
    ports.plot<- unique(names(port.foots))

    for(i in 1:length(ports.plot)) {
      port.use<- ports.plot[i]
      rast.id<- match(port.use, names(port.foots))
      
      port.rast.dat<- as.data.frame(port.foots[[rast.id]], xy = TRUE)
      names(port.rast.dat)[3]<- "z"
      
      port.rast.dat$z<- ifelse(port.rast.dat$z > 0, 1, NA)
      
      if(all(is.na(port.rast.dat$z))){
        print(port.use)
        next
      }
      
      plot.out<- ggplot() + 
        geom_tile(data = port.rast.dat, aes(x = x, y = y, fill = z), color = "grey", show.legend = F) +
        scale_fill_continuous(na.value = 'white') +
        #scale_fill_manual(name = "Proportion of Kept Catch", values = "blue", na.value = "white") +
        geom_map(data = us.states.f, map = us.states.f,
                 aes(x = long, y = lat, map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        geom_map(data = ca.provinces.f, map = ca.provinces.f,
                 aes(x = long, y = lat, map_id = id, group = group),
                 fill = "gray65", color = "gray45", size = 0.15) +
        ylim(ylim.use) + ylab("Lat") +
        scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
        coord_fixed(1.3) + 
        ggtitle(paste(gsub(".JGS.PROPORTION", "", port.use))) +
        theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black")) 
      port.use<- gsub("/", ".", port.use)
      pdf(file = paste("~/Desktop/PortFootprints/", port.use, "NOAASAFE.pdf", sep = ""), width = 8, height = 12)
      plot(plot.out)
      dev.off()
      print(port.use)
    }
  }
  
  # Calculate proportional weighted change within a footprint
  port.foots.stack<- stack(port.foots)
  
  # Function to calculate fish availability within a port
  fish.avail<- function(data, coords, update) {
    
    if(FALSE){
      row.use<- 1
      data = proj.dat$Data[[row.use]]
      coords = proj.dat$Coords[[row.use]]
      update = paste(proj.dat$COMNAME[[row.use]], proj.dat$SEASON[[row.use]], sep = ".")
    }
    
    dat.sp<- data.frame(coords, "z" = data[,3])
    names(dat.sp)[3]<- "z"
    coordinates(dat.sp)<- ~x+y
    proj4string(dat.sp)<- proj.wgs84
    
    pts.rast<- rasterize(dat.sp, port.foots.stack[[1]], field = "z", fun = mean)
    pts.rast<- raster::resample(pts.rast, port.foots.stack[[1]])
    
    res.mean<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
    res.mean[]<- NA
    res.sd<- vector(length = raster::nlayers(port.foots.stack), mode = "double")
    res.sd[]<- NA
    
    for(i in 1:raster::nlayers(port.foots.stack)) {
      lay.use<- port.foots.stack[[i]]
      
      if(all(is.na(values(lay.use)))){
        res.mean[i]<- NA
        res.sd[i]<- NA
      } else {
        m <- c(0, Inf, 1,  -Inf, 0, 0)
        rclmat <- matrix(m, ncol=3, byrow=TRUE)
        lay.bin<- raster::reclassify(lay.use, rclmat)
        
        # Get coordinates of footprint
        foot.pts <- data.frame(rasterToPoints(lay.bin, function(x) x == 1))
        coordinates(foot.pts)<- ~x+y
        proj4string(foot.pts)<- proj.wgs84
        
        # Okay, we need the projected values at those points and then the proportion of catch at each to use as the weights. 
        # Projected p(presence)
        proj.vals<- raster::extract(pts.rast, foot.pts)
        # Proportion of catch
        proj.weights<- raster::extract(lay.use, foot.pts)
        
        # Weighted mean and sd
        if(all(is.na(proj.vals))) {
          res.mean[i]<- NA
          res.sd[i]<- NA
        } else {
          res.mean[i]<- stats::weighted.mean(proj.vals, proj.weights, na.rm = TRUE)
          res.sd[i]<- sqrt(wtd.var(proj.vals, proj.weights, na.rm = TRUE))
        }
      }
    }
    
    out<- data.frame("Port" = names(port.foots.stack), "Mean.Avail" = res.mean, "SD.Avail" = res.sd)
    return(out)
    print(paste(update, " is done!", sep = ""))
    
  }
  
  preds.df<- proj.dat %>%
    dplyr::mutate("Fish.Availability" = purrr::pmap(list(data = Data, coords = Coords, update = paste(COMNAME, SEASON, sep = ".")), possibly(fish.avail, NA)))
  saveRDS(preds.df, paste("./Data/FishAvailability.", mod.out, ".rds", sep = ""))
  
  ### Now that we've got all the port-fish availability. In the end, mostly interested in the mean change in fish availability for each species. So, for each species, this would mean getting the Fall Proj.P.Diff and Spring Proj.P.Diff combos for each species
  keep.rows<- c("Base", "Combo.2025.Diff", "Combo.2040.Diff", "Combo.2055.Diff")
  
  # Filter and then use spread to get fall and spring season data as their own columns to make averaging easier
  preds.df.sub<- preds.df %>%
    dplyr::filter(., Scenario %in% keep.rows) %>%
    dplyr::select(., -Data, -Coords) %>%
    dplyr::mutate(., "Spread.Col" = paste(SEASON, Scenario, sep = "_")) %>%
    dplyr::select(., -SEASON, -Scenario) %>%
    spread(., Spread.Col, Fish.Availability)
  
  # Yearly average difference
  mean_base_func<- function(fall, spring) {
    mean.df<- data.frame("Port" = fall$Port, "Fall" = fall$Mean.Avail, "Spring" = spring$Mean.Avail)
    data.frame("Port" = mean.df$Port, "Difference" = rowMeans(mean.df[,c(2:3)], na.rm = TRUE))
  }
  
  mean_diff_func<- function(fall, spring) {
    mean.df<- data.frame("Port" = fall$Port, "Fall" = fall$Mean.Avail, "Spring" = spring$Mean.Avail)
    data.frame("Port" = mean.df$Port, "Difference" = rowMeans(mean.df[,c(2:3)], na.rm = TRUE))
  }
  
  ports.df<- preds.df.sub %>%
    dplyr::mutate("Port.Avg.Base" = purrr::map2(FALL_Base, SPRING_Base, possibly(mean_base_func, NA)),
           "Port.Avg.Diff.2025" = purrr::map2(FALL_Combo.2025.Diff, SPRING_Combo.2025.Diff, possibly(mean_diff_func, NA)),
           "Port.Avg.Diff.2040" = purrr::map2(FALL_Combo.2040.Diff, SPRING_Combo.2040.Diff, possibly(mean_diff_func, NA)),
           "Port.Avg.Diff.2055" = purrr::map2(FALL_Combo.2055.Diff, SPRING_Combo.2055.Diff, possibly(mean_diff_func, NA)))
  
  # Save the outputs   
  keep.vec<- c("COMNAME", "Port.Avg.Base", "Port.Avg.Diff.2025", "Port.Avg.Diff.2040", "Port.Avg.Diff.2055")
  
  out<- ports.df %>%
    dplyr::select(., one_of(keep.vec))
  saveRDS(out, file = paste("./Data/PortAverageChangesinFishAvailability.", mod.out, ".rds", sep = ""))
}