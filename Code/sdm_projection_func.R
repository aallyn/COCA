sdm_projection_func<- function(mod.dat =  "./Data/dat.fitted.rds", fall.preds = "./Data/fall.rast.preds05222017.rds", spring.preds = "./Data/spring.rast.preds05222017.rds") {
  
  # Details -----------------------------------------------------------
  # The function makes quantitative projections
  
  # Preliminaries -----------------------------------------------------------
  # Source helper functions
  source("./Code/HelperFunctions.R")
  
  # Load libraries, using package_check to download if it is missing
  packages.needed<- c("tidyverse", "sp", "raster", "geosphere", "mgcv")
  package_check(packages.needed)
  
  if(FALSE) {
   mod.dat =  "./Data/dat.fitted.rds"
   fall.preds = "./Data/fall.rast.preds05222017.rds"
   spring.preds = "./Data/spring.rast.preds05222017.rds"
  }
  
  # Fitted model object
  dat.fitted<- readRDS(mod.dat)
  
  # Need to get into a long format
  keep.cols<- c("COMNAME", "SEASON", "mod.fitted.p", "mod.fitted.b")
  dat.grouped.l<- dat.fitted %>%
    dplyr::select(., one_of(keep.cols)) %>%
    gather(., Model.Type, Model.Fitted, -COMNAME, -SEASON) %>%
    dplyr::mutate(., Component = ifelse(grepl(".p", Model.Type), "Presence", "Biomass"))
  
  # Prediction covariate datasets
  fall.preds <- readRDS(fall.preds)
  spring.preds <- readRDS(spring.preds)
  
  # Some spatial stuff
  xlim.use<- c(-77, -65)
  ylim.use<- c(35.05, 45.2)
  
  states <- c("Maine", "New Hampshire", "Massachusetts", "Vermont", "New York", "Rhode Island", "Connecticut", "Delaware", "New Jersey", "Maryland", "Pennsylvania", "Virginia", "North Carolina", "South Carolina", "Georgia", "Florida", "District of Columbia", "West Virgina")
  provinces <- c("Ontario", "Québec", "Nova Scotia", "New Brunswick")
  
  us <- raster::getData("GADM",country="USA",level=1)
  canada <- raster::getData("GADM",country="CAN",level=1)
  
  us.states <- us[us$NAME_1 %in% states,]
  us.states.f<- fortify(us.states, NAME_1)
  ca.provinces <- canada[canada$NAME_1 %in% provinces,]
  ca.provinces.f<- fortify(ca.provinces, NAME_1)
  
  ## Prediction datasets 
  # Spatial projections
  proj.wgs84<- CRS("+init=epsg:4326") #WGS84
  proj.utm<- CRS("+init=epsg:2960") #UTM 19
  
  # Add data to dat.grouped.l
  # First, we need to calculate our SHELF_POS covariate
  fall.sp<- fall.preds
  coordinates(fall.sp)<- ~x+y
  proj4string(fall.sp)<- proj.wgs84
  fall.preds.df<- data.frame(fall.sp)
  fall.preds.df$SHELF_POS<- distCosine(fall.sp, cbind(-75, 35), r=6378137)/1000
  
  spring.sp<- spring.preds
  coordinates(spring.sp)<- ~x+y
  proj4string(spring.sp)<- proj.wgs84
  spring.preds.df<- data.frame(spring.sp)
  spring.preds.df$SHELF_POS<- distCosine(spring.sp, cbind(-75, 35), r=6378137)/1000
  
  # Okay, next we can subset each of these
  fall.preds.sub<- fall.preds.df %>%
    dplyr::select(., Baseline, Fall.2025, Fall.2040, Fall.2055, DEPTH, SHELF_POS, SEASON, x, y) %>%
    filter(complete.cases(.))
  spring.preds.sub<- spring.preds.df %>%
    dplyr::select(., Baseline, Spring.2025, Spring.2040, Spring.2055, DEPTH, SHELF_POS, SEASON, x, y) %>%
    filter(complete.cases(.))
  
  # To make use of mapping functions in purrr library, need these projection datasets in our dat.grouped.l tibble
  fall.preds.list<- vector("list", 1)
  fall.preds.list[]<- replicate(1, fall.preds.sub, simplify = FALSE)
  spring.preds.list<- vector("list", 1)
  spring.preds.list[]<- replicate(1, spring.preds.sub, simplify = FALSE)
  preds.data.merge<- tibble("FALL" = fall.preds.list, "SPRING" = spring.preds.list) %>%
    gather(., "SEASON", "Proj.Df")
  
  # Okay, add it to dat.grouped.l
  dat.grouped.l<- dat.grouped.l %>%
    left_join(., preds.data.merge, by = "SEASON")
  
  # Now, a projection function
  project_mod_func<- function(model.fitted, response, proj.dat){
    
    # Make the predictions
    # Get model type
    mod.type<- class(model.fitted)
    
    # Different scenarios
    scenario.vec<- c("Baseline", "2025", "2040", "2055")
    
    if(any(mod.type == "randomForest")) {
      
      # Baseline, 2025, 2040, 2055 -- All stored as a dataframe, but need to rename each to SEASONALMU.OISST. Save as proj.dat
      for(i in seq_along(scenario.vec)){
        # Renaming column as "SEASONALMU.OISST
        temp.dat<- proj.dat
        colnames(temp.dat)[which(grepl(scenario.vec[i], colnames(temp.dat)))]<- "SEASONALMU.OISST"
        
        # Making projections
        proj.temp<- switch(response,
                           "Presence" = round(predict(model.fitted, newdata = temp.dat, type = "prob")[,2], 3),
                           "Biomass" = round(exp(as.numeric(predict(model.fitted, newdata = temp.dat))), 3))
        # Save em
        if(i == 1){
          proj.out<- data.frame("x" = temp.dat$x, "y" = temp.dat$y, "z" = as.numeric(proj.temp))
          names(proj.out)[3]<- scenario.vec[i]
        } else {
          proj.out<- cbind(proj.out, as.numeric(proj.temp))
          names(proj.out)[2+i]<- scenario.vec[i]
        }
      }
      return(proj.out)
    }
    
    if(any(mod.type == "gam")){
      # Baseline, 2025, 2040, 2055 -- All stored as a dataframe, but need to rename each to SEASONALMU.OISST. Save as proj.dat
      for(i in seq_along(scenario.vec)){
        # Renaming column as "SEASONALMU.OISST
        temp.dat<- proj.dat
        colnames(temp.dat)[which(grepl(scenario.vec[i], colnames(temp.dat)))]<- "SEASONALMU.OISST"
        
        # Making projections
        proj.temp<- switch(response,
                           "Presence" = round(as.numeric(predict.gam(model.fitted, newdata = temp.dat, type = "response", se.fit = TRUE)$fit), 3),
                           "Biomass" = round(exp(as.numeric(predict.gam(model.fitted, newdata = temp.dat, type = "response", se.fit = TRUE)$fit)), 3))
        
        # Save em
        if(i == 1){
          proj.out<- data.frame("x" = temp.dat$x, "y" = temp.dat$y, "z" = as.numeric(proj.temp))
          names(proj.out)[3]<- scenario.vec[i]
        } else {
          proj.out<- cbind(proj.out, as.numeric(proj.temp))
          names(proj.out)[2+i]<- scenario.vec[i]
        }
      }
      return(proj.out)
    }
  }
  
  # Alright, run it on each row of our dat.grouped.l dataframe
  dat.grouped.l<- dat.grouped.l %>%
    mutate(., "Projections" = purrr::pmap(list(model.fitted = Model.Fitted, response = Component, proj.dat = Proj.Df), project_mod_func))
  
  # Save SDM projections
  saveRDS(dat.grouped.l, file = "./Data/sdm.projections.rds")
}