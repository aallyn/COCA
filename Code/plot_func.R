plot_func<- function(gam.mod0, season, base.preds, fut.preds, like.posts, posts.samps) {
  if(FALSE){
    gam.mod0 = gam.mod0
    season = season.use
    base.preds = base.preds
    fut.preds = fut.preds
    likes.posts = likes.df[likes.df$Sample == "Posterior",]
    posts.samps = posts.df
  }
  
  # SDM
  base.map<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = predict.gam(gam.mod0, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "response"))
  fut.map<- data.frame("x" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$x, "y" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$y, "pred" = predict.gam(gam.mod0, newdata = fut.preds$Data[[match(season, fut.preds$SEASON)]], type = "response"))
  
  # SDM+NEVA
  # Get maximimum value and extract row from posts.samps
  mod.ind<- like.posts$Iteration[which.max(like.posts$Value)]
  best.fit<- posts.samps[mod.ind,]
  best.fit.mat<- matrix(as.numeric(best.fit), nrow = 1, ncol = length(best.fit), byrow = T, dimnames = list(NULL, names(coef(gam.mod0))))
  
  # Make predictions with these values
  lpmat.base<- predict.gam(gam.mod0, newdata = base.preds$Data[[match(season, base.preds$SEASON)]], type = "lpmatrix")
  combo.map.base<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y, "pred" = as.numeric(lpmat.base %*% t(best.fit.mat)))
  
  lpmat.fut<- predict.gam(gam.mod0, newdata = fut.preds$Data[[match(season, fut.preds$SEASON)]], type = "lpmatrix")
  combo.map.fut<- data.frame("x" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$x, "y" = fut.preds$Data[[match(season, fut.preds$SEASON)]]$y, "pred" = as.numeric(lpmat.fut %*% t(best.fit)))
  
  # Response scale
  ilink <- family(gam.mod0)$linkinv
  combo.map.base$pred<- ilink(combo.map.base$pred)
  combo.map.fut$pred<- ilink(combo.map.fut$pred)
  
  # Maps
  # Spatial projections
  proj.wgs84<- "+init=epsg:4326" #WGS84
  proj.utm<- "+init=epsg:2960" #UTM 19
  
  # NELME domaine
  nelme<- st_read("~/GitHub/COCA/Data/NELME_clipped.shp")
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
  coords.df<- data.frame("x" = base.preds$Data[[match(season, base.preds$SEASON)]]$x, "y" = base.preds$Data[[match(season, base.preds$SEASON)]]$y)
  pred.df<- na.omit(data.frame("x" = coords.df$x, "y" = coords.df$y, "layer" = rep(0, length(coords.df$x))))
  pred.df.interp<- interp(pred.df[,1], pred.df[,2], pred.df[,3], duplicate = "mean", extrap = TRUE,
                          xo=seq(-87.99457, -57.4307, length = 115),
                          yo=seq(22.27352, 48.11657, length = 133))
  pred.df.interp.final<- data.frame(expand.grid(x = pred.df.interp$x, y = pred.df.interp$y), z = c(pred.df.interp$z))
  pred.sp<- st_as_sf(pred.df.interp.final, coords = c("x", "y"), crs = proj.wgs84)
  
  # Baseline -- SDM first and then SDM and NEVA
  data.use<- base.map
  pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
  pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
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
  
  # Discrete scale
  pred.df.use$breaks<- cut(pred.df.use$z, 
                     breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                     labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
  pred.df.use<- pred.df.use %>%
    drop_na(., breaks)
  
  plot.base.sdm<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
    scale_fill_viridis(option = "viridis", name = "SDM Base", na.value = "white", discrete = T, drop = FALSE, direction = 1) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5, 0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) +
    guides(fill=guide_legend(ncol=2))
  
  # Combo
  data.use<- combo.map.base
  pred.df.base<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
  pred.df.interp<- interp(pred.df.base[,1], pred.df.base[,2], pred.df.base[,3], duplicate = "mean", extrap = TRUE,
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
  
  # Discrete scale
  pred.df.use$breaks<- cut(pred.df.use$z, 
                           breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                           labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
  pred.df.use<- pred.df.use %>%
    drop_na(., breaks)
  
  plot.base.combo<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
    scale_fill_viridis(option = "viridis", name = "SDM + NEVA Base", na.value = "white", discrete = T, drop = FALSE) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) +
    guides(fill=guide_legend(ncol=2))
  
  ## Future
  data.use<- fut.map
  pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
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
  
  # Discrete scale
  pred.df.use$breaks<- cut(pred.df.use$z, 
                           breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                           labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
  pred.df.use<- pred.df.use %>%
    drop_na(., breaks)
  
  plot.fut.sdm<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
    scale_fill_viridis(option = "viridis", name = "SDM Fut", na.value = "white", discrete = T, drop = F) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) +
    guides(fill=guide_legend(ncol=2))
  
  # Combo future
  data.use<- combo.map.fut
  pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
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
  
  # Discrete scale
  pred.df.use$breaks<- cut(pred.df.use$z, 
                           breaks = seq(from = 0.0, to = 1.0, by = 0.1), 
                           labels = c("0.0:0.1", "0.1:0.2", "0.2:0.3", "0.3:0.4", "0.4:0.5", "0.5:0.6", "0.6:0.7", "0.7:0.8", "0.8:0.9", "0.9:1.0"), include.lowest = T)
  pred.df.use<- pred.df.use %>%
    drop_na(., breaks)
  
  plot.fut.combo<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = breaks), show.legend = TRUE) +
    scale_fill_viridis(option = "viridis", name = "SDM+NEVA Fut", na.value = "white", discrete = T, drop = F) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg P(Pres) = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum P(Pres) = ", round(sum(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) +
    guides(fill=guide_legend(ncol=2))
  
  # Differences
  # Getting limits
  sdm<- range(data.frame("x" = fut.map$x, "y" = fut.map$y, "pred" = fut.map$pred - base.map$pred)$pred, na.rm = T)
  combo<- range(data.frame("x" = combo.map.fut$x, "y" = combo.map.fut$y, "pred" = combo.map.fut$pred - combo.map.base$pred)$pred, na.rm = T)
  diff.lim<- c(round(min(sdm, combo), 2), round(max(sdm, combo), 2))
  
  # SDM future difference
  data.use<- data.frame("x" = fut.map$x, "y" = fut.map$y, "pred" = fut.map$pred - base.map$pred)
  pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
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
  
  plot.diff.sdm<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
    scale_fill_gradient2(name = "SDM Diff", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg Diff Pres = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum Abs Diff Pres = ", round(sum(abs(pred.df.use$z), na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) 
  
  # Combo
  data.use<- data.frame("x" = combo.map.fut$x, "y" = combo.map.fut$y, "pred" = combo.map.fut$pred - combo.map.base$pred)
  pred.df<- na.omit(data.frame("x" = data.use$x, "y" = data.use$y, "layer" = data.use$pred))
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
  
  plot.diff.combo<- ggplot() + 
    geom_tile(data = pred.df.use, aes(x = X, y = Y, fill = z), show.legend = TRUE) +
    scale_fill_gradient2(name = "SDM+NEVA Diff", low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", limits = diff.lim) +
    geom_map(data = us.states.f, map = us.states.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_map(data = ca.provinces.f, map = ca.provinces.f,
             aes(map_id = id, group = group),
             fill = "gray65", color = "gray45", size = 0.15) +
    geom_label(aes(x = -73, y = 45), label = paste("Avg Diff Pres = ", round(mean(pred.df.use$z, na.rm = T), 2), sep = " ")) + 
    geom_label(aes(x = -73, y = 43), label = paste("Sum Abs Diff Pres = ", round(sum(abs(pred.df.use$z), na.rm = T), 2), sep = " ")) + 
    ylim(ylim.use) + ylab("Lat") +
    scale_x_continuous("Long", breaks = c(-75.0, -70.0, -65.0), labels = c("-75.0", "-70.0", "-65.0"), limits = xlim.use) +
    coord_fixed(1.3) + 
    theme(panel.background = element_rect(fill = "white", color = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.background = element_rect(fill="white", color = "black"), legend.position = c(0.5,0.25), legend.text=element_text(size=10), legend.title=element_text(size=10)) 
  
  out<- plot_grid(plot.base.sdm, plot.fut.sdm, plot.diff.sdm, plot.base.combo, plot.fut.combo, plot.diff.combo, nrow = 2, ncol = 3, align = "hv", scale = 1)
  return(out)
}
