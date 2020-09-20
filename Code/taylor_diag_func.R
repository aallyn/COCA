#####
## Taylor Diagram Function
#####

# Helper functions: correlation coefficient and bias ----------------------
corrcoeff_func_simp<- function(df){
  df.use<- df %>%
    drop_na(obs, mod)
  mean.obs<- mean(df.use$obs)
  mean.mod<- mean(df.use$mod)
  sd.obs<- sd(df.use$obs)
  sd.mod<- sd(df.use$mod)
  samps<- nrow(df.use)
  
  out<- ((1/samps)*(sum((df.use$mod - mean.mod)*(df.use$obs - mean.obs))))/(sd.obs*sd.mod)
  return(out)
}

bias_func_simp<- function(df){
  df.use<- df %>%
    drop_na(obs, mod)
  out<- sd(df.use$mod)/sd(df.use$obs)
  return(out)
}


# Main Taylor Diagram function --------------------------------------------
taylor_diagram_func<- function(dat, obs = "obs", mod = "mod", group = NULL, out.dir, grad.corr.lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1), pcex = 1, cex.axis = 1, normalize = TRUE, mar = c(5, 4, 6, 6), sd.r = 1, pt.col = NULL, pt.cols = NULL, example = FALSE) {
  
  ## Details
  # This function plots a Taylor Diagram of model prediction accuracy, sumamrizing the root mean square error, the coefficient of determination, and the ratio of standard deviations. 
  
  # Args:
  # dat = data frame with observations and model predictions, as well as group if necessary
  # obs = Column name for the observation response
  # mod = Column name for the modeled response
  # group = Grouping variable, used for comparing different species/ages or stages/models, etc
  # out.dir = Directory where to save the Taylor Diagram plot
  # ... All these other things correspond to some of the aesthetics of the plot. pt.col gives color if just plotting one point (group, model), pt.cols is a vector of colors for plotting multiple points (groups, models) on one plot.

  # Returns: NULL; saves plot to output directory
  
  ## Start function
  # Install libraries
  library(tidyverse)
  
  # Set arguments for debugging -- this will NOT run when you call the function. Though, you can run each line inside the {} and then you will have everything you need to walk through the rest of the function.
  if(example){
    # Create a data set with observations and predictions
    data(trees)
    tree.mod<- lm(Volume ~ Girth, data = trees)
    trees$Volume.pred<- as.numeric(predict(tree.mod))
    dat<- trees
    obs<- "Volume"
    mod<- "Volume.pred"
    group<- NULL
    out.dir<- "~/Desktop/"
    grad.corr.lines = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1)
    pcex = 1
    cex.axis = 1
    normalize = TRUE
    mar = c(5, 4, 6, 6)
    sd.r = 1
    pt.col<- "#006d2c"
  }
  
  # Some house keeping -- rename the obs and mod columns to work with generic functions
  old.names<- c(obs, mod)
  new.names<- c("obs", "mod")
  dat<- dat %>%
    rename_at(vars(old.names), ~ new.names)
  
  # Calculate the correlation coefficient and bias, two stats needed for Taylor Diagram. Flexibility to group and then calculate stats by group (e.g., species, model, etc).
  if(is.null(group)){
    mod.stats<- dat %>%
      nest(data = everything()) %>%
      mutate(., "CorrCoeff" = as.numeric(map(data, corrcoeff_func_simp)),
             "Bias" = as.numeric(map(data, bias_func_simp)))
  } else {
    # Group by group and calculate stats
    mod.stats<- dat %>%
      group_by(group) %>%
      nest() %>%
      mutate(., "CorrCoeff" = as.numeric(map(data, corrcoeff_func_simp)),
             "Bias" = as.numeric(map(data, bias_func_simp)))
  }

  # Now plot creation....
  # Getting maxSD for plotting
  maxsd<- max(mod.stats$Bias, 1)
  
  # Empty plot first
  # Creating empty plot first
  plot.base<- ggplot() + 
    scale_x_continuous(name = "Standard deviation (normalized)", limits = c(0, maxsd+0.01), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0.01)) +
    scale_y_continuous(name = "Standard deviation (normalized)", limits = c(-0.015, maxsd), breaks = seq(from = 0, to = maxsd, by = 0.5), expand = c(0, 0.01)) +
    theme_bw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16, face="bold"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) 
  
  # Coeff D rays 
  for(i in 1:length(grad.corr.lines)){
    x.vec<- c(0, maxsd*grad.corr.lines[i])
    y.vec<- c(0, maxsd*sqrt(1 - grad.corr.lines[i]^2))
    
    if(i ==1){
      coeffd.rays.df<- data.frame("Ray" = rep(1, length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Ray" = rep(i, length(x.vec)), "x" = x.vec, "y" = y.vec)
      coeffd.rays.df<- bind_rows(coeffd.rays.df, temp)
    }
  }
  
  # Add rays
  plot.coeffd<- plot.base +
    geom_line(data = coeffd.rays.df, aes(x = x, y = y, group = Ray), lty = "longdash", col = "lightgray") 
    
  coeffd.labs<- coeffd.rays.df %>%
    group_by(Ray) %>%
    summarize(., 
              "x" = max(x, na.rm = TRUE), 
              "y" = max(y, na.rm = TRUE)) %>%
    data.frame()
  
  coeffd.labs$Label<- grad.corr.lines
  
  plot.coeffd<- plot.coeffd +
    geom_label(data = coeffd.labs, aes(x = x, y = y, label = Label), fill = "white", label.size = NA, size = 5)
  
  # SD arcs
  # Need to add in SD arcs
  sd.arcs<- seq(from = 0, to = maxsd, by = 0.5)
  
  for(i in 1:length(sd.arcs)){
    x.vec<- sd.arcs[i]*cos(seq(0, pi/2, by = 0.03))
    y.vec<- sd.arcs[i]*sin(seq(0, pi/2, by = 0.03))
    
    if(i ==1){
      sd.arcs.df<- data.frame("Arc" = rep(sd.arcs[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Arc" = rep(sd.arcs[i], length(x.vec)), "x" = x.vec, "y" = y.vec)
      sd.arcs.df<- bind_rows(sd.arcs.df, temp)
    }
  }
  
  # Add arcs to plot.base
  plot.sd<- plot.coeffd +
    geom_line(data = sd.arcs.df, aes(x = x, y = y, group = Arc), lty = "dotted", color = "lightgray") 
  
  # Now gamma? -- Standard deviation arcs around the reference point
  gamma<- pretty(c(0, maxsd), n = 4)[-1]
  gamma<- gamma[-length(gamma)]
  labelpos<- seq(45, 70, length.out = length(gamma))
  
  for(gindex in 1:length(gamma)) {
    xcurve <- cos(seq(0, pi, by = 0.03)) * gamma[gindex] + sd.r
    endcurve <- which(xcurve < 0)
    endcurve <- ifelse(length(endcurve), min(endcurve) - 1, 105)
    ycurve <- sin(seq(0, pi, by = 0.03)) * gamma[gindex]
    maxcurve <- xcurve * xcurve + ycurve * ycurve
    startcurve <- which(maxcurve > maxsd * maxsd)
    startcurve <- ifelse(length(startcurve), max(startcurve) + 1, 0)
    x.vec<- xcurve[startcurve:endcurve]
    y.vec<- ycurve[startcurve:endcurve]
    
    if(gindex ==1){
      gamma.df<- data.frame("Gamma" = rep(gamma[1], length(x.vec)), "x" = x.vec, "y" = y.vec)
    } else {
      temp<- data.frame("Gamma" = rep(gamma[gindex], length(x.vec)), "x" = x.vec, "y" = y.vec)
      gamma.df<- bind_rows(gamma.df, temp)
    }
  }
  
  gamma.df$Gamma<- factor(gamma.df$Gamma, levels = unique(gamma.df$Gamma))
  
  # Add em
  plot.gamma<- plot.sd +
    geom_line(data = gamma.df, aes(x = x, y = y, group = Gamma), lty = "solid", col = "lightgray")
  
  # Label...
  gamma.labs<- gamma.df %>%
    group_by(Gamma) %>%
    summarize("x" = mean(x, na.rm = TRUE), 
              "y" = median(y, na.rm = TRUE))
  
  inflection_func<- function(df){
    d1<- diff(df$y)/diff(df$x)    
    pt.id<- which.max(d1)
    pt.out<- df[pt.id,]
    pt.out$y<- rep(0, nrow(pt.out))
    return(pt.out)
  }
  
  gamma.labs<- gamma.df %>%
    group_by(Gamma) %>%
    nest() %>%
    summarize("pt" = map(data, inflection_func)) %>%
    unnest(cols = c(pt))
  
  plot.gamma<- plot.gamma +
    geom_label(data = gamma.labs, aes(x = x, y = y, label = Gamma), fill = "white", label.size = NA)
  
  # Add in reference point
  plot.all<- plot.gamma +
    geom_point(aes(x = sd.r, y = 0), color = "black", size = 3.5) 
  
  # Add in reference points
  mod.td<- mod.stats %>%
    mutate(., "TD.X" = Bias * CorrCoeff,
           "TD.Y" = Bias * sin(acos(CorrCoeff)))
  
  if(is.null(group)){
    plot.td<- plot.all +
      geom_point(data = mod.td, aes(x = TD.X, y = TD.Y), color = pt.col, size = 4) +
      geom_text(aes(label = "Correlation coefficient", x = 0.8, y = 0.75), fontface = "bold", angle = -38, size = 5.25)
  } else {
    plot.td<- plot.all +
      geom_point(data = mod.td, aes(x = TD.X, y = TD.Y, color = pt.cols), size = 3) +
      geom_text(aes(label = "Correlation coefficient", x = 0.8, y = 0.75), fontface = "bold", angle = -38, size = 5.25)
  }
  ggsave(paste(out.dir, "_TaylorDiagram.jpg", sep = ""), plot.td)
}

