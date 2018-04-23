##### Beta distribution
###########
library(mgcv)
library(tidyverse)
library(cowplot)

# Need a function that allows us to parameterize the beta distribution, which has two positive shape parameters: alpha and beta. 
est_beta_params_func<- function(mu, var) {
  alpha<- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta<- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

bmu<- 0.1
bsd<- 0.2
fmu<- 0.8
fsd<- 0.2

base.params<- est_beta_params_func(mu = bmu, var = bsd^2)
fut.params<- est_beta_params_func(mu = fmu, var = fsd^2)

# What do those two curves look like?
x.vec<- seq(from = 0, to = 1, length.out = 1000)
alpha.base<-base.params$alpha
beta.base<- base.params$beta
alpha.fut<- fut.params$alpha
beta.fut<- fut.params$beta
beta.plot<- ggplot(data.frame(x = x.vec), aes(x = x)) +
  stat_function(fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "blue") +
  stat_function(fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "red")

# How can we describe these two curves --- interquartile ranges?
base.vec<- dbeta(x.vec, shape1 = alpha.base, shape2 = beta.base)
fut.vec<- dbeta(x.vec, shape1 = alpha.fut, shape2 = beta.fut)

dir.iqr<- as.numeric(quantile(base.vec, prob = c(0.25, 0.75)))
base.rect <- data.frame(
  "dir" = c("Negative", "Neutral", "Positive"),
  "from" = c(0.0, dir.iqr[1], dir.iqr[2]), 
  "to" = c(dir.iqr[1], dir.iqr[2], 1.0))

# Area of future curve within those bins?
dtmp<- pbeta(dir.iqr, alpha.fut, beta.fut) 
pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2])

beta.plot.dir<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 3, fill = dir), alpha = 0.4) +
  scale_fill_manual(name = "Directional Effect", values = c("#0571b0", "#f7f7f7", "#ca0020")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.1, label = paste("Negative", round(pd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Neutral", round(pd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Positive", round(pd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence") + 
  ylab("Beta density") +
  ylim(c(-1.5, 3)) +
  theme_bw()
#beta.plot.dir

# Vulnerability
base.vuln<- quantile(base.vec, prob = seq(from = 0, to = 1, length.out = 6))
md<- length(base.vuln)/2+1

base.rect <- data.frame(
  "vuln" = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"),
  "from" = c(base.vuln[3], base.vuln[2], base.vuln[4], base.vuln[1], base.vuln[5], 0,  max(c(base.vuln[6], 1))), 
  "to" = c(base.vuln[4], base.vuln[3], base.vuln[5], base.vuln[2],  max(c(base.vuln[6], 1)), base.vuln[1], 1.0))
base.rect$vuln<- factor(base.rect$vuln, levels = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"))

vtmp<- pbeta(base.vuln, alpha.fut, beta.fut) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])

# Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
vd<- rep(NA, md)
vd[1]<- vtmp[md]

for(k in 1:(md-1)){
  vd[k+1]<- vtmp[md+k]+vtmp[md-k]
}

beta.plot.vuln<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 3, fill = vuln), alpha = 0.4) +
  scale_fill_manual(name = "Vulnerability", values = c("#fef0d9", "#fdcc8a", "#fdcc8a", "#fc8d59", "#fc8d59", "#d7301f", "#d7301f")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.1, label = paste("Low", round(vd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Mod", round(vd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("High", round(vd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("VeryHigh", round(vd[1], 2)), y = -1.25), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence") + 
  ylab("Beta density") +
  ylim(c(-1.5, 3)) +
  theme_bw()
#beta.plot.vuln

# Show em both
temp.out<- plot_grid(beta.plot.dir, beta.plot.vuln, align = "hv")
title<- ggdraw() + draw_label(paste("BaseMn = ", bmu, ",", "BaseSD = ", bsd, ",", "FutMn = ", fmu, ",", "FSD = ", fsd), fontface='bold')
out<- plot_grid(title, temp.out, ncol=1, rel_heights=c(0.1, 1))
ggsave("~/Desktop/out.theoretical.pdf", out, width = 11, height = 8)
dev.off()


#####
## What do our actual distributions look like?
mod<- readRDS("~/Desktop/NormalVoting_BetaDist_IncKnots/gamfitamerican lobster_spring.rds")
season<- "SPRING"

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

# Make predictions
base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
base.pred<- na.omit(as.numeric(predict.gam(mod, newdata = base.preds.use, type = "response")))
fut.pred<- na.omit(as.numeric(predict.gam(mod, newdata = fut.preds.use, type = "response")))

# Get beta parameters
bmn<- mean(base.pred, na.rm = T)
bsd<- sd(base.pred, na.rm = T)
fmn<- mean(fut.pred, na.rm = T)
fsd<- sd(fut.pred, na.rm = T)

base.params<- est_beta_params_func(mu = bmn , var = bsd^2)
fut.params<- est_beta_params_func(mu = fmn, var = fsd^2)

# What do those two curves look like?
x.vec<- seq(from = 0, to = 1, length.out = 1000)
alpha.base<-base.params$alpha
beta.base<- base.params$beta
alpha.fut<- fut.params$alpha
beta.fut<- fut.params$beta
beta.plot<- ggplot(data.frame(x = x.vec), aes(x = x)) +
  stat_function(fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "blue") +
  stat_function(fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "red")

# How can we describe these two curves --- interquartile ranges?
base.vec<- dbeta(x.vec, shape1 = alpha.base, shape2 = beta.base)
fut.vec<- dbeta(x.vec, shape1 = alpha.fut, shape2 = beta.fut)

temp<- base.vec - fut.vec

dir.iqr<- as.numeric(quantile(base.vec, prob = c(0.25, 0.75)))
base.rect <- data.frame(
  "dir" = c("Negative", "Neutral", "Positive"),
  "from" = c(0.0, dir.iqr[1], dir.iqr[2]), 
  "to" = c(dir.iqr[1], dir.iqr[2], 1.0))

# Area of future curve within those bins?
dtmp<- pbeta(dir.iqr, alpha.fut, beta.fut) 
pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2])

beta.plot.dir<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 5, fill = dir), alpha = 0.4) +
  scale_fill_manual(name = "Directional Effect", values = c("#0571b0", "#f7f7f7", "#ca0020")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.1, label = paste("Negative", round(pd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Neutral", round(pd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Positive", round(pd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence") + 
  ylab("Beta density") +
  ylim(c(-1.5, 5)) +
  theme_bw()
#beta.plot.dir

# Vulnerability
base.vuln<- quantile(base.vec, prob = seq(from = 0, to = 1, length.out = 6))
md<- length(base.vuln)/2+1

base.rect <- data.frame(
  "vuln" = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"),
  "from" = c(base.vuln[3], base.vuln[2], base.vuln[4], base.vuln[1], base.vuln[5], 0,  max(c(base.vuln[6], 1))), 
  "to" = c(base.vuln[4], base.vuln[3], base.vuln[5], base.vuln[2],  max(c(base.vuln[6], 1)), base.vuln[1], 1.0))
base.rect$vuln<- factor(base.rect$vuln, levels = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"))

vtmp<- pbeta(base.vuln, alpha.fut, beta.fut) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])

# Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
vd<- rep(NA, md)
vd[1]<- vtmp[md]

for(k in 1:(md-1)){
  vd[k+1]<- vtmp[md+k]+vtmp[md-k]
}

beta.plot.vuln<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 5, fill = vuln), alpha = 0.4) +
  scale_fill_manual(name = "Vulnerability", values = c("#fef0d9", "#fdcc8a", "#fdcc8a", "#fc8d59", "#fc8d59", "#d7301f", "#d7301f")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.base, shape2 = beta.base), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dbeta, args = list(shape1 = alpha.fut, shape2 = beta.fut), color = "black", linetype = "dashed") +
  geom_text(aes(x = 0.1, label = paste("Low", round(vd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Mod", round(vd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("High", round(vd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("VeryHigh", round(vd[4], 2)), y = -1.25), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence") + 
  ylab("Beta density") +
  ylim(c(-1.5, 5)) +
  theme_bw()
#beta.plot.vuln

# Show em both
temp.out<- plot_grid(beta.plot.dir, beta.plot.vuln, align = "hv")
title<- ggdraw() + draw_label(paste("BaseMn = ", round(bmn, 2), ",", "BaseSD = ", round(bsd, 2), ",", "FutMn = ", round(fmn, 2), ",", "FutSD = ", round(fsd, 2)), fontface='bold')
out.real<- plot_grid(title, temp.out, ncol=1, rel_heights=c(0.1, 1))
ggsave("~/Desktop/out.real.pdf", out.real, width = 11, height = 8)
dev.off()

#
plot_grid(out, out.real, nrow = T)


#####
## Normal
#####
## What do our actual distributions look like?
mod<- readRDS("~/Desktop/NormalVoting_BetaDist_IncKnots/gamfitamerican lobster_spring.rds")
season<- "SPRING"

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

# Make predictions
base.preds.use<- base.preds$Data[[match(season, base.preds$SEASON)]]
fut.preds.use<- fut.preds$Data[[match(season, fut.preds$SEASON)]]
base.pred<- as.numeric(predict.gam(mod, newdata = base.preds.use, type = "link"))
base.pred.mu<- mean(base.pred, na.rm = T)
base.pred.sd<- sd(base.pred, na.rm = T)
fut.pred<- as.numeric(predict.gam(mod, newdata = fut.preds.use, type = "link"))
fut.pred.mu<- mean(fut.pred, na.rm = T)
fut.pred.sd<- sd(fut.pred, na.rm = T)
norm.diff<- fut.pred-base.pred

# Get mean and sd of difference
diff.mn<- mean(norm.diff, na.rm = T)
diff.sd<- sd(norm.diff, na.rm = T)

# Normal pdf
# What do those two curves look like?
x.vec<- seq(from = -7, to = 7, length.out = 1000)
norm.plot<- ggplot(data.frame(x = x.vec), aes(x = x.vec)) +
  stat_function(fun = dnorm, args = list(mean = diff.mn, sd = diff.sd), color = "green") +
  stat_function(fun = dnorm, args = list(mean = base.pred.mu, sd = base.pred.sd), color = "blue") +
  stat_function(fun = dnorm, args = list(mean = fut.pred.mu, sd = fut.pred.sd), color = "red")

# IQR of difference??
diff.mu<- mean(norm.diff, na.rm = T)
diff.sd<- sd(norm.diff, na.rm = T)

dir.iqr<- as.numeric(quantile(norm.diff, prob = c(0.25, 0.75), na.rm = T))
base.rect <- data.frame(
  "dir" = c("Negative", "Neutral", "Positive"),
  "from" = c(-4, dir.iqr[1], dir.iqr[2]), 
  "to" = c(dir.iqr[1], dir.iqr[2], 4))

# Area of future curve within those bins?
dtmp<- pnorm(dir.iqr, diff.mu, diff.sd) 
pd<- c(dtmp[1], dtmp[2]-dtmp[1], 1-dtmp[2])

norm.plot.dir<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 5, fill = dir), alpha = 0.4) +
  scale_fill_manual(name = "Directional Effect", values = c("#0571b0", "#f7f7f7", "#ca0020")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(mean = base.pred.mu, sd = base.pred.sd), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(mean = fut.pred.mu, sd = fut.pred.sd), color = "black", linetype = "dashed") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(mean = diff.mu, sd = diff.sd), color = "#2ca25f", linetype = "solid") +
  geom_text(aes(x = -3, label = paste("Negative", round(pd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = -3, label = paste("Neutral", round(pd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = -3, label = paste("Positive", round(pd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence - Logit") + 
  ylab("Norm density") +
  ylim(c(-1.5, 5)) +
  theme_bw()
#beta.plot.dir

# Vulnerability
diff.vuln<- quantile(norm.diff, prob = seq(from = 0, to = 1, length.out = 6), na.rm = T)
md<- length(base.vuln)/2+1

base.rect <- data.frame(
  "vuln" = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"),
  "from" = c(base.vuln[3], base.vuln[2], base.vuln[4], base.vuln[1], base.vuln[5], -10,  base.vuln[6]), 
  "to" = c(base.vuln[4], base.vuln[3], base.vuln[5], base.vuln[2], base.vuln[6], base.vuln[1], 10))
base.rect$vuln<- factor(base.rect$vuln, levels = c("Low", "Mod.Left", "Mod.Right", "High.Left", "High.Right", "VeryHigh.Left", "VeryHigh.Right"))

vtmp<- pnorm(base.vuln, diff.mu, diff.sd) # Cumulative distribution function, p(x) <= each value supplied in Vbrks
vtmp<- c(vtmp[1], diff(vtmp), 1-vtmp[length(vtmp)])

# Storage vector --  allows for more than 4 vulnerability categories, returns probability of being in each bin
vd<- rep(NA, md)
vd[1]<- vtmp[md]

for(k in 1:(md-1)){
  vd[k+1]<- vtmp[md+k]+vtmp[md-k]
}

norm.plot.vuln<- ggplot() + 
  geom_rect(data = base.rect, aes(xmin = from, xmax = to, ymin = 0, ymax = 5, fill = vuln), alpha = 0.4) +
  scale_fill_manual(name = "Vulnerability", values = c("#fef0d9", "#fdcc8a", "#fdcc8a", "#fc8d59", "#fc8d59", "#d7301f", "#d7301f")) +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(mean = base.pred.mu, sd = base.pred.sd), color = "black", linetype = "solid") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(fut.pred.mu, fut.pred.sd), color = "black", linetype = "dashed") +
  stat_function(data = (data.frame(x = x.vec)), aes(x = x.vec), fun = dnorm, args = list(mean = diff.mu, sd = diff.sd), color = "#2ca25f", linetype = "solid") +
  geom_text(aes(x = 0.1, label = paste("Low", round(vd[1], 2)), y = -0.5), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("Mod", round(vd[2], 2)), y = -0.75), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("High", round(vd[3], 2)), y = -1), colour = "black", size = 5, hjust = 0) +
  geom_text(aes(x = 0.1, label = paste("VeryHigh", round(vd[4], 2)), y = -1.25), colour = "black", size = 5, hjust = 0) +
  xlab("Probability of Presence - Logit") + 
  ylab("Norm density") +
  ylim(c(-1.5, 5)) +
  theme_bw()


# Show em both
temp.out<- plot_grid(norm.plot.dir, norm.plot.vuln, align = "hv")
title<- ggdraw() + draw_label(paste("BaseMn = ", round(base.pred.mu, 2), ",", "BaseSD = ", round(base.pred.sd, 2), ",", "FutMn = ", round(fut.pred.mu, 2), ",", "FutSD = ", round(fut.pred.sd, 2)), fontface='bold')
out.real<- plot_grid(title, temp.out, ncol=1, rel_heights=c(0.1, 1))
ggsave("~/Desktop/out.real.norm.pdf", out.real, width = 11, height = 8)
dev.off()

#
plot_grid(out, out.real, nrow = T)




