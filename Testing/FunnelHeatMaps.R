### Remove previous variables
rm(list = ls())

#### Libraries, set seed, set cores ----
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)
library(compiler)
library(metafor)
library(ggplot2)
library(reshape2) # For melt function
library(MASS)
library(scales)


#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 1000000

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(20,1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level mean - need good sense of range for SMD
theta = c(-0.76,  -0.12,  0, 0.12, 0.76)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.005, 0.022, 1.676)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 5
Sd.split <- 0.6

# Size of per unit bias increase
Bias.multiple <- c(0, log(0.85)/(-1.81) * 2, log(0.7225)/(-1.81) * 2)


system.time(Normal.Simulation <- readRDS(file = "FunnelOutRDS"))

Normal.Simulation <- data.table(Normal.Simulation)

## label for titles

bias_type <- "Outcome Reporting Bias"

################## Heat maps


# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
# d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
#   coord_flip(xlim = c(0,100)) + geom_hline(yintercept = theta[3]) +
#   #scale_fill_gradient(low="grey", high="black") +
#   geom_smooth(method = "lm", colour = "black", linetype = "dotted") 
# 
# d + geom_point(alpha = 1/100) + geom_smooth(colour = "black", linetype = "dotted") #+ coord_cartesian(xlim = c(0, 20))
# d + geom_density_2d() + geom_point(alpha = 1/100) + coord_cartesian(xlim = c(0, 300)) 
# d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) + coord_cartesian(xlim = c(0, 300)) 

# New maps

# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(x = Study_sd, y = Study_estimate)) + theme_bw()
# d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
#   coord_cartesian(xlim = c(-1.5,1.5), ylim = c(0, 0.5)) + scale_y_reverse() + geom_hline(yintercept = theta[3]) +
#   #scale_fill_gradient(low="grey", high="black") +
#   geom_smooth(method = "lm", colour = "black", linetype = "dotted") 
# 
# asdf <- Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20, .(Study_sd, Study_estimate)]
# xrng = range(c(0, 0.5))
# yrng = range(c(-1.5, 1.5))
# #d1 <- asdf[, .(d1 = kde2d(Study_sd, Study_estimate, lims=c(xrng, yrng), n=200))]
# 
# 
# # ggplot(df,aes(x=x,y=y))+
# #   stat_density2d(aes(alpha=..level..), geom="polygon") +
# #   scale_alpha_continuous(limits=c(0,0.2),breaks=seq(0,0.2,by=0.025))+
# #   geom_point(colour="red",alpha=0.02)+
# #   theme_bw()

# Working version
# ? Could add diagonal lines for 95% etc

# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(y = Study_sd, x = Study_estimate)) 
# d + stat_density_2d(aes(fill=..level.., alpha = ..level..), geom="polygon", contour = TRUE, n = 50, bins = 10, show.legend = FALSE)  +
#   #geom_density_2d() + 
#   #+ theme_bw()
#   #scale_color_continuous(high = "blue") +
#   #scale_alpha_continuous() + 
#   #scale_alpha_continuous(breaks = trans_breaks(identity, identity, n = 15)) + 
#   scale_y_reverse() +
#   geom_vline(xintercept=0) +
#   #geom_smooth(method = "lm", colour = "black", linetype = "dotted")+ 
#   xlab("Mean difference") +
#   ylab("Standard Error") + 
#   #coord_cartesian(xlim = c(-.25,.25), ylim = c(0.055, 0.1)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())
# 
# # Trying with raster
# 
# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(y = Study_sd, x = Study_estimate)) 
# d + stat_density_2d(aes(fill= ..density..), geom="raster")  +
#   #geom_density_2d() + 
#   #+ theme_bw()
#   #scale_color_continuous(high = "blue") +
#   #scale_alpha_continuous() + 
#   #scale_alpha_continuous(breaks = trans_breaks(identity, identity, n = 15)) + 
#   scale_y_reverse() +
#   geom_vline(xintercept=0) +
#   #geom_smooth(method = "lm", colour = "black", linetype = "dotted")+ 
#   xlab("Mean difference") +
#   ylab("Standard Error") + 
#   #coord_cartesian(xlim = c(-.25,.25), ylim = c(0.055, 0.1)) +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank(),
#         axis.text=element_text(size=12),
#         axis.title=element_text(size=14,face="bold")
#         )

plot1 <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(y = Study_sd, x = Study_estimate))+
stat_density_2d(aes(alpha=..level.., fill = ..level..), geom="polygon", contour = TRUE, n = 200, bins = 10, show.legend = FALSE)  +
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  #geom_smooth(method = "lm", colour = "black", linetype = "dotted")+ 
  xlab("Mean difference") +
  ylab("Standard Error") + 
  ggtitle(paste(bias_type)) + 
  #coord_cartesian(xlim = c(-.25,.25), ylim = c(0.055, 0.1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold")
        )

ggsave(paste("FunnelPlot", bias_type, ".png", sep = ""), plot1, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")

plot1 <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(y = Study_sd, x = Study_estimate))+
  stat_density_2d(aes(alpha=..level.., fill = ..level..), geom="polygon", contour = TRUE, n = 200, bins = 10, show.legend = FALSE)  +
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  #geom_smooth(method = "lm", colour = "black", linetype = "dotted")+ 
  xlab("Mean difference") +
  ylab("Standard Error") + 
  coord_cartesian(xlim = c(-.25,.25), ylim = c(0.055, 0.1)) +
  ggtitle(paste(bias_type)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold")
        )

ggsave(paste("FunnelPlot", bias_type, "FixedSize", ".png", sep = ""), plot1, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")

## Crazy idea - scale density by level of y to create proportion of values at a given x in each tile
asdf <- Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20, .(Study_sd, Study_estimate)]
xrng = range(c(min(Normal.Simulation$Study_sd), 0.5))
yrng = range(c(-.5, .5))

d1 <- kde2d(asdf$Study_sd, asdf$Study_estimate, lims=c(xrng, yrng), n=200)
#contour(d1)

rownames(d1$z) = d1$x
colnames(d1$z) = d1$y

d1.m = melt(d1$z, id.var=rownames(d1))
#names(diff12.m) = c("Duration","Waiting","z")

d1.m <- data.table(d1.m)
d1.m[, sum := sum(value), by = .(Var1)]
#d1.m[sum == 0, value := NA] # Prevents weird graphing later
d1.m[, Prop := value/sum, by = .(Var1, Var2)]


### Attempting to get ggplot2 to do density work using factors

brks <- pretty(range(asdf$Study_sd), n = 50, min.n = 1)
test1 <- cut(asdf$Study_sd, breaks = brks)

d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(y = factor(test1), x = Study_estimate)) 
d + stat_density(aes(fill = ..count..), geom = "raster", position = "identity")
d + geom_point(alpha = 0.01)

### is it legitimate to re-smooth?

ggplot(na.omit(d1.m), aes(y = Var1, x = Var2, z=Prop, fill=Prop)) +
  geom_raster(aes(fill = Prop), interpolate = TRUE) +
  #geom_contour(colour = "white") + 
  scale_fill_continuous(low = "white", high = "blue") + 
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())



ggplot(na.omit(d1.m), aes(y = Var1, x = Var2, z=Prop, fill=Prop)) +
  geom_raster(aes(fill = Prop)) +
  #geom_contour(colour = "white") + 
  scale_fill_continuous(low = "white", high = "blue") + 
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

plot2 <- ggplot(na.omit(d1.m), aes(y = Var1, x = Var2, z=Prop, fill=Prop)) +
  geom_raster(aes(fill = Prop), show.legend = FALSE) +
  geom_contour(colour = "white") + 
  scale_fill_continuous(low = "white", high = "blue") + 
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  xlab("Mean difference") +
  ylab("Standard Error") + 
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave(paste("ProportionFunnelPlot", bias_type, ".png", sep = ""), plot2, dpi = 300, device = "png")


## What happens if you square se?

asdf <- Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20, .(Study_sd, Study_estimate)]
asdf <- asdf[, Study_sd := Study_sd^2]
xrng = range(c(min(Normal.Simulation$Study_sd), 0.25))
yrng = range(c(-.5, .5))

d1 <- kde2d(asdf$Study_sd, asdf$Study_estimate, lims=c(xrng, yrng), n=200)
#contour(d1)

rownames(d1$z) = d1$x
colnames(d1$z) = d1$y

d1.m = melt(d1$z, id.var=rownames(d1))
#names(diff12.m) = c("Duration","Waiting","z")

d1.m <- data.table(d1.m)
d1.m[, sum := sum(value), by = .(Var1)]
#d1.m[sum == 0, value := NA] # Prevents weird graphing later
d1.m[, Prop := value/sum, by = .(Var1, Var2)]

ggplot(na.omit(d1.m), aes(y = Var1, x = Var2, z=Prop, fill=Prop)) +
  geom_raster(aes(fill = Prop)) +
  geom_contour(colour = "white") + 
  scale_fill_continuous(low = "white", high = "blue") + 
  scale_y_reverse() +
  geom_vline(xintercept=0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

plot3 <- ggplot(na.omit(d1.m), aes(y = Var1, x = Var2, z=Prop, fill=Prop)) +
  geom_raster(aes(fill = Prop), show.legend = FALSE) +
  geom_contour(colour = "white") + 
  scale_fill_continuous(low = "white", high = "blue") + 
  scale_y_reverse() +
  xlab("Mean difference") +
  ylab("Standard Error Squared") + 
  geom_vline(xintercept=0) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

ggsave("SquaredProportionFunnelPlotNone.png", plot3, dpi = 300, device = "png")



# d1 <- kde2d(asdf$Study_sd, asdf$Study_estimate, lims=c(xrng, yrng), n=200)
# 
# image(d1)
# 
# d2 <- reshape2::melt(d1$z, id.var = (d1$x))
# names(d2) = c("SE","Estimate","z")
# 
# 
# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(x = Study_sd, y = Study_estimate)) + theme_bw()
# d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE, n = 100)  +
#   #scale_fill_gradient(low="grey", high="black") +
#   #geom_smooth(method = "lm", colour = "black", linetype = "dotted")+ 
#   coord_cartesian(xlim = c(-1.5,1.5), ylim = c(0, 0.5)) + 
#   guides(colour=FALSE)  + geom_vline(xintercept=0) + scale_y_reverse() +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())
# 
# d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 20], aes(x = Study_sd, y = Study_estimate)) + theme_bw()
# d + stat_density_2d(geom = "raster", aes(fill = ..density..), contour = FALSE, n = 200) +
#   scale_fill_continuous(low = "white", high = "blue") + stat_density_2d()
# 
# heat.map <- ggplot(d2 , aes(Estimate, SE, z=z, fill=z)) +
#   geom_tile() +
#   stat_contour(aes(colour=..level..), binwidth=0.1) +
#   scale_fill_gradient2()+ #low="red",mid="white", high="blue", midpoint=0) +
#   scale_colour_gradient2()+ #low=muted("red"), mid="white", high=muted("blue"), midpoint=0)  + 
#   coord_cartesian(xlim = c(-1.5,1.5), ylim = c(0, 0.5)) + 
#   guides(colour=FALSE)  + geom_vline(xintercept=0) + scale_y_reverse() +
#   theme_bw() +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         panel.background = element_blank())
# 
# heat.map
# 
# heat.map2 <- ggplot(d2 , aes(x = Estimate, y = SE, z=z, fill=z)) +
#   stat_contour(aes(colour=..level..), binwidth=0.1)
# 
# heat.map2


#### Import data here ----


system.time(Normal.Simulation2 <- readRDS(file = "FunnelBeggRDS"))

Normal.Simulation2 <- data.table(Normal.Simulation2)



# y = Study_estimate, x = Study_sd^(-2)



den1 <- Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2]
den2 <- Normal.Simulation2[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2]

rm(Normal.Simulation2)
rm(Normal.Simulation)

gc()

# Calculate the common x and y range for geyser1 and geyser2
# xrng = range(c(den1$Study_sd, den2$Study_sd))
# yrng = range(c(den1$Study_estimate, den2$Study_estimate))
xrng = range(c(0, 0.5))
yrng = range(c(-1.5, 1.5))


# Calculate the 2d density estimate over the common range
d1 <- kde2d(den1$Study_sd, den1$Study_estimate, lims=c(xrng, yrng), n=400)
gc()
d2 <- kde2d(den2$Study_sd, den2$Study_estimate, lims=c(xrng, yrng), n=400)
gc()

# Confirm that the grid points for each density estimate are identical
identical(d1$x, d2$x) # TRUE
identical(d1$y, d2$y) # TRUE

# Calculate the difference between the 2d density estimates
diff12 = d1 
diff12$z = d2$z - d1$z

## Melt data into long format
# First, add row and column names (x and y grid values) to the z-value matrix
rownames(diff12$z) = diff12$x
colnames(diff12$z) = diff12$y

# Now melt it to long format
diff12.m = melt(diff12$z, id.var=rownames(diff12))
names(diff12.m) = c("SE","Estimate","z")

## Plot difference between geyser2 and geyser1 density
# ggplot(diff12.m, aes(SE, Estimate, z=z, fill=z)) +
#   geom_tile() +
#   stat_contour(aes(colour=..level..), binwidth=0.001) +
#   scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
#   scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0)  +
#   coord_flip(xlim=c(1, 0), ylim=c(-2,2)) +
#   guides(colour=FALSE)  + geom_hline(yintercept=0) + scale_x_reverse()

heat.map <- ggplot(diff12.m, aes(Estimate, SE, z=z, fill=z)) +
  geom_tile() +
  stat_contour(aes(colour=..level..), binwidth=0.1) +
  scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
  scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0)  + 
  coord_cartesian(xlim = c(-1.5,1.5), ylim = c(0, 0.5)) + 
  guides(colour=FALSE)  + geom_vline(xintercept=0) + scale_y_reverse() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

heat.map

ggsave("BeggHeatMap.png", heat.map, dpi = 300, device = "png")

#### Copy from online
## Plot difference between geyser2 and geyser1 density
# ggplot(diff12.m, aes(SE, Estimate, z=z, fill=z)) +
#   geom_tile() +
#   stat_contour(aes(colour=..level..), binwidth=0.001) +
#   scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
#   scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0) +
#   coord_cartesian(xlim=c(0,1), ylim=yrng) +
#   guides(colour=FALSE)

object.size(x=lapply(ls(), get))
print(object.size(x=lapply(ls(), get)), units="Mb")

##### Need to fix for precision
