#### Analysis

### Remove previous variables
rm(list = ls())

#### Libraries, set seed, set cores ----
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)
library(compiler)
library(ggplot2)
enableJIT(3)

set.seed(123)

# Number of cores for parallel
#num.Cores <- detectCores() - 1
num.Cores <- 16
c1 <- makeCluster(num.Cores)

#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 1000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level mean - need good sense of range for SMD
theta = c(-1.53, -0.25, 0, 0.25, 1.53)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.004, 0.067, 1.267)

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

#### Import data here ----


system.time(Normal.Simulation <- readRDS(file = "NSB0V1"))

Normal.Simulation <- data.table(Normal.Simulation)

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
Normal.Simulation <- Normal.Simulation[order(Unique_ID)]

ID =  length(Subj) * length(theta) * length(tau.sq) * Reps * sum(Studies)
Normal.Simulation$Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
Normal.Simulation$Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
Normal.Simulation$Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(tau.sq)))
Normal.Simulation$Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(60, 20, 250, 4.2)
Normal.Simulation$Rep_Subj = rep(Subj2, each = ID / length(Subj))



################ Start checks

### Is Unique ID stable

identical(1:dim(Normal.Simulation)[1], Normal.Simulation$Unique_ID)

### Number of NAs

sum(is.na(Normal.Simulation))

### Overall Summary

summary(Normal.Simulation)

asdf <- Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_NumStudies)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_Subj)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_theta)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_tau.sq)]

### SD(est) is equal to mean(sd) in condition of no heterogeneity
Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sd(Study_estimate) ), by = .(Rep_Subj)]

### This shows that variance is a combination of tau2 and sigma as expected
Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 2.533, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sqrt(var(Study_estimate) - 2.533) ), by = .(Rep_Subj)]

Normal.Simulation[Rep_theta == 1.5 & Rep_tau.sq == 2.533, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sqrt(var(Study_estimate) - 2.533) ), by = .(Rep_Subj)]


#### Bias of estimates of effect
Bias.values <- Normal.Simulation[, .(Bias = mean(Study_estimate) - Rep_theta), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]$Bias
summary(Bias.values)
hist(Bias.values)


#### Plot of curve

### Plots

d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,100)) + geom_hline(yintercept = theta[3]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

d + geom_point(alpha = 1/100) + geom_smooth(colour = "black", linetype = "dotted") #+ coord_cartesian(xlim = c(0, 20))
d + geom_density_2d() + geom_point(alpha = 1/100) + coord_cartesian(xlim = c(0, 300)) 
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE) + coord_cartesian(xlim = c(0, 300)) 

### Increasing tau2

d <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == tau.sq[4] & Rep_Subj == 4.2], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,70)) + geom_hline(yintercept = theta[3]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

## Lowess

d <- ggplot(Normal.Simulation[Rep_theta == theta[3] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 4.2], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,70)) + geom_hline(yintercept = theta[3]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth( colour = "black", linetype = "dotted") 

### Opposing direction

d <- ggplot(Normal.Simulation[Rep_theta == theta[5] & Rep_tau.sq == tau.sq[4] & Rep_Subj == 4.2], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,80)) + geom_hline(yintercept = theta[5]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

## Small opposite effect

d <- ggplot(Normal.Simulation[Rep_theta == theta[4] & Rep_tau.sq == tau.sq[3] & Rep_Subj == 4.2], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,80)) + geom_hline(yintercept = theta[4]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

a <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2], aes(x = Study_estimate, y = (1/(Study_sd^2)) )) 
a + geom_point() + geom_smooth(data = Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2], aes(x = Study_estimate, y = (1/(Study_sd^2))), method = "lm", formula = x~y)

b <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2], aes(y = Study_estimate, x = Study_sd^(-2) )) #+ xlim(0, 100) + ylim(-2,2)
b + geom_point(alpha = 1/10) + geom_smooth(method = "lm")
b + stat_density2d(aes(fill = ..level..), geom = "polygon",  contour = TRUE)  + geom_smooth(method = "lm") + xlim(0, 100) + ylim(-2,2) + coord_flip() + geom_hline(yintercept=0)
b + stat_density2d(aes(fill = ..level..), geom = "polygon",  contour = TRUE)  + geom_smooth(span = 0.8) + xlim(0, 100) + ylim(-2,2) + coord_flip() + geom_hline(yintercept=0)

b + geom_density2d()

## Testing with tau2

b <- ggplot(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 2.533 & Rep_Subj == 4.2], aes(y = Study_estimate, x = Study_sd^(-2) )) #+ xlim(0, 100) + ylim(-2,2)
b + geom_point(alpha = 1/10) + geom_smooth(method = "lm") + coord_flip()
b + stat_density2d(aes(fill = ..level..), geom = "polygon",  contour = TRUE)  + geom_smooth(method = "lm") + xlim(0, 100) + ylim(-2,2) + coord_flip() + geom_hline(yintercept=0)
b + stat_density2d(aes(fill = ..level..), geom = "polygon",  contour = TRUE)  + geom_smooth(span = 0.8) + xlim(0, 100) + ylim(-2,2) + coord_flip() + geom_hline(yintercept=0)


## Testing with effect size difference

b <- ggplot(Normal.Simulation[Rep_theta == -0.3 & Rep_tau.sq == 0.133 & Rep_Subj == 4.2], aes(y = Study_estimate, x = Study_sd^(-2) )) + theme_bw()#+ xlim(0, 100) + ylim(-2,2)
b + geom_point(alpha = 1/10) + geom_smooth(method = "lm") + coord_flip()

b + stat_density2d(aes(fill = ..level..), geom = "raster",  contour = TRUE)  + geom_smooth(method = "lm") + xlim(0, 100) + ylim(-2,2) +
  coord_flip() + geom_hline(yintercept=-0.3) #+ scale_fill_gradient(limits = c(1, 10), low = 'gray80', high = 'gray20')

b + stat_density2d(aes(fill = ..level..), geom = "polygon",  contour = TRUE)  + geom_smooth(span = 0.8) + xlim(0, 100) + ylim(-2,2) + coord_flip() + geom_hline(yintercept=0.3)



summary(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2])
summary(Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2]$Study_sd^(-2) )

#### Get MCE from Outcome estimates

################## Heat maps

#### Import data here ----


system.time(Normal.Simulation2 <- readRDS(file = "NSBOutV1"))

Normal.Simulation2 <- data.table(Normal.Simulation2)

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
Normal.Simulation2 <- Normal.Simulation2[order(Unique_ID)]

ID =  length(Subj) * length(theta) * length(tau.sq) * Reps * sum(Studies)
Normal.Simulation2$Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
Normal.Simulation2$Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
Normal.Simulation2$Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(tau.sq)))
Normal.Simulation2$Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(60, 20, 250, 4.2)
Normal.Simulation2$Rep_Subj = rep(Subj2, each = ID / length(Subj))

# y = Study_estimate, x = Study_sd^(-2)

library(reshape2) # For melt function
library(MASS)
library(scales)

den1 <- Normal.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2]
den2 <- Normal.Simulation2[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2]

rm(Normal.Simulation2)

# Calculate the common x and y range for geyser1 and geyser2
xrng = range(c(den1$Study_sd, den2$Study_sd))
yrng = range(c(den1$Study_estimate, den2$Study_estimate))

# Calculate the 2d density estimate over the common range
d1 = kde2d(den1$Study_sd, den1$Study_estimate, lims=c(xrng, yrng), n=500)
d2 = kde2d(den2$Study_sd, den2$Study_estimate, lims=c(xrng, yrng), n=500)

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

ggplot(diff12.m, aes(Estimate, SE, z=z, fill=z)) +
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

#### Copy from online
## Plot difference between geyser2 and geyser1 density
# ggplot(diff12.m, aes(SE, Estimate, z=z, fill=z)) +
#   geom_tile() +
#   stat_contour(aes(colour=..level..), binwidth=0.001) +
#   scale_fill_gradient2(low="red",mid="white", high="blue", midpoint=0) +
#   scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"), midpoint=0) +
#   coord_cartesian(xlim=c(0,1), ylim=yrng) +
#   guides(colour=FALSE)

