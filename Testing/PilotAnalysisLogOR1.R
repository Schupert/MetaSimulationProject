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
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 1000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.8), log(1), log(1.25), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw)
tau.sq = c(0, 0.008, 0.04, 3.04)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3, 0.5)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias
Tested.outcomes <- 5
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.85

##### Import data here ----

system.time(LogOR.Simulation <- readRDS(file = "LSBOutV1"))

LogOR.Simulation <- data.table(LogOR.Simulation)

LogOR.Simulation[,c("Study_G1O1", "Study_G2O1") := NULL]

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
LogOR.Simulation <- LogOR.Simulation[order(Unique_ID)]

# ID = total number of data points required, also used as an ID number. 
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * sum(Studies)

LogOR.Simulation$Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
LogOR.Simulation$Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
rm(intermediate)
LogOR.Simulation$Rep_ev_freq = rep(rep(EvFreq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(EvFreq)))
LogOR.Simulation$Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)*length(EvFreq)), times = ID/(Reps*sum(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Simulation$Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(100, 20, 250, 4.7)
LogOR.Simulation$Rep_Subj = rep(Subj2, each = ID / length(Subj))



################ Start checks

### Is Unique ID stable

identical(1:dim(LogOR.Simulation)[1], LogOR.Simulation$Unique_ID)

### Number of NAs

sum(is.na(LogOR.Simulation))

### Overall Summary

summary(LogOR.Simulation)

asdf <- LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_tau.sq, Rep_theta, Rep_Subj)]

# LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_NumStudies)]

LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_Subj)]

LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_theta)]

LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_tau.sq)]

LogOR.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_ev_freq)]

### SD(est) is equal to mean(sd) in condition of no heterogeneity
LogOR.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_ev_freq == 0.5, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sd(Study_estimate) ), by = .(Rep_Subj)]

### This shows that variance is a combination of tau2 and sigma as expected
LogOR.Simulation[Rep_theta == 0 & Rep_tau.sq == 3.04 & Rep_ev_freq == 0.5, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sqrt(var(Study_estimate) - 3.04) ), by = .(Rep_Subj)]

LogOR.Simulation[Rep_theta == log(4) & Rep_tau.sq == 3.04 & Rep_ev_freq == 0.5, .(Estimate = mean(Study_estimate), Est.S.E = mean(Study_sd), Var.out.minus.tau2 = sqrt(var(Study_estimate) - 3.04) ), by = .(Rep_Subj)]

### Plots

d <- ggplot(LogOR.Simulation[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.7 & Rep_ev_freq == 0.5], aes(x = Study_sd^(-2), y = Study_estimate)) + theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,70)) + geom_hline(yintercept = theta[3]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

d + geom_density_2d() + geom_point(alpha = 1/100) + coord_cartesian(xlim = c(0,100))

### With increasing tausq

d <- ggplot(LogOR.Simulation[Rep_theta == theta[4] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 4.7 & Rep_ev_freq == 0.5], aes(x = Study_sd^(-2), y = Study_estimate))+ theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,70)) + geom_hline(yintercept = log(0)) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 

### Other effects, event sizes, tausq etc

d <- ggplot(LogOR.Simulation[Rep_theta == theta[2] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 4.7 & Rep_ev_freq == 0.1], aes(x = Study_sd^(-2), y = Study_estimate))+ theme_bw()
d + stat_density_2d(aes(fill = ..level..), geom = "polygon", contour = TRUE)  +
  coord_flip(xlim = c(0,20)) + geom_hline(yintercept = theta[2]) +
  scale_fill_gradient(low="grey", high="black") + geom_smooth(method = "lm", colour = "black", linetype = "dotted") 


d + geom_point(alpha = 1/10)

#### Bias of estimates of effect
Bias.values <- LogOR.Simulation[, .(Bias = mean(Study_estimate) - Rep_theta), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]$Bias
summary(Bias.values)
hist(Bias.values)