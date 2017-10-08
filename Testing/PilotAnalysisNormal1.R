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
True.sd = sqrt(2)

# theta = population level mean - need good sense of range for SMD
theta = c(-1.5, -0.3, 0, 0.3, 1.5)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.007, 0.133, 2.533)

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
Bias.multiple <- c(0, log(0.9)/(-1.81) * 2, log(0.81)/(-1.81) * 2)

system.time(Normal.Simulation <- readRDS(file = "NSB0V1RDS"))
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

identical(1:15840000, Normal.Simulation$Unique_ID)

### Number of NAs

sum(is.na(Normal.Simulation))

### Overall Summary

summary(Normal.Simulation)

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_NumStudies)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_Subj)]

### Why is sd(est) not similar to tau.sq?
Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd), SD.est = var(Study_estimate) - (mean(Study_sd)^2)), by = .(Rep_tau.sq)]

Normal.Simulation[, .(Estimate = mean(Study_estimate), SD = mean(Study_sd)), by = .(Rep_theta)]



################## Slightly confused issue, these are more relevant for actual analysis of results

### Bias

Normal.Simulation[, .(Bias = mean(Study_estimate) - Rep_theta), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]

### MSE - two different ways of calculating

Normal.Simulation[, .(MSE = mean((Study_estimate - Rep_theta)^2)), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]

Normal.Simulation[, .(MSE = (mean(Study_estimate) - Rep_theta) + var(Study_estimate)), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta, Rep_Subj)]


#### Selected for fixed sample size

Normal.Simulation[Rep_Subj == 60, .(Bias = mean(Study_estimate) - Rep_theta), by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta)]

ggplot(Normal.Simulation[Rep_Subj == 60], aes(Rep_tau.sq, (Bias = mean(Study_estimate) - Rep_theta)))
