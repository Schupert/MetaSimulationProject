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
library(stargazer)
library(reshape2)

#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 1000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.8), log(1), log(1.25), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw)
tau.sq = c(0, 0.01777778, 0.04, 3.04)

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

#### Functions -----

CI.betw <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 1, 0)
  return(output)
}

#### Import data here ----

system.time(LogOR.Sim.Results <- readRDS(file = "LSBOutV1An"))

LogOR.Sim.Results <- data.table(LogOR.Sim.Results)

#### Select subset of data for analysis ----
# Currently moving with number of studies in analysis

sig.level <- (1 - 0.05/2)

#An.Cond <- LogOR.Sim.Results[Rep_theta == 0 & Rep_tau.sq == 0 & Rep_Subj == 4.2 & Rep_NumStudies == 100]
An.Cond <- LogOR.Sim.Results[Rep_theta == log(1) & Rep_tau.sq == 0 & Rep_Subj == 4.70 & Rep_ev_freq == 0.5]

## There may be cleaner way to do this with data.table
An.Cond$FE_CIlb <- An.Cond$FE_Estimate - qnorm(sig.level) * An.Cond$FE_se
An.Cond$FE_CIub <- An.Cond$FE_Estimate + qnorm(sig.level) * An.Cond$FE_se

An.Cond$REML_CIlb <- An.Cond$REML_Estimate - qnorm(sig.level) * An.Cond$REML_se
An.Cond$REML_CIub <- An.Cond$REML_Estimate + qnorm(sig.level) * An.Cond$REML_se

An.Cond$DL_CIlb <- An.Cond$DL_Estimate - qnorm(sig.level) * An.Cond$DL_se
An.Cond$DL_CIub <- An.Cond$DL_Estimate + qnorm(sig.level) * An.Cond$DL_se

### Does Moreno use z-score?

An.Cond$Moreno_CIlb <- An.Cond$Moreno_Estimate - qnorm(sig.level) * An.Cond$Moreno_se
An.Cond$Moreno_CIub <- An.Cond$Moreno_Estimate + qnorm(sig.level) * An.Cond$Moreno_se

An.Cond$Mult_CIlb <- An.Cond$FE_Estimate - qnorm(sig.level) * An.Cond$Mult_se
An.Cond$Mult_CIub <- An.Cond$FE_Estimate + qnorm(sig.level) * An.Cond$Mult_se

#### Bias ----

An.Cond[, .(Bias = mean(FE_Estimate) - Rep_theta, MSE1 = mean((FE_Estimate - Rep_theta)^2), 
            MSE2 = (mean(FE_Estimate) - Rep_theta) + var(FE_Estimate), 
            Coverage = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub))), by = .(Rep_theta, Rep_NumStudies)]

Bias.values <- An.Cond[, .(FE = mean(FE_Estimate) - Rep_theta, REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                           DL = mean(DL_Estimate) - Rep_theta, Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
                       by = .(Rep_theta, Rep_NumStudies)]


Bias.values <- melt(Bias.values, id = c("Rep_theta", "Rep_NumStudies"))

### Bias plot

bias.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Bias.values) + xlab("Number of Studies") + ylab("Bias")
bias.plot <- ggplot(Bias.values, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + ylab("Bias") + scale_colour_grey()
bias.plot

#### MSE ----

MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                           DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                       by = .(Rep_theta, Rep_NumStudies)]

MSE1.values <- melt(MSE1.values, id = c("Rep_theta", "Rep_NumStudies"))

MSE1.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = MSE1.values) + coord_cartesian(ylim = c(0, 0.1)) + xlab("Number of Studies") + ylab("MSE")
MSE1.plot <- ggplot(MSE1.values, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + ylab("MSE") + scale_colour_grey()
MSE1.plot

MSE2.values <- An.Cond[, .(FE = (mean(FE_Estimate) - Rep_theta) + var(FE_Estimate), REML = (mean(REML_Estimate, na.rm = TRUE) - Rep_theta) + var(REML_Estimate, na.rm = TRUE),
                           DL = (mean(DL_Estimate, na.rm = TRUE) - Rep_theta) + var(DL_Estimate, na.rm = TRUE), Moreno = (mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta) + var(Moreno_Estimate, na.rm = TRUE) ),
                       by = .(Rep_theta, Rep_NumStudies)]

MSE2.values <- melt(MSE2.values, id = c("Rep_theta", "Rep_NumStudies"))

MSE2.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = MSE2.values) #+ coord_cartesian(ylim = c(-0.1, 0.1))
MSE2.plot

#### Coverage ---- 

Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                               DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), HC.DL = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                               HC.REML = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                               KH.DL = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                               KH.REML = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                               Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                               Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies)]

Coverage.values2 <- melt(Coverage.values, id = c("Rep_theta", "Rep_NumStudies"))

Coverage.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) + xlab("Coverage") + ylab("Number of studies")
Coverage.plot <- ggplot(Coverage.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Coverage") + ylab("Number of studies") + scale_colour_grey()
Coverage.plot


#### Testing stargazer -----


stargazer(Coverage.values[, -1, with = FALSE], summary = FALSE, rownames = FALSE)
