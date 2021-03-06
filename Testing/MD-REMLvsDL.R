#### Looped Analysis UMD

### Remove previous variables
#rm(list = ls())
rm(list=setdiff(ls(), "UMD.Complete"))

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

#### Functions -----

CI.betw <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 1, 0)
  return(output)
}

CI.updown <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 0, ifelse(a < b, -1, 1))
  return(output)
}

#### UMD Declare variables ----

# Reps = number of repetitions of experiment
Reps = 10000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level mean - need good sense of range for SMD
theta = c( -0.76,  -0.12,  0, 0.12, 0.76)

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


#### UMD Import data here ----

system.time(Normal.Sim.Results <- readRDS(file = "NSTotalV2RDS"))

Normal.Sim.Results <- data.table(Normal.Sim.Results)

##### Confidence intervals ----

sig.level <- (1 - 0.05/2)

Normal.Sim.Results[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
Normal.Sim.Results[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]

Normal.Sim.Results[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
Normal.Sim.Results[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]

Normal.Sim.Results[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
Normal.Sim.Results[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]

### Does Moreno use z-score?

Normal.Sim.Results[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
Normal.Sim.Results[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]

Normal.Sim.Results[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
Normal.Sim.Results[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]


Normal.Sim.Results[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
Normal.Sim.Results[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]

#### Testing large scale summaries ----

# Bias.values <- Normal.Sim.Results[, .(Bias_FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, Bias_REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
#                            Bias_DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Bias_Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
#                        by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
# 
# MSE.values <- Normal.Sim.Results[, .(MSE_FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), MSE_REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
#                            MSE_DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "MSE Moreno" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
#                        by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
# 
# Coverage.values <- Normal.Sim.Results[, .(Coverage_FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Coverage_REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
#                                Coverage_DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Coverage_HC_DL = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
#                                Coverage_HC_REML = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
#                                Coverage_KH_DL = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
#                                Coverage_KH_REML = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
#                                Coverage_IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
#                                Coverage_Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
#                                Coverage_Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
# ),
# by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
# 
# Outside.values <- Normal.Sim.Results[, .(Outside_FE = mean(CI.updown(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Outside_REML = mean(CI.updown(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
#                               Outside_DL = mean(CI.updown(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Outside_HC_DL = mean(CI.updown(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
#                               Outside_HC_REML = mean(CI.updown(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
#                               Outside_KH_DL = mean(CI.updown(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
#                               Outside_KH_REML = mean(CI.updown(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
#                               Outside_IVHet = mean(CI.updown(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
#                               Outside_Moreno = mean(CI.updown(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
#                               Outside_Mult = mean(CI.updown(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
# ),
# by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
# 
# 
# 
# Intermediate1 <- merge(MSE.values, Bias.values, by = c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))
# Intermediate2 <- merge(Intermediate1, Outside.values, by = c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))
# Total.results <- merge(Intermediate2, Coverage.values, by = c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

Total.results <- Normal.Sim.Results[, .(Bias_FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, Bias_REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                                        Bias_DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Bias_Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta,
                                        MSE_FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), MSE_REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                                        MSE_DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "MSE Moreno" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ,
                                        Coverage_FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Coverage_REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                        Coverage_DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Coverage_HC_DL = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                        Coverage_HC_REML = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                        Coverage_KH_DL = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                        Coverage_KH_REML = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                        Coverage_IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                        Coverage_Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                        Coverage_Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                                        Outside_FE = mean(CI.updown(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Outside_REML = mean(CI.updown(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                        Outside_DL = mean(CI.updown(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Outside_HC_DL = mean(CI.updown(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                        Outside_HC_REML = mean(CI.updown(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                        Outside_KH_DL = mean(CI.updown(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                        Outside_KH_REML = mean(CI.updown(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                        Outside_IVHet = mean(CI.updown(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                        Outside_Moreno = mean(CI.updown(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                        Outside_Mult = mean(CI.updown(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                                        DL_I2 = mean(DL_I2, na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]



#### REML versus DL Bias ----

summary(Total.results$Bias_REML)
summary(Total.results$Bias_DL)

a <- Normal.Sim.Results[, .(RE_diff = REML_Estimate - DL_Estimate)]
summary(a$RE_diff)
hist(a$RE_diff)

Normal.Sim.Results[which.max(a$RE_diff)]

b <- Normal.Sim.Results[, .(RE_diff = mean(REML_Estimate - DL_Estimate, na.rm = TRUE),
                           Better = ifelse(abs(mean(REML_Estimate, na.rm = TRUE)) < abs(mean(DL_Estimate, na.rm = TRUE)), "REML", "DL")),
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
summary(b$RE_diff)
summary(as.factor(b$Better))

plot(b$Rep_theta, b$RE_diff)
plot(b$Rep_Subj, b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point(alpha = 0.5) + facet_grid(Rep_Subj~Rep_theta)
c

#### REML versus DL MSE ----

a <- Normal.Sim.Results[, .(RE_diff = (REML_Estimate - theta)^2 - (DL_Estimate - theta)^2)]
summary(a$RE_diff)
hist(a$RE_diff)

b <- Normal.Sim.Results[, .(RE_diff = mean((REML_Estimate - theta)^2 - (DL_Estimate - theta)^2, na.rm = TRUE) ), 
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
summary(b$RE_diff)

b <- Normal.Sim.Results[, .(RE_diff = mean((REML_Estimate - theta)^2 - (DL_Estimate - theta)^2, na.rm = TRUE),
                           Better = ifelse( mean((REML_Estimate - theta)^2, na.rm = T) < mean((DL_Estimate - theta)^2, na.rm = T), "REML", "DL") ), 
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

summary(as.factor(b$Better))

plot(b$Rep_theta, b$RE_diff)
plot(as.factor(b$Rep_Subj), b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)


c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point(alpha = 0.5) + facet_grid(Rep_Subj~Rep_theta)
c

#### REML versus DL Coverage ----
summary(Total.results$Coverage_REML)
summary(Total.results$Coverage_DL)

a <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL)]
# Plain CI
summary(a$RE_diff)

Total.results[which.max(a$RE_diff)]

b <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL), by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

plot(b$Rep_theta, b$RE_diff)
plot(as.factor(b$Rep_Subj), b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_REML - 0.95) < abs(Coverage_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

summary(as.factor(d$Better))

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e


# HC CI
summary(a$HC_diff)

Total.results[which.max(a$HC_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = HC_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_HC_REML - 0.95) < abs(Coverage_HC_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

summary(as.factor(d$Better))

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e


# KH CI
summary(a$KH_diff)

Total.results[which.max(a$KH_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = KH_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_KH_REML - 0.95) < abs(Coverage_KH_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

summary(as.factor(d$Better))


e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e


# ? No individual level differences possible - wouldn't be meaningful



########################################################### Depreciated

#### REML versus DL Bias ----

summary(Total.results$Bias_REML)
summary(Total.results$Bias_DL)

a <- Normal.Sim.Results[, .(RE_diff = REML_Estimate - DL_Estimate)]
summary(a$RE_diff)
hist(a$RE_diff)

Normal.Sim.Results[which.max(a$RE_diff)]

b <- Normal.Sim.Results[, .(RE_diff = mean(REML_Estimate - DL_Estimate, na.rm = TRUE)), by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
summary(b$RE_diff)

plot(b$Rep_theta, b$RE_diff)
plot(as.factor(b$Rep_Subj), b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)

#### REML versus DL MSE ----

a <- Normal.Sim.Results[, .(RE_diff = (REML_Estimate - theta)^2 - (DL_Estimate - theta)^2)]
summary(a$RE_diff)
hist(a$RE_diff)

b <- Normal.Sim.Results[, .(RE_diff = mean((REML_Estimate - theta)^2 - (DL_Estimate - theta)^2, na.rm = TRUE)), by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
summary(b$RE_diff)

#### REML versus DL Coverage ----
summary(Total.results$Coverage_REML)
summary(Total.results$Coverage_DL)

a <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL)]
# Plain CI
summary(a$RE_diff)

Total.results[which.max(a$RE_diff)]

b <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                      HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                      KH_diff = Coverage_KH_REML - Coverage_KH_DL), by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

# plot(b$Rep_theta, b$RE_diff)
# plot(b$Rep_Subj, b$RE_diff)
# plot(b$Rep_NumStudies, b$RE_diff)
# plot(b$Rep_tau.sq, b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c

# HC CI
summary(a$HC_diff)

Total.results[which.max(a$HC_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = HC_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c


# KH CI
summary(a$KH_diff)

Total.results[which.max(a$KH_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = KH_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(~Rep_theta)
c

# ? No individual level differences possible - wouldn't be meaningful