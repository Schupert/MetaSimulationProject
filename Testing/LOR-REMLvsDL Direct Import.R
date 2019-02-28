#### Looped Analysis UMD

### Remove previous variables
#rm(list = ls())
rm(list=setdiff(ls(), "LogOR.Complete"))

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

### Direct import

setwd("D:/Stats/AFP/Results/2019/LOR")

Total.results <- data.table(read.csv("LogORResults_Looped.csv"))

#### REML versus DL Bias ----

summary(Total.results$Bias_REML)
summary(Total.results$Bias_DL)

a <- LogOR.Sim.Results[, .(RE_diff = REML_Estimate - DL_Estimate)]
summary(a$RE_diff)
hist(a$RE_diff)

LogOR.Sim.Results[which.max(a$RE_diff)]

b <- LogOR.Sim.Results[, .(RE_diff = mean(REML_Estimate - DL_Estimate, na.rm = TRUE),
                           Better = ifelse(abs(mean(REML_Estimate, na.rm = TRUE)) < abs(mean(DL_Estimate, na.rm = TRUE)), "REML", "DL")),
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]
summary(b$RE_diff)
summary(as.factor(b$Better))

plot(b$Rep_theta, b$RE_diff)
plot(b$Rep_Subj, b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)
plot(b$Rep_ev_freq, b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
c

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point(alpha = 0.5) + facet_grid(Rep_Subj~Rep_theta)
c

#### REML versus DL MSE ----

a <- LogOR.Sim.Results[, .(RE_diff = (REML_Estimate - theta)^2 - (DL_Estimate - theta)^2)]
summary(a$RE_diff)
hist(a$RE_diff)

b <- LogOR.Sim.Results[, .(RE_diff = mean((REML_Estimate - theta)^2 - (DL_Estimate - theta)^2, na.rm = TRUE) ), 
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]
summary(b$RE_diff)

b <- LogOR.Sim.Results[, .(RE_diff = mean((REML_Estimate - theta)^2 - (DL_Estimate - theta)^2, na.rm = TRUE),
                           Better = ifelse( mean((REML_Estimate - theta)^2, na.rm = T) < mean((DL_Estimate - theta)^2, na.rm = T), "REML", "DL") ), 
                       by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

summary(as.factor(b$Better))

plot(b$Rep_theta, b$RE_diff)
plot(as.factor(b$Rep_Subj), b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)
plot(b$Rep_ev_freq, b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
c

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
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL), by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

plot(b$Rep_theta, b$RE_diff)
plot(as.factor(b$Rep_Subj), b$RE_diff)
plot(b$Rep_NumStudies, b$RE_diff)
plot(as.factor(b$Rep_tau.sq), b$RE_diff)
plot(b$Rep_ev_freq, b$RE_diff)

c <- ggplot(b, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_REML - 0.95) < abs(Coverage_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

summary(as.factor(d$Better))

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e


### Extra idea - consider excluding small differences

summary(as.factor(d[abs(RE_diff) > 0.01]$Better))

e <- ggplot(d[abs(RE_diff) > 0.01], aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e


# HC CI
summary(a$HC_diff)

Total.results[which.max(a$HC_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = HC_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_HC_REML - 0.95) < abs(Coverage_HC_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

summary(as.factor(d$Better))

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e

### Extra idea - consider excluding small differences

summary(as.factor(d[abs(RE_diff) > 0.01]$Better))

e <- ggplot(d[abs(RE_diff) > 0.01], aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e

# KH CI
summary(a$KH_diff)

Total.results[which.max(a$KH_diff)]

c <- ggplot(b, aes(x = Rep_NumStudies, y = KH_diff, colour = as.factor(Rep_Subj), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
c

d <- Total.results[, .(RE_diff = Coverage_REML - Coverage_DL,
                       HC_diff = Coverage_HC_REML - Coverage_HC_DL,
                       KH_diff = Coverage_KH_REML - Coverage_KH_DL,
                       Better = ifelse(abs(Coverage_KH_REML - 0.95) < abs(Coverage_KH_DL - 0.95), "REML", "DL")), 
                   by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

summary(as.factor(d$Better))

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e

e <- ggplot(d, aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_Subj~Rep_theta)
e

### Extra idea - consider excluding small differences

summary(as.factor(d[abs(RE_diff) > 0.01]$Better))

e <- ggplot(d[abs(RE_diff) > 0.01], aes(x = Rep_NumStudies, y = RE_diff, colour = as.factor(Better), shape = as.factor(Rep_tau.sq))) +
  geom_point() + facet_grid(Rep_ev_freq~Rep_theta)
e


# ? No individual level differences possible - wouldn't be meaningful