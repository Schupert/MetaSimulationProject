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

#### Functions -----

CI.betw <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 1, 0)
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

#### LOR Declare variables ----

# Reps = number of repetitions of experiment
Reps = 10000

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

#### UMD Import data here ----

system.time(Normal.Sim.Results <- readRDS(file = "NSTotalV2RDS"))

Normal.Sim.Results <- data.table(Normal.Sim.Results)

#### UMD Select subset of data for analysis ----


sig.level <- (1 - 0.05/2)

An.Cond <- Normal.Sim.Results[Rep_theta == theta[3] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 250]

#### LOR Import data here ----

system.time(LogOR.Sim.Results <- readRDS(file = "LSTotalBeggRDS"))

LogOR.Sim.Results <- data.table(LogOR.Sim.Results)
#### LOR Select subset of data for analysis ----


sig.level <- (1 - 0.05/2)

An.Cond <- LogOR.Sim.Results[Rep_theta == theta[5] & Rep_tau.sq == tau.sq[4] & Rep_Subj == 20 & Rep_ev_freq == 0.5]

#### Calculate CI values -----

## There may be cleaner way to do this with data.table
An.Cond$FE_CIlb <- An.Cond$FE_Estimate - qnorm(sig.level) * An.Cond$FE_se
An.Cond$FE_CIub <- An.Cond$FE_Estimate + qnorm(sig.level) * An.Cond$FE_se

An.Cond$REML_CIlb <- An.Cond$REML_Estimate - qnorm(sig.level) * An.Cond$REML_se
An.Cond$REML_CIub <- An.Cond$REML_Estimate + qnorm(sig.level) * An.Cond$REML_se

An.Cond$DL_CIlb <- An.Cond$DL_Estimate - qnorm(sig.level) * An.Cond$DL_se
An.Cond$DL_CIub <- An.Cond$DL_Estimate + qnorm(sig.level) * An.Cond$DL_se

An.Cond$Doi_CIlb <- An.Cond$FE_Estimate - qnorm(sig.level) * An.Cond$HC_DL_se
An.Cond$Doi_CIub <- An.Cond$FE_Estimate + qnorm(sig.level) * An.Cond$HC_DL_se

### Does Moreno use z-score?

An.Cond$Moreno_CIlb <- An.Cond$Moreno_Estimate - qnorm(sig.level) * An.Cond$Moreno_se
An.Cond$Moreno_CIub <- An.Cond$Moreno_Estimate + qnorm(sig.level) * An.Cond$Moreno_se

An.Cond$Mult_CIlb <- An.Cond$FE_Estimate - qnorm(sig.level) * An.Cond$Mult_se
An.Cond$Mult_CIub <- An.Cond$FE_Estimate + qnorm(sig.level) * An.Cond$Mult_se


#### Bias ----

# An.Cond[, .(Bias = mean(FE_Estimate) - Rep_theta, MSE1 = mean((FE_Estimate - Rep_theta)^2), 
#             MSE2 = (mean(FE_Estimate) - Rep_theta) + var(FE_Estimate), 
#             Coverage = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub))), by = .(Rep_theta, Rep_NumStudies)]

Bias.values <- An.Cond[, .(FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                           DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
                       by = .(Rep_theta, Rep_NumStudies)]


Bias.values2 <- melt(Bias.values, id = c("Rep_theta", "Rep_NumStudies"))

### Bias plot

bias.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Bias.values2) + xlab("Number of Studies") + ylab("Bias")
bias.plot <- ggplot(Bias.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + 
  geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + 
  ylab("Bias") + scale_colour_grey() + scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") #+ coord_cartesian(ylim = c(-0.5, 0.25))

bias.plot


#### MSE ----

MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                           DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                       by = .(Rep_theta, Rep_NumStudies)]

MSE1.values2 <- melt(MSE1.values, id = c("Rep_theta", "Rep_NumStudies"))

MSE1.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = MSE1.values2) + xlab("Number of Studies") + ylab("MSE")+ coord_cartesian(ylim = c(0, 2)) 
MSE1.plot <- ggplot(MSE1.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + 
  xlab("Number of Studies") + ylab("MSE") + scale_colour_grey() + 
  scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") + coord_cartesian(ylim = c(0, 0.75))

MSE1.plot

MSE1.plot <- ggplot(MSE1.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + 
  xlab("Number of Studies") + ylab("MSE") + scale_colour_grey() + coord_cartesian(ylim = c(0, 0.5)) + geom_point(aes(shape = variable))
MSE1.plot

MSE2.values <- An.Cond[, .(FE = (mean(FE_Estimate) - Rep_theta)^2 + var(FE_Estimate), REML = (mean(REML_Estimate, na.rm = TRUE) - Rep_theta)^2 + var(REML_Estimate, na.rm = TRUE),
                           DL = (mean(DL_Estimate, na.rm = TRUE) - Rep_theta)^2 + var(DL_Estimate, na.rm = TRUE), Moreno = (mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta)^2 + var(Moreno_Estimate, na.rm = TRUE) ),
                       by = .(Rep_theta, Rep_NumStudies)]

MSE2.values <- melt(MSE2.values, id = c("Rep_theta", "Rep_NumStudies"))

MSE2.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = MSE2.values) #+ coord_cartesian(ylim = c(0, 0.1))
MSE2.plot


#### Coverage ---- 

Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                               DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                               "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                               "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                               "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                               IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                               Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                               Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_theta", "Rep_NumStudies"))

Coverage.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) +
  ylab("Coverage") + xlab("Number of studies") + coord_cartesian(ylim = c(0, 1)) + scale_colour_discrete(name = "Estimator")
Coverage.plot

Coverage.plot <- ggplot(Coverage.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() +
  ylab("Coverage") + xlab("Number of studies") + scale_colour_grey() + coord_cartesian(ylim = c(0.9, 1))
Coverage.plot

Coverage.plot <- ggplot(Coverage.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() +
  ylab("Coverage") + xlab("Number of studies") + scale_colour_grey() + coord_cartesian(ylim = c(0.9, 1)) + geom_point(aes(shape = variable))
Coverage.plot



#### Testing stargazer -----

Coverage.values <- data.table(Coverage.values[,2, with = FALSE], Coverage.values[, -1:-2, with = FALSE ]*100)
stargazer(Coverage.values, summary = FALSE, rownames = FALSE, digits = 2)

stargazer(Bias.values[, -1, with = FALSE ], summary = FALSE, rownames = FALSE)

stargazer(MSE1.values[, -1, with = FALSE ], summary = FALSE, rownames = FALSE)



#### Above and below coverage -----

CI.updown <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 0, ifelse(a < b, -1, 1))
  return(output)
}

Outside.values <- An.Cond[, .(FE = mean(CI.updown(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.updown(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                              DL = mean(CI.updown(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.updown(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                              "HC REML" = mean(CI.updown(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                              "KH DL" = mean(CI.updown(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                              "KH REML" = mean(CI.updown(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                              IVHet = mean(CI.updown(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                              Moreno = mean(CI.updown(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                              Mult = mean(CI.updown(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies)]

Outside.values2<- melt(Outside.values, id = c("Rep_theta", "Rep_NumStudies"))
#Outside.values

stargazer(Outside.values[, -1, with = FALSE ], summary = FALSE, rownames = FALSE, digits = 2)


###### Saving out ----

ggsave("LOBeggBiasTh0Tau004ER5Emp.png", bias.plot, dpi = 300, device = "png")
ggsave("LOBeggMSETh0Tau004ER5Emp.png", MSE1.plot, dpi = 300, device = "png")
ggsave("LOBeggCovTh0Tau004ER5Emp.png", Coverage.plot, dpi = 300, device = "png")



################################# Setting up 95% intervals ----


# Bias.values <- An.Cond[, .(FE = mean(FE_Estimate) - Rep_theta, REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
#                            DL = mean(DL_Estimate) - Rep_theta, Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
#                        by = .(Rep_theta, Rep_NumStudies)]

Bias.values <- An.Cond[, .(FE = mean(FE_Estimate) - Rep_theta, FE95 = quantile(FE_Estimate, probs = 0.95) - Rep_theta, FE5 = quantile(FE_Estimate, probs = 0.05)- Rep_theta) ,
                       by = .(Rep_theta, Rep_NumStudies)]


Bias.values <- melt(Bias.values, id = c("Rep_theta", "Rep_NumStudies"))

### Bias plot

bias.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Bias.values) + xlab("Number of Studies") + ylab("Bias")
bias.plot <- ggplot(Bias.values, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + ylab("Bias") + scale_colour_grey()
bias.plot

#### UMD AOV Attempt -----
# Potentially problematic as it creates large ~1.6Gb element
# Consider not saving output? Just printing

Normal.Sim.Results$FE_Bias <- Normal.Sim.Results$FE_Estimate - Normal.Sim.Results$Rep_theta
system.time(aov1 <- aov(FE_Bias ~ as.factor(Rep_NumStudies) + 
                          as.factor(Rep_tau.sq)+ as.factor(Rep_Subj) + as.factor(Rep_theta), data = Normal.Sim.Results))

summary(aov1)

# Bartlett Test of Homogeneity of Variances
bartlett.test(FE_Bias ~ as.factor(Rep_NumStudies), data=Normal.Sim.Results)

# Figner-Killeen Test of Homogeneity of Variances
fligner.test(FE_Bias ~ as.factor(Rep_NumStudies) , data=Normal.Sim.Results)

rm(aov1)
model1 <- lm(FE_Bias ~ Rep_NumStudies + 
              Rep_tau.sq + Rep_Subj + Rep_theta, data = Normal.Sim.Results)

summary(model1)

rm(model1)
## This breaks R

layout(matrix(c(1,2,3,4),2,2)) # optional layout 
plot(aov1) # diagnostic plots

##### UMD Model attempt -----

## FE Bias

Normal.Sim.Results$FE_Bias <- Normal.Sim.Results$FE_Estimate - Normal.Sim.Results$Rep_theta

model1 <- lm(FE_Bias ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta), data = Normal.Sim.Results)

summary(model1)
stargazer(model1)

rm(model1)

## DL Bias

Normal.Sim.Results$DL_Bias <- Normal.Sim.Results$DL_Estimate - Normal.Sim.Results$Rep_theta

model1 <- lm(DL_Bias ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta), data = Normal.Sim.Results)

summary(model1)
stargazer(model1)

rm(model1)

## FE MSE

Normal.Sim.Results$FE_MSE <- (Normal.Sim.Results$FE_Estimate - Normal.Sim.Results$Rep_theta)^2

model2 <- lm(FE_MSE ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta), data = Normal.Sim.Results)

summary(model2)
stargazer(model2, no.space=TRUE, single.row=TRUE)

rm(model2)

## DL MSE

Normal.Sim.Results$DL_MSE <- (Normal.Sim.Results$DL_Estimate - Normal.Sim.Results$Rep_theta)^2

model2 <- lm(DL_MSE ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta), data = Normal.Sim.Results)

summary(model2)
stargazer(model2, no.space=TRUE, single.row=TRUE)

rm(model2)


### Attempt including type of metric
Normal.Sim.Results$FE_MSE <- (Normal.Sim.Results$FE_Estimate - Normal.Sim.Results$Rep_theta)^2
Normal.Sim.Results$DL_MSE <- (Normal.Sim.Results$DL_Estimate - Normal.Sim.Results$Rep_theta)^2
Normal.Sim.Results$REML_MSE <- (Normal.Sim.Results$REML_Estimate - Normal.Sim.Results$Rep_theta)^2
Normal.Sim.Results$Moreno_MSE <- (Normal.Sim.Results$Moreno_Estimate - Normal.Sim.Results$Rep_theta)^2

dummy.table <- Normal.Sim.Results[,.(Rep_NumStudies, Rep_theta, Rep_Subj, Rep_tau.sq, FE_MSE, DL_MSE, REML_MSE, Moreno_MSE)]
dummy.table <- melt(dummy.table, id = c("Rep_theta", "Rep_NumStudies", "Rep_Subj", "Rep_tau.sq"))

model2 <- lm(value ~ as.factor(Rep_NumStudies) + as.factor(variable) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta), data = dummy.table)

summary(model2)
stargazer(model2, no.space=TRUE, single.row=TRUE)

rm(model2)
rm(dummy.table)

##### LOR Model attempt -----

## Bias

LogOR.Sim.Results$FE_Bias <- LogOR.Sim.Results$FE_Estimate - LogOR.Sim.Results$Rep_theta

model1 <- lm(FE_Bias ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta) + as.factor(Rep_ev_freq), data = LogOR.Sim.Results)

summary(model1)
stargazer(model1)

rm(model1)

## MSE

LogOR.Sim.Results$FE_MSE <- (LogOR.Sim.Results$FE_Estimate - LogOR.Sim.Results$Rep_theta)^2

model2 <- lm(FE_MSE ~ as.factor(Rep_NumStudies) + 
               as.factor(Rep_tau.sq) + as.factor(Rep_Subj) + as.factor(Rep_theta) + as.factor(Rep_ev_freq), data = LogOR.Sim.Results)

summary(model2)
stargazer(model2)

#### Moreno estimate with other standard error ----

An.Cond$MorenoAlt_CIlb <- An.Cond$Moreno_Estimate - qnorm(sig.level) * An.Cond$HC_DL_se
An.Cond$MorenoAlt_CIub <- An.Cond$Moreno_Estimate + qnorm(sig.level) * An.Cond$HC_DL_se

An.Cond$MorenoAlt2_CIlb <- An.Cond$DL_Estimate - qnorm(sig.level) * An.Cond$Moreno_se
An.Cond$MorenoAlt2_CIub <- An.Cond$DL_Estimate + qnorm(sig.level) * An.Cond$Moreno_se

Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), #REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                               DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                               #"HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                               "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                               #"KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                               #IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                               Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                               #Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                              "Moreno est HC se" = mean(CI.betw(Rep_theta, MorenoAlt_CIlb, MorenoAlt_CIub), na.rm = TRUE),
                              "Moreno se DL est" = mean(CI.betw(Rep_theta, MorenoAlt2_CIlb, MorenoAlt2_CIub), na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_theta", "Rep_NumStudies"))

Coverage.plot.extra <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) +
  ylab("Coverage") + xlab("Number of studies") + coord_cartesian(ylim = c(0, 1)) + scale_colour_discrete(name = "Estimator")
Coverage.plot.extra

Coverage.values <- data.table(Coverage.values[,2, with = FALSE], Coverage.values[, -1:-2, with = FALSE ]*100)
stargazer(Coverage.values, summary = FALSE, rownames = FALSE, digits = 2, type = "text")


CI.updown <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 0, ifelse(a < b, -1, 1))
  return(output)
}

Outside.values <- An.Cond[, .(FE = mean(CI.updown(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), #REML = mean(CI.updown(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                              DL = mean(CI.updown(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.updown(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                              #"HC REML" = mean(CI.updown(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                              "KH DL" = mean(CI.updown(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                              #"KH REML" = mean(CI.updown(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                              #IVHet = mean(CI.updown(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                              Moreno = mean(CI.updown(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                              #Mult = mean(CI.updown(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                              "Moreno est HC se" = mean(CI.updown(Rep_theta, MorenoAlt_CIlb, MorenoAlt_CIub), na.rm = TRUE),
                              "Moreno se DL est" = mean(CI.updown(Rep_theta, MorenoAlt2_CIlb, MorenoAlt2_CIub), na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies)]

Outside.values2<- melt(Outside.values, id = c("Rep_theta", "Rep_NumStudies"))
#Outside.values

stargazer(Outside.values[, -1, with = FALSE ], summary = FALSE, rownames = FALSE, digits = 2, type = "text")


#### Checking I2 ----

summary(Normal.Sim.Results[Rep_tau.sq == tau.sq[4]]$DL_I2)
hist(Normal.Sim.Results[Rep_tau.sq == tau.sq[4]]$DL_I2)

summary(LogOR.Sim.Results[Rep_tau.sq == tau.sq[4]]$DL_I2)
hist(LogOR.Sim.Results[Rep_tau.sq == tau.sq[4]]$DL_I2)

##### Early testing ----

