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

#### UMD Select subset of data for analysis ----

sig.level <- (1 - 0.05/2)

#An.Cond <- Normal.Sim.Results[Rep_theta == theta[2] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 4.2]

#### Loop to get images

#### Specify type of bias and directory
getwd()
Type.of.bias <- "Outcome"
mainDir <- file.path(getwd(), Type.of.bias)


Subj <- c(4.2, 60, 20, 250)

### Summary Bias across all theta - assumes equal effect across theta

dir.create(file.path(mainDir), recursive = TRUE, showWarnings = FALSE)

subDir <- file.path("AcrossAllTheta")

subDir2 <- file.path("Summary")

dir.create(file.path(mainDir, subDir, subDir2), recursive = TRUE, showWarnings = FALSE)

setwd(file.path(mainDir, subDir, subDir2))

##

Bias.values.total <- Normal.Sim.Results[, .(FE = FE_Estimate - Rep_theta,
                                 DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                             by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("Bias across all theta")

p16

ggsave(paste("MD", Type.of.bias, "BiasAcrossAllThetabySubj", "Size", c, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")

## MSE across all theta

MSE1.values <- Normal.Sim.Results[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                           DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                       by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

#MSE1.values <- MSE1.values[, Best := ifelse( (Moreno < DL & Moreno < FE), "Moreno", ifelse( (DL < FE), "DL", "FE"  ) )]

MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p15 <- ggplot(MSE1.values2[variable != "Moreno"], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point( size = 2, position=position_dodge(width = .5)) +   #theme_bw()+
  xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("MSE across all thetas")#+ coord_cartesian(ylim = c(0, 1))

p15

ggsave(paste("MD", Type.of.bias, "MSEAcrossAllTheta", "Size", c, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
       height = 5.67, units = "in")

## Coverage across all theta


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

Coverage.values <- Normal.Sim.Results[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                               DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                               "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                               "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                               "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                               IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                               Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                               Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1))  #+theme_bw()+  scale_colour_grey()

p19

ggsave(paste("MD", Type.of.bias, "CovAcrossAllTheta", "Size", c, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


## Without REML

p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 0.9, size = 2, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1)) #+  #+theme_bw()+  scale_colour_grey()
  # theme(axis.line = element_line(colour = "black"),
  #       #panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       #panel.border = element_blank(),
  #       panel.background = element_blank())
  
p20


ggsave(paste("MD", Type.of.bias, "CovAcrossAllThetaNoREML", "Size", c, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


############################ Not completed - need to consider actually which graphs have value


for (c in Subj){
  
  dir.create(file.path(mainDir), recursive = TRUE, showWarnings = FALSE)
  
  subDir <- file.path("AcrossAllTheta")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  ## once established that theta has no effect, we can use graphs averaged across all thetas
  
  An.Cond <- Normal.Sim.Results[Rep_Subj == c]
  
  ## Bias across all theta
  
  #### Removing REML
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies, Rep_tau.sq)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(.~ Rep_tau.sq) + ggtitle("Bias across all theta")
  
  ggsave(paste("MD", Type.of.bias, "BiasAcrossAllTheta", "Size", c, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
 ## MSE across all theta
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies, Rep_tau.sq)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(aes(shape = variable), size = 1)  +  #theme_bw()+
    xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(~ Rep_tau.sq) + ggtitle("MSE across all thetas")#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("MD", Type.of.bias, "MSEAcrossAllTheta", "Size", c, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ## Coverage across all theta
  
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
  An.Cond[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  ### Does Moreno use z-score?
  
  An.Cond[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
  An.Cond[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_NumStudies, Rep_tau.sq)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid( ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1))  #+theme_bw()+  scale_colour_grey()
  
  p19
  
  ggsave(paste("MD", Type.of.bias, "CovAcrossAllTheta", "Size", c, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  ## Without REML
  
  p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid( ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1))  #+theme_bw()+  scale_colour_grey()
  
  p20
  
  
  ggsave(paste("MD", Type.of.bias, "CovAcrossAllThetaNoREML", "Size", c, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  
  for (d in Studies){
  
  ## Bias by tau2 and theta
  
  An.Cond <- Normal.Sim.Results[Rep_Subj == c & Rep_NumStudies == d]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta, REML = REML_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_theta, Rep_tau.sq)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta", "Rep_tau.sq"))
  
  p11 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_theta), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Theta") + ylab("Bias") + 
  theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(.~ Rep_tau.sq)
  
  setwd(file.path(mainDir))
  
  ggsave(paste("MD", Type.of.bias, "BiasBoxplotTheta", "Size", c, "Studies", d, ".png", sep = ""), p11, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  ## MSE by tau2 and theta
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_theta, Rep_tau.sq)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_theta", "Rep_tau.sq"))
  
  # MSE1.plot <- ggplot(MSE1.values2, aes(x = as.factor(Rep_theta), y = value, group = variable, colour = as.factor(Rep_tau.sq))) +
  #   geom_point(aes(shape = variable), size = 1)  +  #theme_bw()+
  #   xlab("Theta") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
  #   scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") + facet_grid(~ as.factor(tau.sq))#+ coord_cartesian(ylim = c(0, 1))

  
  p13 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_theta), y = value, colour = variable)) +
    geom_point(aes(shape = variable), size = 1)  +  #theme_bw()+
    xlab("Theta") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(~ as.factor(Rep_tau.sq))#+ coord_cartesian(ylim = c(0, 1))
  
  
  ggsave(paste("MD", Type.of.bias, "MSEThetabytau2", "Size", c, "Studies", d, ".png", sep = ""), p13, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  
  ## Coverage by tau2 and theta
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
  An.Cond[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  ### Does Moreno use z-score?
  
  An.Cond[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
  An.Cond[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_theta, Rep_tau.sq)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_theta", "Rep_tau.sq"))
  
  # Coverage.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) +
  #   ylab("Coverage") + xlab("Number of studies") + coord_cartesian(ylim = c(0, 1)) + scale_colour_discrete(name = "Estimator")
  # #Coverage.plot
  
  
  p17 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_theta), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + 
    xlab("Theta")  + ylab("Coverage") + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(~ as.factor(Rep_tau.sq))  + coord_cartesian(ylim = c(0, 1))  # + theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "Cov", "Size", c, "Studies", d, ".png", sep = ""), p17, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  }
  
for (a in theta){
  
  ## Bias by tau2 and number of studies
  
  An.Cond <- Normal.Sim.Results[Rep_Subj == c & Rep_theta == a]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta, REML = REML_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies, Rep_tau.sq)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p12 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(.~ Rep_tau.sq)
  
  subDir <- file.path(paste("Theta", a, sep =""))
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  ggsave(paste("MD", Type.of.bias, "BiasBoxplotStudies", "Size", c, "Theta", a, ".png", sep = ""), p12, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  ## MSE by tau2 and number of studies (separated by theta)
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies, Rep_tau.sq)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p14 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(aes(shape = variable), size = 1)  +  #theme_bw()+
    xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(~ Rep_tau.sq)#+ coord_cartesian(ylim = c(0, 1))
  
  #p14
  
  ggsave(paste("MD", Type.of.bias, "MSENumberbyTau2", "Size", c, "Theta", a, ".png", sep = ""), p14, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ## Coverage by tau2 and number of studies (separated by theta)
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
  An.Cond[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  ### Does Moreno use z-score?
  
  An.Cond[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
  An.Cond[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_NumStudies, Rep_tau.sq)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq"))
  
  p18 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  +  
    xlab("Theta")  + ylab("Number of Studies") + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1))  # + theme_bw()+ scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CovNumberbyTau2", "Size", c, "Theta", a, ".png", sep = ""), p18, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  for (b in tau.sq){
    
      An.Cond <- Normal.Sim.Results[Rep_theta == a & Rep_tau.sq == b & Rep_Subj == c]
      
      An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
      An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
      
      An.Cond[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
      An.Cond[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]
      
      An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
      An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
      
      ### Does Moreno use z-score?
      
      An.Cond[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
      An.Cond[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]
      
      An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
      An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
      
      
      An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
      An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
      
      # An.Cond[, .(Bias = mean(FE_Estimate) - Rep_theta, MSE1 = mean((FE_Estimate - Rep_theta)^2), 
      #             MSE2 = (mean(FE_Estimate) - Rep_theta) + var(FE_Estimate), 
      #             Coverage = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub))), by = .(Rep_theta, Rep_NumStudies)]
      # 
      Bias.values <- An.Cond[, .(FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                                 DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
                             by = .(Rep_theta, Rep_NumStudies)]
      
      
      Bias.values2 <- melt(Bias.values, id = c("Rep_theta", "Rep_NumStudies"))
      
      ### Bias plot
      
      #bias.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Bias.values2) + xlab("Number of Studies") + ylab("Bias")
      bias.plot <- ggplot(Bias.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + 
        geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + 
        ylab("Bias") + scale_colour_grey() + scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") + coord_cartesian(ylim = c(-0.75, 0.75))
      
      # Bias distribution
      
      Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta, REML = REML_Estimate - Rep_theta,
                                       DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                                   by = .(Rep_theta, Rep_NumStudies)]
      
      
      Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta", "Rep_NumStudies"))
      
      
      
      p10 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable)) +
        #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
        geom_boxplot(alpha=1, outlier.shape = NA) +
        ggtitle("Bias") + xlab("Number of studies") + ylab("Bias") + 
        theme_bw() +
        coord_cartesian(ylim = c(-0.75,0.75))
      
      
      MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                                 DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                             by = .(Rep_theta, Rep_NumStudies)]
      
      MSE1.values2 <- melt(MSE1.values, id = c("Rep_theta", "Rep_NumStudies"))
      
      #MSE1.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = MSE1.values2) + xlab("Number of Studies") + ylab("MSE")+ coord_cartesian(ylim = c(0, 2)) 
      MSE1.plot <- ggplot(MSE1.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + geom_line(aes(linetype = variable), size = 1) + theme_bw() + 
        xlab("Number of Studies") + ylab("MSE") + scale_colour_grey() + 
        scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") + coord_cartesian(ylim = c(0, 1))
      
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
      
      # Coverage.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) +
      #   ylab("Coverage") + xlab("Number of studies") + coord_cartesian(ylim = c(0, 1)) + scale_colour_discrete(name = "Estimator")
      # #Coverage.plot
      
      p21 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
        geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
        xlab("Number of studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
         coord_cartesian(ylim = c(0, 1))  #+theme_bw()+  scale_colour_grey()
      
      
      
      subDir <- file.path(paste("Theta", a, sep =""), paste("Tau", b, sep =""))
      
      dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
      
      setwd(file.path(mainDir, subDir))
      
      ggsave(paste("MD", Type.of.bias, "BiasBoxplot", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), p10, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      ggsave(paste("MD", Type.of.bias, "Bias", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), bias.plot, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      ggsave(paste("MD", Type.of.bias, "MSE", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), MSE1.plot, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      ggsave(paste("MD", Type.of.bias, "Cov", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), p21, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      
      
    }
    
  }
  
}




#### Extra checking

summary(Bias.values.total)

Av.Bias <- Normal.Sim.Results[, .(FE = mean(FE_Estimate - Rep_theta, na.rm = T),
                                  DL = mean(DL_Estimate - Rep_theta, na.rm = T), Moreno = mean(Moreno_Estimate - Rep_theta, na.rm = T)),
                              by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

Av.Bias2 <- melt(Av.Bias, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p16 <- ggplot(Av.Bias2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_point(position=position_dodge(width = .5)) + 
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.025,0.025)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("Bias across all theta")

p16
