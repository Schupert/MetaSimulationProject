#### Analysis

### Remove previous variables
rm(list = ls())

Global_list_directory <- c("D:/Stats/AFP/Pilot/NormalAltOutcome1",
                           "D:/Stats/AFP/Pilot/NormalBegg5",
                           "D:/Stats/AFP/Pilot/NormalMethods5",
                           "D:/Stats/AFP/Pilot/NormalNoBias6",
                           "D:/Stats/AFP/Pilot/NormalOutcome5",
                           "D:/Stats/AFP/Pilot/NormalStep5")

Global_list_filenames <- c("NSTotalOutRDS",
                           "NSBeggRDS",
                           "NSTotalMethRDS",
                           "NSTotalV2RDS",
                           "NSTotalOutRDS",
                           "NSStepRDS")

Global_list_label <- c("AltOut",
                       "ModBegg",
                       "Method",
                       "None",
                       "Out",
                       "Step")

Global_list_formalname <- c("Adjusted Outcome Bias",
                            "Modified Begg Bias",
                            "Methodological Bias",
                            "No Bias",
                            "Outcome Bias",
                            "Step Bias")

###  Loop ----

for (xi in 1:1){
  
  rm(list=setdiff(ls(), c("Global_list_directory", "Global_list_filenames", "Global_list_formalname", "xi")))

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

sig.level <- (1 - 0.05/2)
  
FileLoc <- Global_list_directory[xi]
setwd(paste(FileLoc))
  
  
  


#### UMD Import data here ----

system.time(Normal.Sim.Results <- readRDS(file = paste(Global_list_filenames[xi]) ))

Normal.Sim.Results <- data.table(Normal.Sim.Results)

#### Renaming variables ----

# Theta in log conditions

# Subject in all conditions
Subj <- c("Empirical", "Small", "Fixed", "Large")
Normal.Sim.Results[, Rep_Subj := factor(Rep_Subj, labels = c("Empirical", "Small", "Fixed", "Large"))]


#### Specify type of bias and directory ----
getwd()
Type.of.bias <- Global_list_formalname[xi]
mainDir <- file.path(getwd(), Type.of.bias)

dir.create(file.path(mainDir), recursive = TRUE, showWarnings = FALSE)


for (c in Subj){
  
  subDir <- file.path("SubjectHeld")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  An.Cond <- Normal.Sim.Results[Rep_Subj == c]
  
  
  ### Bias
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_theta"))
  
  p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(Rep_theta~ Rep_tau.sq) + 
    ggtitle(paste("Bias in", c, "studies:", Type.of.bias, sep = " "))
  
  ggsave(paste("MD", Type.of.bias, "BiasByAll", "Size", c, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ### MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_theta"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(aes(shape = variable))  +  #theme_bw()+
    xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(Rep_theta~ Rep_tau.sq) + ggtitle(paste("MSE in", c, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  
  ggsave(paste("MD", Type.of.bias, "MSEByAll", "Size", c, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  
  ### Coverage
  
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
  by = .(Rep_NumStudies, Rep_tau.sq, Rep_theta)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_theta"))
  
  p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage in", c, "studies:", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_theta ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAll", "Size", c, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  ## Without REML
  
  p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage in", c, "studies:", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_theta ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAllNoREML", "Size", c, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
}


for (a in theta){
  
  subDir <- file.path("ThetaHeld")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  An.Cond <- Normal.Sim.Results[Rep_theta == a]
  
  
  ### Bias
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))
  
  p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(Rep_Subj~ Rep_tau.sq) + 
    ggtitle(paste("Bias with theta", a, ":", Type.of.bias, sep = " "))
  
  ggsave(paste("MD", Type.of.bias, "BiasByAll", "Theta", a, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ### MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(aes(shape = variable))  +  #theme_bw()+
    xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with theta", a, ":", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  
  ggsave(paste("MD", Type.of.bias, "MSEByAll", "Theta", a, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  
  ### Coverage
  
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
  by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))
  
  p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with theta", a, ":", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAll", "Theta", a, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  ## Without REML
  
  p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with theta", a, ":", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAllNoREML", "Theta", a, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
}

for (e in Studies){
  
  subDir <- file.path("NumberofStudiesHeld")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  An.Cond <- Normal.Sim.Results[Rep_NumStudies == e]
  
  
  ### Bias
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_theta, Rep_tau.sq, Rep_Subj)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta", "Rep_tau.sq", "Rep_Subj"))
  
  p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_theta), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Theta") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(Rep_Subj~ Rep_tau.sq) + 
    ggtitle(paste("Bias with", e, "studies:", Type.of.bias, sep = " "))
  
  ggsave(paste("MD", Type.of.bias, "BiasByAll", "NumStudies", e, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ### MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_theta, Rep_tau.sq, Rep_Subj)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_theta", "Rep_tau.sq", "Rep_Subj"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_theta), y = value, colour = variable)) +
    geom_point(aes(shape = variable))  +  #theme_bw()+
    xlab("Theta") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  
  ggsave(paste("MD", Type.of.bias, "MSEByAll", "NumStudies", e, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  
  ### Coverage
  
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
  by = .(Rep_theta, Rep_tau.sq, Rep_Subj)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_theta", "Rep_tau.sq", "Rep_Subj"))
  
  p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_theta), y = value, colour = variable)) +
    geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with", e, "studies:", Type.of.bias, sep = " ")) + 
    xlab("Theta")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAll", "NumStudies", e, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  ## Without REML
  
  p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_theta), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with", e, "studies:", Type.of.bias, sep = " ")) + 
    xlab("Theta")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAllNoREML", "NumStudies", e, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
}

for (b in tau.sq){
  
  subDir <- file.path("Tau2Held")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  An.Cond <- Normal.Sim.Results[Rep_tau.sq == b]
  
  
  ### Bias
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies, Rep_theta, Rep_Subj)]
  
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_theta", "Rep_Subj"))
  
  p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75))+
    facet_grid(Rep_Subj~ Rep_theta) + 
    ggtitle(paste("Bias with tau2", b, ":", Type.of.bias, sep = " "))
  
  ggsave(paste("MD", Type.of.bias, "BiasByAll", "Tau2", b, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  ### MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies, Rep_theta, Rep_Subj)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_theta", "Rep_Subj"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(aes(shape = variable))  +  #theme_bw()+
    xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
    facet_grid(Rep_Subj~ Rep_theta) + ggtitle(paste("MSE with tau2", b, ":", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  
  ggsave(paste("MD", Type.of.bias, "MSEByAll", "Tau2", b, ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
         height = 5.67, units = "in")
  
  
  ### Coverage
  
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
  by = .(Rep_NumStudies, Rep_theta, Rep_Subj)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_theta", "Rep_Subj"))
  
  p19 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with tau2", b, ":", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_theta)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAll", "Tau2", b, ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
  
  ## Without REML
  
  p20 <-  ggplot(Coverage.values2[!(variable %like% "REML"),], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
    geom_point(alpha = 0.9, position=position_dodge(width=0.5))  + ggtitle(paste("Coverage with tau2", b, ":", Type.of.bias, sep = " ")) + 
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    facet_grid(Rep_Subj ~ as.factor(Rep_theta)) + coord_cartesian(ylim = c(0.5, 1))  #+theme_bw()+  scale_colour_grey()
  
  
  ggsave(paste("MD", Type.of.bias, "CoverageByAllNoREML", "Tau2", b, ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
  
}

}
