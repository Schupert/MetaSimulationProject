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

system.time(Normal.Sim.Results <- readRDS(file = "NSTotalMethRDS"))

Normal.Sim.Results <- data.table(Normal.Sim.Results)

#### UMD Select subset of data for analysis ----

sig.level <- (1 - 0.05/2)

#An.Cond <- Normal.Sim.Results[Rep_theta == theta[2] & Rep_tau.sq == tau.sq[1] & Rep_Subj == 4.2]

#### Loop to get images

#### Specify type of bias and directory
getwd()
Type.of.bias <- "Method"
mainDir <- file.path(getwd(), Type.of.bias)

Subj <- c(4.2, 60, 20, 250)

for (a in theta){
  for (b in tau.sq){
    for (c in Subj){
      An.Cond <- Normal.Sim.Results[Rep_theta == a & Rep_tau.sq == b & Rep_Subj == c]
      
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
      
      An.Cond[, .(Bias = mean(FE_Estimate) - Rep_theta, MSE1 = mean((FE_Estimate - Rep_theta)^2), 
                  MSE2 = (mean(FE_Estimate) - Rep_theta) + var(FE_Estimate), 
                  Coverage = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub))), by = .(Rep_theta, Rep_NumStudies)]
      
      Bias.values <- An.Cond[, .(FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                                 DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta),
                             by = .(Rep_theta, Rep_NumStudies)]
      
      
      Bias.values2 <- melt(Bias.values, id = c("Rep_theta", "Rep_NumStudies"))
      
      ### Bias plot
      
      #bias.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Bias.values2) + xlab("Number of Studies") + ylab("Bias")
      bias.plot <- ggplot(Bias.values2, aes(x = Rep_NumStudies, y = value, group = variable)) + 
        geom_line(aes(linetype = variable), size = 1) + theme_bw() + xlab("Number of Studies") + 
        ylab("Bias") + scale_colour_grey() + scale_linetype_manual(values=c("solid", "dotted", "dotdash", "longdash"), name = "Estimator") + coord_cartesian(ylim = c(-0.75, 0.75))
      
      
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
      
      Coverage.plot <- qplot(Rep_NumStudies, value, colour = variable, geom = "line", data = Coverage.values2) +
        ylab("Coverage") + xlab("Number of studies") + coord_cartesian(ylim = c(0, 1)) + scale_colour_discrete(name = "Estimator")
      #Coverage.plot
      
      subDir <- file.path(paste("Theta", a, sep =""), paste("Tau", b, sep =""))
      
      dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
      
      setwd(file.path(mainDir, subDir))
      
      ggsave(paste("MD", Type.of.bias, "Bias", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), bias.plot, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      ggsave(paste("MD", Type.of.bias, "MSE", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), MSE1.plot, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      ggsave(paste("MD", Type.of.bias, "Cov", "Th", a, "Tau", b, "Size", c, ".png", sep = ""), Coverage.plot, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
      
      
    }
    
  }
  
}