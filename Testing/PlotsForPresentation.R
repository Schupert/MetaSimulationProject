##### Images for presentation
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
library(RColorBrewer)

################## UMD section -----

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

Subj <- c(4.2, 60, 20, 250)

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
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("Bias across all theta")  +
  scale_fill_brewer(type = "qual", palette = "Dark2")

#p16

ggsave(paste("MD", "None", "BiasAcrossAllThetabySubj", ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")

## MSE across all theta

MSE1.values <- Normal.Sim.Results[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                                      DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), Moreno = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                                  by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

#MSE1.values <- MSE1.values[, Best := ifelse( (Moreno < DL & Moreno < FE), "Moreno", ifelse( (DL < FE), "DL", "FE"  ) )]

MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point( size = 2, position=position_dodge(width = .5)) +   #theme_bw()+
  xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("MSE across all thetas")#+ coord_cartesian(ylim = c(0, 1))

#p15

ggsave(paste("MD", "None", "MSEAcrossAllTheta", ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
       height = 5.67, units = "in")

p15 <- ggplot(MSE1.values2[variable != "Moreno"], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point( size = 2, position=position_dodge(width = .5)) +   #theme_bw()+
  xlab("Number of Studies") + ylab("MSE")  + scale_y_log10() + # scale_colour_grey()+
  facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle("MSE across all thetas")#+ coord_cartesian(ylim = c(0, 1))

#p15

ggsave(paste("MD", "None", "MSEAcrossAllTheta", "MorenoExcluded", ".png", sep = ""), p15, dpi = 300, device = "png", width = 8.01, 
       height = 5.67, units = "in")


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
  facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.7, 1))  #+theme_bw()+  scale_colour_grey()

p19

ggsave(paste("MD", "None", "CovAcrossAllTheta", ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")

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

## KH vs HC

Coverage.values.red <- Normal.Sim.Results[, .("HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                          "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                          "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                          "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE)
),
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

Coverage.values2.red<- melt(Coverage.values.red, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj"))

p21 <- ggplot(Coverage.values2.red, aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.85, 1))  #+theme_bw()+  scale_colour_grey()
p21

ggsave(paste("MD", "None", "CovAcrossAllThetaKHvsHC", ".png", sep = ""), p21, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


#### FE vs REML Bias ----
## Uses NO BIAS
# Empirical

An.cond <- Normal.Sim.Results[Rep_Subj == Subj[1]]

Bias.values.total <- Normal.Sim.Results[, .(FE = FE_Estimate - Rep_theta,
                                            DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                                        by = .(Rep_NumStudies, Rep_tau.sq)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq"))

p16 <- ggplot(Bias.values.total2[variable != "Moreno"], aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(~ Rep_tau.sq) + ggtitle("Bias across all theta with empirical study size")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  scale_fill_brewer(type = "qual", palette = "Dark2")

p16

ggsave(paste("MD", "None", "BiasFEDLAcrossAllTheta", "Subj", 4.2, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")

# Including Moreno

p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(~ Rep_tau.sq) + ggtitle("Bias across all theta with empirical study size")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  scale_fill_brewer(type = "qual", palette = "Dark2")

#p16

ggsave(paste("MD", "None", "BiasFEDLMorenoAcrossAllTheta", "Subj", 4.2, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")


#### FE vs REML Coverage ----
## Uses NO BIAS
# Empirical

An.Cond <- Normal.Sim.Results[Rep_Subj == Subj[1]]

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

p19 <-  ggplot(Coverage.values2[variable %in% c("FE", "DL", "KH DL")], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 1, position=position_dodge(width=0.5), size = 3)  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  facet_grid(~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1)) +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[1:3])
  #scale_colour_brewer(type = "qual", palette = "Dark2") #+theme_bw()+  scale_colour_grey()

p19

ggsave(paste("MD", "None", "CovAcrossAllThetaDLFEKH", ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


#### HC IVHet good coverage ----

p20 <-  ggplot(Coverage.values2[variable %in% c("KH DL", "HC DL", "IVHet")], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 1, position=position_dodge(width=0.5), size = 3)  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  facet_grid(~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0, 1))  +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[c(4,3,5)])
  #scale_colour_brewer(type = "qual", palette = "Dark2") #+theme_bw()+  scale_colour_grey()

p20

ggsave(paste("MD", "None", "CovAcrossAllThetaHCIVHet", ".png", sep = ""), p20, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")





################## LOR section -----

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
library(RColorBrewer)

#### Functions -----

CI.betw <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 1, 0)
  return(output)
}

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

#### LOR Import data here ----

system.time(LogOR.Sim.Results <- readRDS(file = "LSTotalStepRDS"))

LogOR.Sim.Results <- data.table(LogOR.Sim.Results)
#### LOR Select subset of data for analysis ----

sig.level <- (1 - 0.05/2)

Subj <- c(4.7, 100, 20, 250)


#An.Cond <- LogOR.Sim.Results[Rep_theta == theta[3] & Rep_tau.sq == tau.sq[3] & Rep_Subj == 4.7 & Rep_ev_freq == 0.5]
An.Cond <- LogOR.Sim.Results[Rep_Subj == 4.7]


Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta, REML = REML_Estimate - Rep_theta,
                                 DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                             by = .(Rep_theta, Rep_ev_freq, Rep_tau.sq)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta", "Rep_ev_freq", "Rep_tau.sq"))

p11 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_theta), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Theta") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75))+
  facet_grid(Rep_tau.sq~ Rep_ev_freq)

p11

ggsave(paste("OR", "None", "Bias", "ThetabyERandTau2", "Size", "4.7", ".png", sep = ""), p11, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


## Coverage across all theta


LogOR.Sim.Results[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
LogOR.Sim.Results[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]

LogOR.Sim.Results[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
LogOR.Sim.Results[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]

LogOR.Sim.Results[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
LogOR.Sim.Results[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]

### Does Moreno use z-score?

LogOR.Sim.Results[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
LogOR.Sim.Results[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]

LogOR.Sim.Results[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
LogOR.Sim.Results[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]


LogOR.Sim.Results[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
LogOR.Sim.Results[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]

Coverage.values <- LogOR.Sim.Results[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                               DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                               "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                               "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                               "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                               IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                               Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                               Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq"))

# Not controlled for theta.

p19 <-  ggplot(Coverage.values2[Rep_ev_freq == 0.5], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  facet_grid(Rep_Subj ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.7, 1))  #+theme_bw()+  scale_colour_grey()

p19


ggsave(paste("OR", "None", "Cov", "AcrossAllTheta", "ER", "0.5", ".png", sep = ""), p19, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")



Coverage.values <- LogOR.Sim.Results[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                         DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                         "HC REML" = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                         "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                         "KH REML" = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                         IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                         Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                         Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
),
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq, Rep_theta)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq", "Rep_theta"))


Coverage.values[, Rep_theta := factor(Rep_theta, labels = c("- log(4)", "- log(5/4)", "log(1)", "log (5/4)", "log(4)"))]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq", "Rep_theta"))


summary(Coverage.values$Rep_theta)

p19 <-  ggplot(Coverage.values2[Rep_ev_freq == 0.5], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(alpha = 0.7, position=position_dodge(width=0.5))  + ggtitle("Coverage across all theta") + 
  xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
  facet_grid(Rep_theta ~ as.factor(Rep_tau.sq)) + coord_cartesian(ylim = c(0.7, 1))  #+theme_bw()+  scale_colour_grey()

p19



##### Selected REML vs DL coverage example -----
#### USES STEP BIAS
An.Cond <- LogOR.Sim.Results[Rep_theta == theta[5] & Rep_Subj == 4.7 & Rep_ev_freq == 0.5]

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
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq"))

# Not controlled for theta.

p19 <-  ggplot(Coverage.values2[variable %in% c("HC DL", "HC REML")], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(size = 2)  + ggtitle("Coverage for step publication bias, ER 0.5, empirical studies, theta log(4)") + 
  xlab("Number of Studies")  + ylab("Coverage")  + 
  geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  facet_grid(~ as.factor(Rep_tau.sq)) #+ 
  #+theme_bw()+  scale_colour_grey()

p19


ggsave(paste("OR", "Step", "Cov", "NumStudiesbyTausq", "ER", "0.5", "Empirical", ".png", sep = ""), p19, dpi = 300, 
       device = "png", width = 9, height = 5.67, units = "in")

##### Selected REML vs DL Bias/MSE example ----
#### USES OUTCOME BIAS


An.Cond <- LogOR.Sim.Results[Rep_theta == theta[1] & Rep_Subj == 4.7 & Rep_ev_freq == 0.5 & Rep_tau.sq == tau.sq[4]]

Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta, REML = REML_Estimate - Rep_theta,
                                 DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                             by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq"))



p11 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias in outcome bias, ER 0.5, empirical, high tau2, theta log(1/4)") + xlab("Number of studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  )#+
  #facet_grid(~ Rep_theta)

p11

ggsave(paste("OR", "Out", "Bias", "NumStudies", "ER", "0.5", "Empirical", "tau2", tau.sq[4], "theta", theta[1], ".png", sep = ""), p11, dpi = 300, 
       device = "png", width = 9, height = 5.67, units = "in")

#### FE vs REML with bias ----
### Uses step log bias

An.Cond <- LogOR.Sim.Results[Rep_theta == theta[2] & Rep_Subj == Subj[1] & Rep_ev_freq == 0.5]

Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                 DL = DL_Estimate - Rep_theta),
                             by = .(Rep_NumStudies, Rep_tau.sq)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq"))

p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(~ Rep_tau.sq) + ggtitle("Bias with step model, ER 0.5, theta log(0.8), empirical")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  )

p16

ggsave(paste("OR", "Step", "BiasFEDLTheta0.8ER0.5", "Subj", 4.2, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")


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
by = .(Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq"))

# Not controlled for theta.

p19 <-  ggplot(Coverage.values2[variable %in% c("KH DL", "FE", "REML")], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(size = 3, position = position_dodge(0.5))  + ggtitle("Coverage with step model, ER 0.5, theta log(0.8), empirical") + 
  xlab("Number of Studies")  + ylab("Coverage")  + 
  geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  facet_grid(~ as.factor(Rep_tau.sq)) +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[1:3])
  #scale_colour_brewer(type = "qual", palette = "Dark2")#+ 
#+theme_bw()+  scale_colour_grey()

p19


ggsave(paste("OR", "Step", "Cov", "FEDLTheta0.8ER0.5", "Subj", 4.2, ".png", sep = ""), p19, dpi = 300, 
       device = "png", width = 8.01, height = 5.67, units = "in")

p20 <-  ggplot(Coverage.values2[variable %in% c("KH DL", "HC DL", "IVHet")], aes(x = as.factor(Rep_NumStudies), y = value, colour = variable)) +
  geom_point(size = 3, position = position_dodge(0.5))  + ggtitle("Coverage with step model, ER 0.5, theta log(0.8), empirical") + 
  xlab("Number of Studies")  + ylab("Coverage")  + 
  geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") +
  coord_cartesian(ylim = c(0, 1)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  facet_grid(~ as.factor(Rep_tau.sq))  +
  scale_colour_manual(values = brewer.pal(8, "Dark2")[c(4,3,5)])
  #scale_colour_brewer(type = "qual", palette = "Dark2")#+ 
#+theme_bw()+  scale_colour_grey()

p20


ggsave(paste("OR", "Step", "Cov", "HCIVHetTheta0.8ER0.5", "Subj", 4.2, ".png", sep = ""), p20, dpi = 300, 
       device = "png", width = 8.01, height = 5.67, units = "in")


##### Moreno versus other bias -----
# Uses Log Step bias

#### FE vs REML with bias ----
### Uses step log bias

An.Cond <- LogOR.Sim.Results[Rep_theta == theta[2] & Rep_Subj == Subj[1] & Rep_ev_freq == 0.5]

Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                 DL = DL_Estimate - Rep_theta, Moreno = Moreno_Estimate - Rep_theta),
                             by = .(Rep_NumStudies, Rep_tau.sq)]


Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies", "Rep_tau.sq"))

p16 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
  #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
  geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
  ggtitle("Bias") + xlab("Number of Studies") + ylab("Bias") + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.75,0.75)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  facet_grid(~ Rep_tau.sq) + ggtitle("Bias with step model, ER 0.5, theta log(0.8), empirical")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        title=element_text(size=14,face="bold"),
        legend.text = element_text(size=12)
  ) +
  scale_fill_brewer(type = "qual", palette = "Dark2") 

p16

ggsave(paste("OR", "Step", "BiasMorenoFEDLTheta0.8ER0.5", "Subj", 4.2, ".png", sep = ""), p16, dpi = 300, device = "png", width = 8.01, 
       height = 7, units = "in")
