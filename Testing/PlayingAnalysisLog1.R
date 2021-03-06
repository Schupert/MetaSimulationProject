### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(metafor)

library(compiler)
enableJIT(3)

## Set working directory

LogOR.Simulation <- read.csv("LogORSimulation1.csv")
LogOR.Simulation <- data.table(LogOR.Simulation)
setkey(LogOR.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies", "Rep_ev_freq")



# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 15

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- c(100, 20, 250, 4.7)

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.75), log(1), log(1.5), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.01777778, 0.04, 3.04)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3, 0.5)

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
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.9


StartTime <- proc.time()

registerDoParallel(c1)

r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "metafor"), 
              .export = c("Studies", "Subj", "True.sd", "Reps",
                          "theta", "tau.sq", "EvFreq", "LogOR.Simulation")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "metafor"), 
               .export = c("Studies", "Subj", "True.sd", "Reps",
                           "theta", "tau.sq", "EvFreq", "LogOR.Simulation")
) %dopar% {
  # ID different for analysis
  ID = length(tau.sq) * Reps * length(Studies) * length(EvFreq)
  
  LogOR.Sim.Results <- data.table(
    Unique_ID = integer(length = ID),
#     Rep_Number = integer(length = ID),
#     Rep_Subj = numeric(length = ID),
#     Rep_theta = numeric(length = ID),
#     Rep_tau.sq = numeric(length = ID),
#     Rep_NumStudies = numeric(length = ID),
#     Rep_ev_freq = numeric(length = ID),
    FE_Estimate = numeric(length = ID),
    FE_se = numeric(length = ID),
    REML_Estimate = numeric(length = ID),
    REML_se = numeric(length = ID),
    REML_tau2 = numeric(length = ID),
    DL_Estimate = numeric(length = ID),
    DL_se = numeric(length = ID),
    DL_tau2 = numeric(length = ID),
    DL_I2 = numeric(length = ID),
    HC_Estimate = numeric(length = ID),
    HC_se = numeric(length = ID),
    KH_REML_CIlb = numeric(length = ID),
    KH_REML_CIub = numeric(length = ID),
    KH_REML_se = numeric(length = ID),
    KH_DL_CIlb = numeric(length = ID),
    KH_DL_CIub = numeric(length = ID),
    KH_DL_se = numeric(length = ID),
    Doi_var = numeric(length = ID),
    Moreno_Estimate = numeric(length = ID),
    Mult_se = numeric(length = ID),
    Num_exc = integer(length = ID)
  )
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (j in EvFreq){
      
      for (n in Studies){
        
        for (m in 1:Reps){
          
          
          counter <- as.integer((match(i, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                  (match(k, theta)-1) * length(tau.sq) * length(EvFreq) * length(Studies) * Reps +
                                  (match(l, tau.sq)-1) * length(EvFreq) * length(Studies) * Reps +
                                  (match(j, EvFreq)-1) * length(Studies) * Reps +
                                  (match(n, Studies)-1) * Reps +
                                  m
          )
          
          
          ### Temporary data.table
          temp.data <- LogOR.Simulation[J(m, i, k, l, n, j)]
          
          temp.data <- temp.data[ !((Study_G1O1 == 0.5 & Study_G2O1 == 0.5) | ( ((Study_n/2 + 1) - Study_G1O1) == 0.5 & ((Study_n/2 + 1) - Study_G2O1) == 0.5 ))]
          Excluded.studies <- n - dim(temp.data)[1]
          
          if (dim(temp.data)[1] < 2){
            LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                   FE_Estimate = NA,
                                                   FE_se = NA,
                                                   REML_Estimate = NA,
                                                   REML_se = NA,
                                                   REML_tau2 = NA,
                                                   DL_Estimate = NA,
                                                   DL_se = NA,
                                                   DL_tau2 = NA,
                                                   DL_I2 = NA,
                                                   HC_Estimate = NA,
                                                   HC_se = NA,
                                                   KH_REML_CIlb = NA,
                                                   KH_REML_CIub = NA,
                                                   KH_REML_se = NA,
                                                   KH_DL_CIlb = NA,
                                                   KH_DL_CIub = NA,
                                                   KH_DL_se = NA,
                                                   Doi_var = NA,
                                                   Moreno_Estimate = NA,
                                                   Mult_se = NA,
                                                   Num_exc = Excluded.studies
            )]
            dummy.counter <- dummy.counter + 1
          } else {

          
            #Fixed and random effects
            tryCatch({
            ma.fe <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
            },
            error = function(e){
              return(list(b = paste(i,k,l,j,n,m, sep = ""),  se = NA))
            },
            warning = function(w){
              return(list(list(b = NA,  se = NA)))
            }
            )
            
            ma.reml <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA))
            },
            warning = function(w){
              return(list(list(b = NA, tau2 = NA, se = NA)))
            }
            )
            
            try({

            ma.DL <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")

            
            # Henmi Copas
            

            ma.hc.DL <- hc(ma.DL)
            },silent = TRUE)

            # Knapp Hartung
            
            ma.reml.kh <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE, control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            },
            warning = function(w){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            }
            )
            
            try({
            
            # ma.reml.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE)

            ma.DL.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)

            
            ## Doi
            # estimate is equal to fixed effect, as are weights

            doi.var <- sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) )
            
            
            ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a

            ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
            moreno.est <- ma.moren$fit[[5]][1]

            
            
            ## Mawdesley
            # Mean of weighted residuals closer to H2 than MSE of unweighted residuals
            
            #           ma.fe$H2
            #           ma.reml$H2
            #ma.DL$H2
            mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd))
            sm.mawd.lm <- summary(mawd.lm)
            #mean(sm.mawd.lm$residuals^2)
            # mean(mawd.lm$residuals^2)
            
            ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
            
            ma.mult <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd * sqrt(phi.est) , method = "FE")
            
            },silent = TRUE)
          
          
          LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                  FE_Estimate = ma.fe[[1]][1],
                                                  FE_se = ma.fe$se,
                                                  REML_Estimate = ma.reml$b[1],
                                                  REML_se = ma.reml$se,
                                                  REML_tau2 = ma.reml$tau2,
                                                  DL_Estimate = ma.DL[[1]][1],
                                                  DL_se = ma.DL$se,
                                                  DL_tau2 = ma.DL$tau2,
                                                  DL_I2 = ma.DL$I2,
                                                  HC_Estimate = ma.hc.DL$b,
                                                  HC_se = ma.hc.DL$se,
                                                  KH_REML_CIlb = ma.reml.kh$ci.lb,
                                                  KH_REML_CIub = ma.reml.kh$ci.ub,
                                                  KH_REML_se = ma.reml.kh$se,
                                                  KH_DL_CIlb = ma.DL.kh$ci.lb,
                                                  KH_DL_CIub = ma.DL.kh$ci.ub,
                                                  KH_DL_se = ma.DL.kh$se,
                                                  Doi_var = doi.var,
                                                  Moreno_Estimate = moreno.est,
                                                  Mult_se = ma.mult$se,
                                                 Num_exc = Excluded.studies
          )]
          
          ###### Can potentially only take estimate, standard error, and possibly I2
          
          #         LogOR.Sim.Results[dummy.counter, `:=` (Rep_Number= m, Rep_Subj = i, 
          #                                                 Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
          #                                                 FE_Estimate = ma.fe[[1]],
          #                                                 FE_Est_Low_CI = ma.fe[[5]],
          #                                                 FE_Est_Up_CI = ma.fe[[6]],
          #                                                 REML_Estimate = ma.reml[[1]],
          #                                                 REML_Est_Low_CI = ma.reml[[5]],
          #                                                 REML_Est_Up_CI = ma.reml[[6]],
          #                                                 REML_tau2 = ma.reml[[8]],
          #                                                 REML_I2 = ma.reml[[22]],
          #                                                 DL_Estimate = ma.DL[[1]],
          #                                                 DL_Est_Low_CI = ma.DL[[5]],
          #                                                 DL_Est_Up_CI = ma.DL[[6]],
          #                                                 DL_tau2 = ma.DL[[8]],
          #                                                 DL_I2 = ma.DL[[22]])]
          
#           LogOR.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
#                                                  FE_Estimate = ma.fe[[1]],
#                                                  FE_Est_Low_CI = ma.fe[[5]],
#                                                  FE_Est_Up_CI = ma.fe[[6]],
#                                                  REML_Estimate = ma.reml[[1]],
#                                                  REML_Est_Low_CI = ma.reml[[5]],
#                                                  REML_Est_Up_CI = ma.reml[[6]],
#                                                  REML_tau2 = ma.reml[[8]],
#                                                  REML_I2 = ma.reml[[22]],
#                                                  DL_Estimate = ma.DL[[1]],
#                                                  DL_Est_Low_CI = ma.DL[[5]],
#                                                  DL_Est_Up_CI = ma.DL[[6]],
#                                                  DL_tau2 = ma.DL[[8]],
#                                                  DL_I2 = ma.DL[[22]])]
#           
          
          dummy.counter <- dummy.counter + 1
          
          }
          
          
        }
      }
    }
  }
  LogOR.Sim.Results
}

LogOR.Sim.Results <- r[order(Unique_ID)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * length(Studies)

LogOR.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
LogOR.Sim.Results$Rep_NumStudies = rep(rep(Studies, each = Reps), times = ID/(Reps*length(Studies)))
LogOR.Sim.Results$Rep_ev_freq = rep(rep(EvFreq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(EvFreq)))
LogOR.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)*length(EvFreq)), times = ID/(Reps*length(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))

# ### Create keyable vector for Subj
# Subj2 <- c(100, 30, 250, 4.7)
LogOR.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTaken <- proc.time() - StartTime

write.csv(LogOR.Sim.Results, file = "LogORSimAnalysis1.csv")

sum(is.na(LogOR.Sim.Results))
head(LogOR.Sim.Results[is.na(REML_Estimate)])

# asdf <- LogOR.Simulation[Rep_Number == 8 & Rep_NumStudies == 5 & Rep_ev_freq == 0.1 & Rep_tau.sq == 0.04 & Rep_theta == log(4) & Rep_Subj == 4.7]
# 
# ma.fe <- rma.uni(asdf$Study_estimate, asdf$Study_sd^2 , method = "FE")

rm(r)

tau.sq = c(0, 0.01777778, 0.04, 3.04)

# mean(as.vector(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 60 & Rep_tau.sq == 3.04,.(DL_I2)]))
mean(LogOR.Sim.Results[Rep_theta == log(1) & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 3.04]$DL_I2)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.04]$DL_I2)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.01777778]$DL_I2)

mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 3.04]$DL_tau2)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.04]$DL_tau2)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.01777778]$DL_tau2)

mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 3.04]$DL_Estimate)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.04]$DL_Estimate)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.01777778]$DL_Estimate)

mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 3.04]$DL_se)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.04]$DL_se)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0.01777778]$DL_se)
mean(LogOR.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 100 & Rep_tau.sq == 0]$DL_se)

mean(LogOR.Sim.Results[Rep_theta == 0]$DL_Est)

LogOR.Sim.Results[Rep_theta == 0, mean(DL_Estimate), by = .(Rep_tau.sq, Rep_Subj)]
LogOR.Sim.Results[Rep_theta == 0, .(mean(DL_Estimate), mean(DL_tau2)), by = .(Rep_tau.sq, Rep_Subj)]

#setkey(LogOR.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies", "Rep_ev_freq")

asdf1 <- LogOR.Simulation[J(16, 100, log(0.25), 0, 3, 0.1)]
ma.fe <- rma.uni(asdf1$Study_estimate, asdf1$Study_sd^2 , method = "FE")

dim(LogOR.Simulation[J(9, 4.7, log(1), 0.04, 3, 0.1)])[1] < 2

temp.data <- asdf1[ !((Study_G1O1 == 0.5 & Study_G2O1 == 0.5) | ( ((Study_n/2 + 1) - Study_G1O1) == 0.5 & ((Study_n/2 + 1) - Study_G2O1) == 0.5 ))]
Excluded.studies <- n - dim(temp.data)[1]
