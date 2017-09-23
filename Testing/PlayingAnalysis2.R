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

Normal.Simulation <- read.csv("NormalSimulation1.csv")
Normal.Simulation <- data.table(Normal.Simulation)
setkey(Normal.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies")



# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 50

# k = number of studies in series
Studies = c(3,5,10,30,50,100)
#Studies = c(3,5,10,30)

# subj = number of subjects in study, changed to first part of list for easier keying
#Subj <- list(as.integer(c(60,60)), as.integer(c(30,40)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))
Subj <- c(60, 30, 250, 4.2)

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level mean - need good sense of range for SMD
theta = c(-0.5, 0, 0.5, 1)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.06904763, 0.2761905, 5.247619)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 10
Chosen.outcomes <- 1
Sd.split <- 0.5

# Size of per unit bias increase
Bias.multiple <- 1/0.9


StartTime <- proc.time()

registerDoParallel(c1)

r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "metafor"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "Normal.Simulation")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "metafor"), 
               .export = c("Studies", "Subj", "True.sd",
                           "theta", "tau.sq", "Normal.Simulation")
) %dopar% {
  # ID different for analysis
  ID = length(tau.sq) * Reps * length(Studies)
  
  Normal.Sim.Results <- data.table(
    Unique_ID = integer(length = ID),
#     Rep_Number = integer(length = ID),
#     Rep_Subj = numeric(length = ID),
#     Rep_theta = numeric(length = ID),
#     Rep_tau.sq = numeric(length = ID),
#     Rep_NumStudies = numeric(length = ID),
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
PEESE_Estimate = numeric(length = ID),
Mult_se = numeric(length = ID)
  )
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (n in Studies){
      
      for (m in 1:Reps){
        
        
#         counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
#                                 (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
#                                 (match(l, tau.sq)-1) * sum(Studies) * Reps +
#                                 (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
#                                 m
#         )
        
        ## Counter without number of studies
        counter <- as.integer((match(i, Subj)-1) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                (match(k, theta)-1) * length(tau.sq) * length(Studies) * Reps +
                                (match(l, tau.sq)-1) * length(Studies) * Reps +
                                (match(n, Studies)-1) * Reps + 
                                m
        )
        
        ### Temporary data.table
        temp.data <- Normal.Simulation[J(m, i, k, l, n)]
        
        #Fixed and random effects
        
        ma.fe <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
       
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
        
        ma.DL <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
        
        # Henmi Copas
        
        ma.hc.DL <- hc(ma.DL)
        
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
        
       # ma.reml.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE)
        ma.DL.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)
        
        ## Doi
        # estimate is equal to fixed effect, as are weights
        doi.var <- sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) )
        
        
        ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a
        ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
        moreno.est <- ma.moren$fit[[5]][1]
        
        ## Stanley v2 (PET-PEESE) simply takes Egger value if test value > 0.05, Moreno value if < 0.05
        
        ma.egger <- regtest(ma.fe , predictor = "sei", model = "lm")
        
        stan.2.est <- ifelse(ma.egger$pval < 0.05, ma.moren$fit[[5]][1], ma.egger$fit[[5]][1])
      
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
        
        ###### Can potentially only take estimate, standard error, and possibly I2
        
        
        Normal.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
                                                FE_Estimate = ma.fe[[1]],
                                                FE_se = ma.fe$se,
                                                REML_Estimate = ma.reml$b,
                                                REML_se = ma.reml$se,
                                                REML_tau2 = ma.reml$tau2,
                                                DL_Estimate = ma.DL[[1]],
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
                                                PEESE_Estimate = stan.2.est,
                                                Mult_se = ma.mult$se
                                                  )]
        
        
        dummy.counter <- dummy.counter + 1
        
      }
    }
  }
  Normal.Sim.Results
}

Normal.Sim.Results <- r[order(Unique_ID)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * Reps * length(Studies)

Normal.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
Normal.Sim.Results$Rep_NumStudies = rep(Studies, times = ID/(Reps*length(Studies)))
Normal.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(tau.sq)))
Normal.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)), times = length(Subj))
Normal.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTaken <- proc.time() - StartTime

write.csv(Normal.Sim.Results, file = "NormalSimAnalysis1.csv")

rm(r)

