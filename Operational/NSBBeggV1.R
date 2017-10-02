### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)

library(compiler)
enableJIT(3)

set.seed(234)

# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 50

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(2)

# theta = population level mean - need good sense of range for SMD
theta = c(-1.5, -0.3, 0, 0.3, 1.5)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.007, 0.133, 2.533)

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
Bias.multiple <- c(0, log(0.9)/(-1.81) * 2, log(0.81)/(-1.81) * 2)

### Unstandardised mean differrence function

UMD <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ControlGroup <- rnorm(Group1Size, -StudyUMD/2, sd)
  TreatmentGroup <- rnorm(Group2Size, StudyUMD/2, sd)
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( (var(ControlGroup) * (Group1Size - 1) + var(TreatmentGroup) * (Group2Size-1))/ (Group1Size + Group2Size -2) * (1/Group1Size + 1/Group2Size))
  return(c(Studymean, Studysd))
}

### UMD function with multiple outcome bias with frac being sd in first level, num.times = number of outcomes simulated
# outputs vectors ordered by p-val

UMD.mult.out <- function(StudySize, Theta, Heterogeneity, Control_Prop, total.sd, frac, num.times){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  z <- normalCopula(param = frac, dim = num.times)
  Z <- rCopula(Group1Size, z)
  ControlGroup <- qnorm(Z, mean = -StudyUMD/2, sd = total.sd)
  y <- normalCopula(param = frac, dim = num.times)
  Y <- rCopula(Group1Size, y)
  TreatmentGroup <- qnorm(Y, mean = StudyUMD/2, sd = total.sd)
  Studymean <- apply(TreatmentGroup,2,mean) - apply(ControlGroup, 2, mean)
  Studysd <- sqrt( (apply(TreatmentGroup, 2, var) * (Group1Size - 1) + apply(TreatmentGroup, 2, var) * (Group2Size-1))/ (Group1Size + Group2Size -2) )
  Begg_p <- pnorm(Studymean/Studysd)
  return(list(Studymean[order(Begg_p)], Studysd[order(Begg_p)]))
}


##### Can't use doRNG on nested loops, see work around in vignette later

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(123)
r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "copula"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
                          "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                          "Bias.multiple", "UMD.mult.out", "Begg_c")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "copula"), 
               .export = c("Studies", "Subj", "True.sd",
                           "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
                           "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                           "Bias.multiple", "UMD.mult.out", "Begg_c")
) %dopar% {
  
  ID = length(tau.sq) * Reps * sum(Studies)
  
  Normal.Simulation <- data.table(
    Unique_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID)
    #     Study_rejectedMeans = list(length = ID),
    #     Study_rejectedSDs = list(length = ID),
    #     Study_Number.of.biases = integer(length = ID)
  )
  
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (n in Studies){
      
      for (o in 1:n){
        
        for (m in 1:Reps){
          
          
          counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * sum(Studies) * Reps +
                                  (match(k, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                  (match(l, tau.sq)-1) * sum(Studies) * Reps +
                                  (sum(Studies[0:(match(n, Studies)-1)]) + o -1) * Reps +
                                  m
          )
          
          #Select sample size
          if (is.integer(i[1]) == TRUE){
            #Study_patientnumber <- i[1]
            Study_patientnumber <- as.integer( (runif(1, sqrt(i[1]), sqrt(i[2])))^2 )
          } else {
            Study_patientnumber <- round(rlnorm(1, meanlog = i[1], sdlog = i[2]) + 4)
          }
          
          ### Implement Begg and Mazumdar publication bias
          repeat{
            
            Study_summary <- UMD(Study_patientnumber, k, l, controlProp, True.sd)
            Study_mean <- Study_summary[1]
            Study_StanDev <- Study_summary[2]
            Study.n <- as.integer(0.5*Study_patientnumber) * 2
            
            # Study mean changed to absolute - not described in paper
            Begg_weight <-exp(
              -Begg_b * (
                (Begg_sided * pnorm(Study_mean/(Study_StanDev))) 
                ^Begg_a  ) * (Study.n ^ Begg_c)
            ) 
            
            if(rbinom(1,1, Begg_weight) == 1 ){break}
            
          }
          
          
          
          Normal.Simulation[dummy.counter, `:=` (Unique_ID = counter, Study_estimate = Study_mean, Study_sd = Study_StanDev,
                                                 Study_n = Study.n)]
          
          dummy.counter <- dummy.counter + 1
          
        }
        
      }
    }
  }
  Normal.Simulation
}

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
Normal.Simulation <- r[order(Unique_ID)]

ID =  length(Subj) * length(theta) * length(tau.sq) * Reps * sum(Studies)
Normal.Simulation$Rep_Number =  rep(1:Reps, times = ID/Reps)
intermediate <- integer()
for (i in Studies){intermediate <- append(intermediate, rep(i, times = i*Reps))}
Normal.Simulation$Rep_NumStudies = rep(intermediate, times = ID/(Reps*sum(Studies)))
Normal.Simulation$Rep_tau.sq = rep(rep(tau.sq, each = Reps * sum(Studies)), times = ID/(Reps*sum(Studies)*length(tau.sq)))
Normal.Simulation$Rep_theta = rep( rep(theta, each = Reps * sum(Studies) * length(tau.sq)), times = length(Subj))

### Create keyable vector for Subj
Subj2 <- c(60, 20, 250, 4.2)
Normal.Simulation$Rep_Subj = rep(Subj2, each = ID / length(Subj))

TimeTaken <- proc.time() - StartTime

stopCluster(c1)

rm(r)

write.csv(Normal.Simulation, file = "NSBBeggV1.csv")

# df.Normal.Simulation <- as.data.frame(Normal.Simulation)
# df.Normal.Simulation$Study_rejectedMeans <- vapply(df.Normal.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
# df.Normal.Simulation$Study_rejectedSDs <- vapply(df.Normal.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
# write.csv(df.Normal.Simulation, file = "NormalSimulation1.csv")

###### Analysis

library(metafor)

Subj <- c(60, 20, 250, 4.2)

#Normal.Simulation <- read.csv("NormalSimulation1.csv")
#Normal.Simulation <- data.table(Normal.Simulation)
setkey(Normal.Simulation, "Rep_Number", "Rep_Subj", "Rep_theta", "Rep_tau.sq", "Rep_NumStudies")


StartTime <- proc.time()
c1 <- makeCluster(num.Cores)
registerDoParallel(c1)

r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table", "metafor"), 
              .export = c("Studies", "Subj", "True.sd", "Reps",
                          "theta", "tau.sq", "Normal.Simulation")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table", "metafor"), 
               .export = c("Studies", "Subj", "True.sd", "Reps",
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
        ma.fe <- tryCatch({
          rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
        },
        error = function(e){
          return(list(b = NA,  se = NA))
        },
        warning = function(w){
          return(list(list(b = NA,  se = NA)))
        }
        )
        
        #ma.fe <- tryCatch(rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE"), error=function(err) NA)
        
        #ma.fe <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
        
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
        
        ma.DL <- tryCatch({
          rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
        },error = function(e){
          return(list(b = NA, tau2 = NA, se = NA, I2 = NA))
        },warning = function(w){
          return(list(list(b = NA, tau2 = NA, se = NA, I2 = NA)))
        })
        
        ma.hc.DL <- tryCatch({
          hc(ma.DL)
        },error = function(e){
          return(list(b = NA, se = NA))
        },warning = function(w){
          return(list(list(b = NA, se = NA)))
        })
        
        # ma.DL <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
        # # Henmi Copas
        # ma.hc.DL <- hc(ma.DL)
        
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
        
        # try({
        #   # ma.reml.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE)
        #   ma.DL.kh <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)
        #   ## Doi
        #   # estimate is equal to fixed effect, as are weights
        #   doi.var <- sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) )
        #   ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a
        #   ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
        #   moreno.est <- ma.moren$fit[[5]][1]
        #   ## Mawdesley
        #   # Mean of weighted residuals closer to H2 than MSE of unweighted residuals
        #   mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd))
        #   sm.mawd.lm <- summary(mawd.lm)
        #   ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
        #   ma.mult <- rma.uni(temp.data$Study_estimate, temp.data$Study_sd * sqrt(phi.est) , method = "FE")
        # },silent = TRUE)
        
        # ma.reml.kh <- tryCatch({rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE)}
        #                        , error = function(e){return(list(se = NA, ci.lb = NA, ci.ub = NA))
        #                        })
        ma.DL.kh <- tryCatch({rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)}
                             , error = function(e){return(list(se = NA, ci.lb = NA, ci.ub = NA))
                             })
        
        ## Doi
        # estimate is equal to fixed effect, as are weights
        doi.var <- tryCatch(sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) ), error=function(err) NA)
        ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a
        
        moreno.est <- tryCatch({ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
        ma.moren$fit[[5]][1]
        }, error=function(err) NA)
        
        ## Mawdesley
        # Mean of weighted residuals closer to H2 than MSE of unweighted residuals
        
        ma.mult <- tryCatch({
          mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd^2))
          sm.mawd.lm <- summary(mawd.lm)
          ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
          rma.uni(temp.data$Study_estimate, temp.data$Study_sd * sqrt(phi.est) , method = "FE")}
          , error = function(e){return(list(se = NA))
          })
        
        
        Normal.Sim.Results[dummy.counter, `:=` (Unique_ID = counter,
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
Normal.Sim.Results$Rep_NumStudies = rep(rep(Studies, each = Reps), times = ID/(Reps*length(Studies)))
Normal.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(tau.sq)))
Normal.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)), times = length(Subj))
Normal.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTakenAn <- proc.time() - StartTime

write.csv(Normal.Sim.Results, file = "NSBBeggV1An.csv")

rm(r)

### Checking values

sum(is.na(Normal.Sim.Results))
Normal.Sim.Results[is.na(Normal.Sim.Results$REML_Est)]

mean(Normal.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 60 & Rep_tau.sq == 2.533]$DL_Est)
mean(Normal.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 60 & Rep_tau.sq == 0.133]$DL_Est)
mean(Normal.Sim.Results[Rep_theta == 0 & Rep_NumStudies == 100 & Rep_Subj == 60 & Rep_tau.sq == 0.007]$DL_Est)