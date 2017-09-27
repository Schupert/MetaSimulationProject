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
Bias.multiple <- c(log(0.9)/(-1.81) * 2, log(0.81)/(-1.81) * 2)

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

write.csv(Normal.Simulation, file = "NSBBeggV1.csv")

# df.Normal.Simulation <- as.data.frame(Normal.Simulation)
# df.Normal.Simulation$Study_rejectedMeans <- vapply(df.Normal.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
# df.Normal.Simulation$Study_rejectedSDs <- vapply(df.Normal.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
# write.csv(df.Normal.Simulation, file = "NormalSimulation1.csv")