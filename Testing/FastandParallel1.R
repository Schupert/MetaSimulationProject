### Remove previous variables
rm(list = ls())

#### Libraries
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)

library(compiler)
#enableJIT(3)

set.seed(123)

# Number of cores for parallel
num.Cores <- detectCores() - 1
c1 <- makeCluster(num.Cores)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 100

# k = number of studies in series
Studies = c(3,5,10,30,50,100)
#Studies = c(3,5,10,30)

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(30,40)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = 2

# theta = population level mean - need good sense of range for SMD
theta = c(-0.5, 0, 0.5, 1)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(1,2,3)

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

### Unstandardised mean differrence function
## Altered with pooled standard deviation formula and sample sizes fixed to be equal

UMD <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ControlGroup <- rnorm(Group1Size, -StudyUMD/2, sd)
  TreatmentGroup <- rnorm(Group2Size, StudyUMD/2, sd)
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( (var(ControlGroup) * (Group1Size - 1) + var(TreatmentGroup) * (Group2Size-1))/ (Group1Size + Group2Size -2) )
  return(c(Studymean, Studysd))
}

### UMD function with multiple outcome bias with frac being sd in first level, num.times = number of outcomes simulated
# outputs vectors ordered by p-val
## Altered to have sample sizes equal and pooled standard deviation

UMD.mult.out <- function(StudySize, Theta, Heterogeneity, Control_Prop, total.sd, frac, num.times){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  ControlGroup1 <- rnorm(Group1Size, -StudyUMD/2, sqrt(frac) * total.sd)
  TreatmentGroup1 <- rnorm(Group2Size, mean = StudyUMD/2, sqrt(frac) * total.sd)
  ControlGroupAll <- replicate(num.times, rnorm(Group1Size, ControlGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
  TreatmentGroupAll <- replicate(num.times, rnorm(Group2Size, TreatmentGroup1, sqrt(1-frac) * total.sd), simplify = FALSE)
  Studymean <- sapply(TreatmentGroupAll, mean) - sapply(ControlGroupAll, mean)
  Studysd <- sqrt( ( sapply(ControlGroupAll, var) * (Group1Size - 1) + sapply(TreatmentGroupAll, var)* (Group2Size-1) ) / (Group1Size + Group2Size -2))
  Begg_p <- pnorm(-Studymean/Studysd)
  return(list(Studymean[order(Begg_p)], Studysd[order(Begg_p)]))
}


##### Can't use doRNG on nested loops, see work around in vignette later

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(123)
r <- foreach (i = Subj, .combine=rbind, .packages = c("data.table"), 
              .export = c("Studies", "Subj", "True.sd",
                          "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
                          "Begg_b", "Begg_sided", "Tested.outcomes", "Chosen.outcomes", "Sd.split",
                          "Bias.multiple", "UMD.mult.out")
) %:% foreach (k = theta, .combine=rbind, .packages = c("data.table"), 
           .export = c("Studies", "Subj", "True.sd",
                       "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
                       "Begg_b", "Begg_sided", "Tested.outcomes", "Chosen.outcomes", "Sd.split",
                       "Bias.multiple", "UMD.mult.out")
  ) %dopar% {
  
  ID = length(tau.sq) * Reps * sum(Studies)
  
  Normal.Simulation <- data.table(
    Unique_ID = integer(length = ID),
    # Rep_Number = integer(length = ID),
    # Rep_Subj = list(length = ID),
    # Rep_theta = numeric(length = ID),
    # Rep_tau.sq = numeric(length = ID),
    # Rep_NumStudies = numeric(length = ID),
    # Study_ID = integer(length = ID),
    Study_estimate = numeric(length = ID),
    Study_sd = numeric(length = ID),
    Study_n = integer(length = ID),
    Study_rejectedMeans = list(length = ID),
    Study_rejectedSDs = list(length = ID),
    Study_Number.of.biases = integer(length = ID)
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
            Study_patientnumber <- round(runif(1, i[1], i[2]))
          } else {
            Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1) + 4)
          }
          
          Study_summary <- UMD(Study_patientnumber, k, l, controlProp, True.sd)
          Study_mean <- Study_summary[1]
          Study_StanDev <- Study_summary[2]
          Study.n <- as.integer(0.5*Study_patientnumber) * 2
          
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
Subj2 <- c(60, 30, 250, 4.2)
Normal.Simulation$Rep_Subj = rep(Subj2, each = ID / length(Subj))

TimeTaken <- proc.time() - StartTime

stopCluster(c1)

#write.csv(Normal.Simulation, file = "NormalSimulation1.csv")

df.Normal.Simulation <- as.data.frame(Normal.Simulation)
df.Normal.Simulation$Study_rejectedMeans <- vapply(df.Normal.Simulation$Study_rejectedMeans, paste, collapse = ", ", character(1L))
df.Normal.Simulation$Study_rejectedSDs <- vapply(df.Normal.Simulation$Study_rejectedSDs, paste, collapse = ", ", character(1L))
write.csv(df.Normal.Simulation, file = "NormalSimulation1.csv")
