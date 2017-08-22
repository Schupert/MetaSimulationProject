#### Libraries
library(data.table)
library(metafor)

#### Set working directory - this is a problem with transferring code around computers. 
####    Current system works when opened using R project - sets wd to project wd/Results
#setwd("D:\\Stats\\AFP\\R Code")
if (length(regmatches(getwd(), gregexpr("/Results", getwd()))) == 0) 
{workingDirectory <- paste(getwd(), "/Results", sep = "")}
if (length(regmatches(workingDirectory, gregexpr("/Results", workingDirectory)))) {setwd(workingDirectory)}

#### set seed for reproduceability
set.seed(1234)

#### Declare variables

# Reps = number of repetitions of experiment
Reps = 10

# k = number of studies in series
Studies = c(1)

# subj = number of subjects in study, likely to be distributed
Subj = 100

# sd = study level standard deviation
True.sd = 1

# theta = population level mean
theta = 0

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(2)

# ?need to state I.sq in advance?

# Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)


# ID = total number of data points required, also used as an ID number. WILL NEED UPDATING
ID =  length(Subj) * length(True.sd) * length(theta) * length(tau.sq) * Reps * sum(Studies) 

### Set up data.table to store results

Normal.Simulation <- data.table(
  Unique_ID = c(1:ID),
  Rep_Number = integer(length = ID),
  Rep_Subj = integer(length = ID),
  Rep_sd = numeric(length = ID),
  Rep_theta = numeric(length = ID),
  Rep_tau.sq = numeric(length = ID),
  Rep_NumStudies = numeric(length = ID),
  Study_ID = integer(length = ID),
  Study_estimate = numeric(length = ID),
  Study_sd = numeric(length = ID),
  Study_n = integer(length = ID),
  Study_Number.of.biases = integer(length = ID)
)

### Simulate
# Order of loops is based on calculation for ID above
# If number in each study is being drawn from distribution will need extra coding

StartTime <- proc.time()

counter <- 1

for (i in Subj){
  
  for (j in True.sd){
    
    for (k in theta){
      
      for (l in tau.sq){
        
        for (m in 1:Reps){
          
          for (n in Studies){
            
            for (o in 1:n){
              
              #Statement left in case of varying number of subjects later
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1)+0.5)
              
              ### Publication bias by step function - currently one sided as per Hedges
              
              repeat{
                
                Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
                Study_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = j)
                Study_mean <- mean(Study_values)
                Study_StanDev <- sd(Study_values)
                
                Begg_p <- pnorm(- Study_mean/(Study_StanDev))
                
                Step_weight <- ifelse(Begg_p < Severity.boundary[1], 1, ifelse(Begg_p < Severity.boundary[2], 0.75, 0.25))
                
                if(rbinom(1,1, Step_weight) == 1 ){break}
                
              }
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean, Study_sd = Study_StanDev, 
                                                            Study_n = Study_patientnumber)]
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

#write.csv(Normal.Simulation, file = "NormalSimulationStepBias1.csv")

TimeTaken <- proc.time() - StartTime

hist(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_estimate)
sd(Normal.Simulation$Study_estimate)
mean(Normal.Simulation$Study_sd)
hist(Normal.Simulation$Study_sd)

plot(Normal.Simulation$Study_estimate ~ Normal.Simulation$Study_n)
plot(Normal.Simulation$Study_sd ~ Normal.Simulation$Study_estimate)

testRes <- rma.uni(Study_estimate, Study_sd^2, data=Normal.Simulation, method="FE")
funnel(testRes)