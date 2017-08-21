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
Reps = 1000

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

# Set up within study reporting bias
Begg_a <- 1.5
Begg_b <- 4
Begg_sided <- 2

Tested.outcomes <- 1
Chosen.outcomes <- 1


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
  Study_rejectedMeans = list(length = ID),
  Study_rejectedSDs = list(length = ID)
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
              Study_patientnumber <- round(rlnorm(1, meanlog = 4.2, sdlog = 1.1))
              
              ### Implement Within study multiple outcomes bias - split variance by simulating values
              # for each individual which represent the between person variability, then sample from these
              # with sd to simulate testing variability
              
              
              Study_mu <- rnorm(1, mean = k, sd = sqrt(l))
              
              # Person level 'true' values with half the total sd
              Person_values <- rnorm(Study_patientnumber, mean = Study_mu, sd = (1/sqrt(2))*j)
              
              # Sample multiple outcome measures from same set of patients, using Person_values as mean
              Study_values <- as.vector(replicate(Tested.outcomes, rnorm(Study_patientnumber, mean = Person_values, sd = (1/sqrt(2))*j)))
              
              Study_mean <- numeric(length = Tested.outcomes)
              Study_StanDev <- numeric(length = Tested.outcomes)
              Begg_p <- numeric(length = Tested.outcomes)
              
              for (z in 1:Tested.outcomes) {
                Study_mean[z] <- mean( Study_values[((z-1)*Study_patientnumber + 1): (z*Study_patientnumber)] )
                Study_StanDev[z] <- sd( Study_values[((z-1)*Study_patientnumber + 1): (z*Study_patientnumber)] )
                Begg_p[z] <- pnorm(-Study_mean[z]/(Study_StanDev[z]))
              }
              
              lv <- which.min(Begg_p)
              
              
              
              
              Normal.Simulation[Unique_ID == counter, `:=` (Rep_Number= m, Rep_Subj = i, Rep_sd = j,
                                                            Rep_theta = k, Rep_tau.sq = l, Rep_NumStudies = n,
                                                            Study_ID = o, 
                                                            Study_estimate = Study_mean[lv], 
                                                            Study_sd = Study_StanDev[lv], 
                                                            Study_n = Study_patientnumber,
                                                            Study_rejectedMeans = list(list(Study_mean[-lv])),
                                                            Study_rejectedSDs = list(list(Study_StanDev[-lv]))
              )]
              
              counter <- counter + 1
            }
          }
        }
      }
    }
  }
}

#write.csv(Normal.Simulation, file = "NormalSimulationOutcomeBias1.csv")

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


#### Problems with list assignment
Normal.Simulation[Unique_ID == 1000, Study_rejectedSDs := Study_StanDev[-lv]]
typeof(Normal.Simulation[1000]$Study_rejectedSDs)

DT = data.table(a=1:2, b=list(1:5,1:10))
DT


sapply(DT$b, length)
typeof(DT$b)
typeof(DT[1]$b)
