a <- numeric()

for (i in Subj){
  
  for (k in theta){
    
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
            
            a <- append(a, counter)
            
          }
        }
      }
    }
  }
}

identical(a, 1:95040 )
hist(counter)



UMD.new <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(length(Theta), Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  return(matrix(c(StudyUMD, Group1Size, Group2Size, rep(sd, times = length(Theta))), ncol = 4))
}

UMD.part2 <- function(v){
  ControlGroup <- rnorm(v[2], -v[1]/2, v[4])
  TreatmentGroup <- rnorm(v[3], v[1]/2, v[4])
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( var(ControlGroup)/v[2] + var(TreatmentGroup)/v[3] )
  return(c(Studymean, Studysd))
}

UMD.new2 <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(length(Theta), Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  asdf <- matrix(c(StudyUMD, Group1Size, Group2Size), ncol = 3)
  Grp1 <- apply(asdf, 1, function(x) {rnorm(x[2], -x[1]/2, sd)})
  Grp2 <- apply(asdf, 1, function(x) {rnorm(x[3], x[1]/2, sd)})
  Grp1mean <- apply(Grp1, 2, mean)
  Grp2mean <- apply(Grp2, 2, mean)
  Grp1var <- apply(Grp1, 2, var)
  Grp2var <- apply(Grp2, 2, var)
  Studymean <- Grp2mean - Grp1mean
  Studysd <- sqrt( Grp1var/Group1Size + Grp2var/Group2Size )
  return(c(Studymean, Studysd))
}


UMD <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- as.integer(StudySize - Group1Size)
  ControlGroup <- rnorm(Group1Size, -StudyUMD/2, sd)
  TreatmentGroup <- rnorm(Group2Size, StudyUMD/2, sd)
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( var(ControlGroup)/Group1Size + var(TreatmentGroup)/Group2Size )
  return(c(Studymean, Studysd))
}

system.time({
test1 <- UMD.new(rep(10, times = 10000), rep(0,times = 10000), rep(1,times = 10000), 0.5,2)
apply(test1, 1, UMD.part2)
})

system.time({
  for (i in 1:100000){
    UMD(10, 0, 1, 0.5, 2)
  }
}
)

system.time({
  test1 <- UMD.new2(rep(10, times = 100000), rep(0,times = 100000), rep(1,times = 100000), 0.5,2)
})


library(microbenchmark)

mbm = microbenchmark(
  normal = for (i in 1:10000){UMD(10, 0, 1, 0.5, 2)},
  new = UMD.new2(rep(10, times = 10000), rep(0,times = 10000), rep(1,times = 10000), 0.5,2)
  
)

library(ggplot2)
autoplot(mbm)