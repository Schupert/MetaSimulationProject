dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
    
    Begg_weight <-exp(
      - 3  * (b ^ -0.3)* (
        (a)) 
      ^.5 ) 
    
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
    
    #Begg_weight <- 0.6 *  (1 - (1/ (b^0.5) * (a )  ))   + 0.4 / (1 +  (a / 0.07)^10)
    #Begg_weight <- 0.5  *  (1 - (1/ (b^0.5)))  + 0.5 *  (1 - (1/ (b^0.5))) / (1 +  (a / 0.07)^3   )
    Begg_weight <- 0.5  *  (1 - (1/ (b^0.5)))  + 0.5  / (1 +  (a / 0.07)^2   )
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
    
    #Begg_weight <- 0.6 *  (1 - (1/ (b^0.5) * (a )  ))   + 0.4 / (1 +  (a / 0.07)^10)
    Begg_weight <- 0.5 *(1 - (1/ ( sqrt(b))  ))  +  0.5 / ( 1 + exp(- 100 * -(a-0.06) )) ^ (1/1)
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)

plot(dummy1, dummy2)
