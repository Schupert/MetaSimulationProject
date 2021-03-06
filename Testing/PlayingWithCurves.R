library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)

display.brewer.all()

#### Final altered BEgg and Mazumdar

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
    
    Begg_weight <-exp(
      - 5  * (a^0.5)* b^-.3 
      )
    
    #Begg_weight <- 0.5 * (1 - a)^17 + 0.5
    
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)

asdf <- data.table(dummy1, dummy2, dummy3)

# ggplotColours <- function(n=6, h=c(0, 360) +15){
#   if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
#   hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
# }

brewerColours <- brewer.pal(6, "Dark2")

a <- ggplot(asdf, aes(x = dummy1, y = dummy2, colour = factor(dummy3))) + geom_line(size = 3) + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05) + theme(axis.text=element_text(size=12),
                                                                   axis.title=element_text(size=14,face="bold"),
                                                                   legend.text=element_text(size=12),
                                                                   legend.title=element_text(size=14)) +
  scale_colour_manual(values = brewerColours) + 
  guides(colour = guide_legend(reverse=T))

b <- ggplot(asdf[dummy3 %in% c(200) ], aes(x = dummy1, y = dummy2, colour = factor(dummy3))) + geom_line(size = 3) + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05) + theme(axis.text=element_text(size=12),
                                                                     axis.title=element_text(size=14,face="bold"),
                                                                     legend.text=element_text(size=12),
                                                                     legend.title=element_text(size=14))  +
  scale_colour_manual(values = brewerColours[4]) + guides(colour = guide_legend(reverse=T))

c <- ggplot(asdf[dummy3 %in% c(1000,600, 200) ], aes(x = dummy1, y = dummy2, colour = factor(dummy3))) + geom_line(size = 3) + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05) + theme(axis.text=element_text(size=12),
                                                                     axis.title=element_text(size=14,face="bold"),
                                                                     legend.text=element_text(size=12),
                                                                     legend.title=element_text(size=14))+
  scale_colour_manual(values = brewerColours[4:6]) + guides(colour = guide_legend(reverse=T))




ggsave("AlteredBeggModelPublicationFull.png", a, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
ggsave("AlteredBeggModelPublication0.png", b, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")
ggsave("AlteredBeggModelPublication1.png", c, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


#### Step bias

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    
    Begg_weight <- ifelse( a >= 0.2, 0.25, ifelse(a >=0.05, 0.75, 1))
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)

asdf <- data.frame(dummy1, dummy2, dummy3)

a <- ggplot(asdf, aes(x = dummy1, y = dummy2)) + geom_line(size = 3) + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05) + theme(axis.text=element_text(size=12),
                                                                     axis.title=element_text(size=14,face="bold"),
                                                                     legend.text=element_text(size=12),
                                                                     legend.title=element_text(size=14))

ggsave("HedgesAndVeveaStepModelPublication.png", a, dpi = 300, device = "png", width = 8.01, height = 5.67, units = "in")


#######




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

asdf <- data.frame(dummy1, dummy2, dummy3)

ggplot(asdf, aes(x = dummy1, y = dummy2, colour = factor(dummy3))) + geom_line() + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05)

dummy1 <- numeric()
dummy2 <- numeric()
dummy3 <- numeric()
for (a in seq(0,1, 0.01) ){
  for (b in c(1000,600, 200, 100, 50 , 30)){
    #for (b in c(1, 0.6, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01)){  
    
    #Begg_weight <- 0.6 *  (1 - (1/ (b^0.5) * (a )  ))   + 0.4 / (1 +  (a / 0.07)^10)
    Begg_weight <- 0.6 *(1 - (1/ ( (b/30) +1)) * (a^0.4))  +  0.4 / ( 1 + exp(- 100 * (0.06-a) )) ^ (1/1)
    
    dummy1 <- append(dummy1, a)
    dummy2 <- append(dummy2, Begg_weight)
    dummy3 <- append(dummy3, b)
    
  }
}

plot(dummy1, dummy2, ylim = c(0,1))
abline(v = 0.05)
abline(h = 0.95)

asdf <- data.frame(dummy1, dummy2, dummy3)

ggplot(asdf, aes(x = dummy1, y = dummy2, colour = factor(dummy3))) + geom_line() + xlim(0,1) + ylim(0,1) + xlab("p-value") + ylab("Likelihood of publication") + 
  labs(colour = "Sample size") + geom_vline(xintercept=0.05)



### likelihood of getting bias based on sample size

dummy1 <- numeric()
dummy2 <- numeric()

for (a in seq(1,1000, 1)){
  Weight2 <- 1/((a)^0.06)
  
  dummy1 <- append(dummy1, a)
  dummy2 <- append(dummy2, Weight2)
  
}

plot(dummy1, dummy2, ylim = c(0,1))
abline (h = 0.7)

asdf <- data.frame(dummy1, dummy2)

ggplot(asdf, aes(x = dummy1, y = dummy2 )) + geom_line(size = 2) + ylim(0,1) + xlim(0,1000) + ylab("Likelihood of bias") + xlab("Sample size")

mean(dummy2)


### Empirical distributions

qlnorm(c(0.25, 0.5, 0.75), 4.2, 1.1)
qlnorm(c(0.25, 0.5, 0.75), 4.7, 1.2)
