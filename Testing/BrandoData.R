## Remove variables
rm(list = ls())


library(foreign)
library(ggplot2)
library(data.table)

### Import Stata brando data
BRANDO <- read.dta("E:/AFP/R Project/AFPproject1/brando.smallset.dta")

summary(BRANDO)

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

BRANDO$year <- as.numeric(substrRight(BRANDO$trial_name, 4))


##### Excluding study size NAs

BRANDO.2 <- BRANDO[is.na(BRANDO$studysize) == FALSE, ]
BRANDO.2 <- data.table(BRANDO.2)

summary(BRANDO.2)

ggplot(BRANDO.2, aes(x = sizedecile, y = blind)) + geom_bar()
ggplot(BRANDO.2, aes(x = log(studysize))) + geom_histogram()
ggplot(BRANDO.2, aes(x = studysize)) + geom_histogram() + scale_x_log10()

decile.av.sizes <- BRANDO.2[, mean(studysize), by = sizedecile][[2]]

asdf <- BRANDO.2[blind != "NA", summary(blind)[2] / (summary(blind)[1] + summary(blind)[2]), by = sizedecile]
ggplot(asdf, aes(x = decile.av.sizes, y = asdf$V1)) + geom_point() + scale_x_log10()

asdf <- BRANDO.2[alloc != "NA", summary(alloc)[2] / (summary(alloc)[1] + summary(alloc)[2]), by = sizedecile]
ggplot(asdf, aes(x = decile.av.sizes, y = asdf$V1)) + geom_point() + scale_x_log10()

asdf <- BRANDO.2[random != "NA", summary(random)[2] / (summary(random)[1] + summary(random)[2]), by = sizedecile]
ggplot(asdf, aes(x = decile.av.sizes, y = asdf$V1)) + geom_point() + scale_x_log10()

###### By quartile

asdf <- BRANDO.2[blind != "NA", summary(blind)[2] / (summary(blind)[1] + summary(blind)[2]), by = sizequart]
ggplot(asdf, aes(x = sizequart, y = V1)) + geom_point()

asdf <- BRANDO.2[alloc != "NA", summary(alloc)[2] / (summary(alloc)[1] + summary(alloc)[2]), by = sizequart]
ggplot(asdf, aes(x = sizequart, y = V1)) + geom_point()

asdf <- BRANDO.2[random != "NA", summary(random)[2] / (summary(random)[1] + summary(random)[2]), by = sizequart]
ggplot(asdf, aes(x = sizequart, y = V1)) + geom_point()


BRANDO.2[, median(studysize), by = sizequart]
summary(BRANDO.2$studysize)

### Logistic regression on study size

model1 <- glm(blind ~ scale(studysize), family = "binomial", data = BRANDO.2)
summary(model1)

model2 <- glm(alloc ~ scale(studysize), family = "binomial", data = BRANDO.2)
summary(model2)

model3 <- glm(random ~ scale(studysize), family = "binomial", data = BRANDO.2)
summary(model3)

### Logistic regression on year

model1 <- glm(blind ~ scale(year), family = "binomial", data = BRANDO.2)
summary(model1)

model2 <- glm(alloc ~ scale(year), family = "binomial", data = BRANDO.2)
summary(model2)

model3 <- glm(random ~ scale(year), family = "binomial", data = BRANDO.2)
summary(model3)


### Select recent studies - can now do this with year variable
toMatch <- seq(1990, 2017, 1)
BRANDO.3 <- BRANDO.2[grep(paste(toMatch, collapse="|"), trial_name)]

summary(BRANDO.3)

asdf <- BRANDO.3[blind != "NA", summary(blind)[2] / (summary(blind)[1] + summary(blind)[2]), by = sizedecile]
ggplot(asdf, aes(x = sizedecile, y = V1)) + geom_point()

asdf <- BRANDO.3[alloc != "NA", summary(alloc)[2] / (summary(alloc)[1] + summary(alloc)[2]), by = sizedecile]
ggplot(asdf, aes(x = sizedecile, y = V1)) + geom_point()

asdf <- BRANDO.3[random != "NA", summary(random)[2] / (summary(random)[1] + summary(random)[2]), by = sizedecile]
ggplot(asdf, aes(x = sizedecile, y = V1)) + geom_point()

#### Bin by sample size

BRANDO.2[blind != "NA", summary(blind)[2] / (summary(blind)[1] + summary(blind)[2]), by = between(studysize, 500, 1000)]
BRANDO.2[alloc != "NA", summary(alloc)[2] / (summary(alloc)[1] + summary(alloc)[2]), by = between(studysize, 500, 1000)]
BRANDO.2[random != "NA", summary(random)[2] / (summary(random)[1] + summary(random)[2]), by = between(studysize, 500, 1000)]

