library(MuMIn)
library(pROC)
library(psych)
library(ggplot2)
library(lme4)
library(mrds)
library(Distance)
library(raster)
library(nlme)
library(lmtest)
library(maptools) #load package for sunrise/sunset calcs
library(rgdal) #load package for utm to lat long conversion
library(lubridate) #load package for date and time manipulation
library(base) #load for subset funxtion (data truncation)
library(ggpubr)
library(mfp)
library(yhat)
library(caret)
library(ICC)
library(dplyr)

setwd("C:/Users/Daniel/Desktop/Analysis/Distance Estimation/CONI/CNN/")
CONI.R <- read.csv(file="CONI-Live_convnetBC7_0.2_0.1_detections_distancecorrected.csv", header=TRUE)

str(CONI.R)
CONI.R1 <- CONI.R[complete.cases(CONI.R[ , 17]),]
CONI.R1 <- CONI.R1[CONI.R1$distance.c <= 501, ]

####housekeeping###
coni2 <- CONI.R1
coni2$dist.c<-coni2$distance.c

#Calculate time since sunset

#Convert utm to lat long
utm <- SpatialPoints(cbind(coni2$UTME, coni2$UTMN),proj4string=CRS("+proj=utm +zone=10"))
latlong <- as.data.frame(spTransform(utm,CRS("+proj=longlat")))

#makesure DateTime is in Posix format
#conv.y$DateTime <- gsub("/", "-", conv.y$DateTime)
coni2$DateTime <- as.POSIXct(coni2$DateTime, tz="Canada/Mountain")

#rename columns as lat and long
colnames(latlong) <- c("long", "lat")
coni2 <- cbind(coni2, latlong) #bind lat long to data frame

#define function for calculating sunset
convunset <- function(long, lat, DateTime)
{
  siteX <- SpatialPoints(matrix(c(long, lat), nrow=1), proj4string=CRS("+proj=longlat +datum=WGS84"))
  dateX <- as.POSIXct(DateTime, tz="Canada/Mountain")
  duskX <- crepuscule(siteX, dateX, solarDep=0, direction="dusk", POSIXct.out=TRUE)
  duskX <- duskX$time #keep sunset time only
  return(duskX)
}

coni2$Sunset <- mapply(convunset, coni2$long, coni2$lat, coni2$DateTime) # returns number of seconds since 1970-01-01
coni2$SunsetR <- as.POSIXlt(coni2$Sunset, tz="Canada/Mountain", origin="1970-01-01") # convert to readable format
coni2$SunsetT <- coni2$DateTime - coni2$SunsetR +6 #calculate time since sunset

#euclidean distance
coni2$dist.euc <- sqrt((coni2$dist.c^2)+(coni2$vertical^2))


####scale/centre variables####
coni2$level.st <- scale(coni2$level, center=TRUE, scale = TRUE)
coni2$score.st <- scale(coni2$score, center=TRUE, scale = TRUE)
coni2$vertical.st <- scale(coni2$vertical, center=TRUE, scale = TRUE)
coni2$Wind.st <- scale(coni2$Wind, center=TRUE, scale = TRUE)
coni2$Temp.st <- scale(coni2$Temp, center=TRUE, scale = TRUE)
coni2$Humidity.st <- scale(coni2$Humidity, center=TRUE, scale = TRUE)
coni2$SunsetT.st <- scale(coni2$SunsetT, center=TRUE, scale = TRUE)

coni2$ID <- as.factor(coni2$ID)
coni2$transect <- as.factor(coni2$transect)

###note ID = singing event, transect = individual (i.e. one individual per transect)
####test random effect structure
transect <- glmer(dist.c ~ level.st + (1|transect), data = coni2, family = gaussian(link=log))
event <- glmer(dist.c ~ level.st + (1|ID), data = coni2, family = gaussian(link=log))
nested <- glmer(dist.c ~ level.st + (1|transect/ID), data = coni2, family = gaussian(link=log))
#best random effect structure -> transect only


####explore polynomials
p1 <- glmer(dist.c ~ level.st + (1|transect), data = coni2, family = gaussian(link=log))
p2 <- glmer(dist.c ~ level.st + I(level.st^2) + (1|transect), data = coni2, family = gaussian(link=log))

###model selection
null.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2), data = coni2, family = gaussian(link=log))
m1.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + score.st, data = coni2, family = gaussian(link=log))
m2.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + vertical.st, data = coni2, family = gaussian(link=log))
m3.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + SunsetT.st, data = coni2, family = gaussian(link=log))
m4.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st, data = coni2, family = gaussian(link=log))
m5.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + score.st + SunsetT.st, data = coni2, family = gaussian(link=log))
m6.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + vertical.st + SunsetT.st, data = coni2, family = gaussian(link=log))
m7.CONIR.full <- glm(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + SunsetT.st, data = coni2, family = gaussian(link=log))
m8.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + (1|transect), data = coni2, family = gaussian(link=log))
m9.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Temp.st + (1|transect), data = coni2, family = gaussian(link=log))
m10.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m11.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Temp.st + (1|transect), data = coni2, family = gaussian(link=log))
m12.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m13.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m14.CONIR.full <- glmer(dist.c ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))

aic.CONIR.full <- model.sel(null.CONIR.full, m1.CONIR.full, m2.CONIR.full, m3.CONIR.full, m4.CONIR.full, m5.CONIR.full,
                            m6.CONIR.full, m7.CONIR.full, m8.CONIR.full, m9.CONIR.full, m10.CONIR.full, m11.CONIR.full,
                            m12.CONIR.full, m13.CONIR.full, m14.CONIR.full)

####best model
summary(m2.CONIR.full)
####CONICNN
#Randomly shuffle the data
coni2<-coni2[sample(nrow(coni2)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(coni2)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
coni2.r2 <- data.frame()
coni2.predicted <- data.frame()
for(i in 1:10){
  #Segement data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- coni2[testIndexes, ]
  trainData <- coni2[-testIndexes, ]
  
  #Fit using training data
  #fitted_model <- glmer(dist.c ~ level.st + I(level.st^2) + vertical.st + (1|transect), data = trainData, family = gaussian(link=log))
  fitted_model <- glm(dist.c ~ level.st + I(level.st^2) + vertical.st, data = trainData, family = gaussian(link=log))
  
  #Use fitted model to predict on test data
  predicted_testData <- cbind(testData, fitted = predict(fitted_model, newdata=testData, type="response", se.fit=FALSE, allow.new.levels=TRUE, re.form=NA))
  
  #Estimate prediction error
  colnames(predicted_testData)[colnames(predicted_testData)=="fitted"] <- "Predicted.dist"
  predicted_testData$Error <- abs(predicted_testData$Predicted.dist - predicted_testData$dist.c)
  predicted_testData$Error2 <- (predicted_testData$Predicted.dist - predicted_testData$dist.c)
  
  coni2.predicted <- rbind(coni2.predicted, predicted_testData)
  
  #Validate model
  validate_model <- lm(dist.c ~ Predicted.dist, data = predicted_testData)
  adj.r2 <- summary(validate_model)$adj.r.squared
  coni2.r2 <- rbind(coni2.r2, adj.r2)
  
}

describe(coni2.r2)

colnames(coni2.predicted)[colnames(coni2.predicted)=="distance"] <- "d1"
colnames(coni2.predicted)[colnames(coni2.predicted)=="dist.c"] <- "distance"
colnames(coni2.predicted)[colnames(coni2.predicted)=="distance.c"] <- "d2"

####test euclidean distance
null.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + (1|transect), data = coni2, family = gaussian(link=log))
m1.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + (1|transect), data = coni2, family = gaussian(link=log))
m2.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + vertical.st + (1|transect), data = coni2, family = gaussian(link=log))
m3.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m4.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + (1|transect), data = coni2, family = gaussian(link=log))
m5.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m6.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + vertical.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m7.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m8.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + (1|transect), data = coni2, family = gaussian(link=log))
m9.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Temp.st + (1|transect), data = coni2, family = gaussian(link=log))
m10.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m11.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Temp.st + (1|transect), data = coni2, family = gaussian(link=log))
m12.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m13.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m14.CONIR.full.euc <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + Wind.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))

aic.CONIR.full.euc <- model.sel(null.CONIR.full.euc, m1.CONIR.full.euc, m2.CONIR.full.euc, m3.CONIR.full.euc, m4.CONIR.full.euc, m5.CONIR.full.euc,
                            m6.CONIR.full.euc, m7.CONIR.full.euc, m8.CONIR.full.euc, m9.CONIR.full.euc, m10.CONIR.full.euc, m11.CONIR.full.euc,
                            m12.CONIR.full.euc, m13.CONIR.full.euc, m14.CONIR.full.euc)

####CONICNN
#Randomly shuffle the data
coni2<-coni2[sample(nrow(coni2)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(coni2)),breaks=10,labels=FALSE)

#Perform 10 fold cross validation
coni2.r2 <- data.frame()
coni2.predicted <- data.frame()
for(i in 1:10){
  #Segement data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- coni2[testIndexes, ]
  trainData <- coni2[-testIndexes, ]
  
  #Fit using training data
  fitted_model <- glmer(dist.euc ~ level.st + I(level.st^2) + score.st + vertical.st + (1|transect), data = trainData, family = gaussian(link=log))
  
  #Use fitted model to predict on test data
  predicted_testData <- cbind(testData, fitted = predict(fitted_model, newdata=testData, type="response", se.fit=FALSE, allow.new.levels=TRUE))
  
  #Estimate prediction error
  colnames(predicted_testData)[colnames(predicted_testData)=="fitted"] <- "Predicted.dist"
  predicted_testData$Error <- abs(predicted_testData$Predicted.dist - predicted_testData$dist.c)
  predicted_testData$Error2 <- (predicted_testData$Predicted.dist - predicted_testData$dist.c)
  
  coni2.predicted <- rbind(coni2.predicted, predicted_testData)
  
  #Validate model
  validate_model <- lm(dist.euc ~ Predicted.dist, data = predicted_testData)
  adj.r2 <- summary(validate_model)$adj.r.squared
  coni2.r2 <- rbind(coni2.r2, adj.r2)
  
}

describe(coni2.r2)

####commonality analysis
apsOut <- aps(coni2, "dist.c", list("level.st", "vertical.st"))
domF <- dominance(apsOut)
comF <- commonality(apsOut)



####BEGIN BOOTSTRAP####
#FULL
CONIR.FULL.BOOT.a <- list()
CONIR.FULL.BOOT.a.CvM <- list()
CONIR.FULL.BOOT.p <- list()
CONIR.FULL.BOOT.p.CvM <- list()
i=1
while(i < 1001){
  boot <- data.frame()
  a <- coni2.predicted
  
  ####RESAMPLE EVENTS WITHOUT REPLACEMENT####
  for(j in 1:130){
    b <- a[sample(nrow(a), size = 1, replace = FALSE), ]
    c <- as.character(b$ID)
    boot <- rbind(boot, b)
    a<-a[a$ID != c, ]
  }
  
  ####BEGIN DISTANCE SAMPLING####
  boot["Sample.Label"] <- NA
  boot$Sample.Label <- boot$transect
  boot["Effort"] <- 1
  boot["Region.Label"] <- "SBT"
  boot["Area"] <- 0.2073764 ####in km2, 13.0072ha, 130072 m2
  
  boot.a <- boot
  colnames(boot.a)[colnames(boot.a)=="DISTANCE"] <- "distance"
  
  boot.p <- boot
  colnames(boot.p)[colnames(boot.p)=="distance"] <- "actual"
  colnames(boot.p)[colnames(boot.p)=="Predicted.dist"] <- "distance"
  
  ds_model.a <- try(ds(boot.a, key = "hn", adjustment = "cos", transect = "point", convert.units=0.01))
  ds_model.p <- try(ds(boot.p, key = "hn", adjustment = "cos", transect = "point", convert.units=0.01))
  
  if (class(ds_model.a) != "try-error"){
    if (class(ds_model.p) != "try-error"){ 
      if(ds_model.a$dht$individuals$D$Estimate < 8){
      
        CONIR.FULL.BOOT.a[[i]] <- as.numeric(ds_model.a$dht$individuals$D$Estimate)
        CONIR.FULL.BOOT.p[[i]] <- as.numeric(ds_model.p$dht$individuals$D$Estimate)
      
        CvM.a <- gof_ds(ds_model.a, plot=FALSE)
        CvM.p <- gof_ds(ds_model.p, plot=FALSE)
      
        CONIR.FULL.BOOT.a.CvM[[i]] <- as.numeric(CvM.a$dsgof$CvM$p)
        CONIR.FULL.BOOT.p.CvM[[i]] <- as.numeric(CvM.p$dsgof$CvM$p)
        i=i+1
      }
    }
  }
  
  print(i)
}

CONIR.FULL.BOOT.a_mat <- do.call(rbind, CONIR.FULL.BOOT.a)
CONIR.FULL.BOOT.p_mat <- do.call(rbind, CONIR.FULL.BOOT.p)
CONIR.FULL.BOOT.a.CvM_mat <- do.call(rbind, CONIR.FULL.BOOT.a.CvM)
CONIR.FULL.BOOT.p.CvM_mat <- do.call(rbind, CONIR.FULL.BOOT.p.CvM)

d.CONIR.FULL.BOOT.a <- apply(CONIR.FULL.BOOT.a_mat, 2, quantile, c(0.025, 0.975))
d.CONIR.FULL.BOOT.p <- apply(CONIR.FULL.BOOT.p_mat, 2, quantile, c(0.025, 0.975))

sum(CONIR.FULL.BOOT.a.CvM_mat >= 0.05)
sum(CONIR.FULL.BOOT.p.CvM_mat >= 0.05)









####prelim distance sampling####
CONIR.FULL <- coni2.predicted


CONIR.FULL["Sample.Label"] <- NA
CONIR.FULL$Sample.Label <- CONIR.FULL$transect
CONIR.FULL["Effort"] <- 1
CONIR.FULL["Region.Label"] <- "SBT"
CONIR.FULL["Area"] <- 0.207376 ####in km2, 13.0072ha, 130072 m2

CONIR.full.o_ds <- CONIR.FULL
colnames(CONIR.full.o_ds)[colnames(CONIR.full.o_ds)=="distance"] <- "distance"

CONIR.full.p_ds <- CONIR.FULL
colnames(CONIR.full.p_ds)[colnames(CONIR.full.p_ds)=="distance"] <- "actual"
colnames(CONIR.full.p_ds)[colnames(CONIR.full.p_ds)=="Predicted.dist"] <- "distance"

####ds actual distance
CONIR.full.halfnorm1.o_ds <- ds(CONIR.full.o_ds, key="hn", adjustment="cos", transect="point", convert.units=0.001)
CONIR.full.halfnorm2.o_ds <- ds(CONIR.full.o_ds, key="hn", adjustment="herm", transect="point", convert.units=0.001)
CONIR.full.unifcos.o_ds <- ds(CONIR.full.o_ds, key="unif", adjustment="cos", mono="strict", transect="point", convert.units=0.001)
CONIR.full.hazard1.o_ds <- ds(CONIR.full.o_ds, key="hr", adjustment="cos", transect="point", convert.units=0.001)
CONIR.full.hazard2.o_ds <- ds(CONIR.full.o_ds, key="hr", adjustment="poly", transect="point", convert.units=0.001)
###best model: halfnormal cosine adjustment

####ds predicted distance
CONIR.full.halfnorm1.p_ds <- ds(CONIR.full.p_ds, key="hn", adjustment="cos", transect="point", convert.units=0.001)
CONIR.full.halfnorm2.p_ds <- ds(CONIR.full.p_ds, key="hn", adjustment="herm", transect="point", convert.units=0.001)
CONIR.full.unifcos.p_ds <- ds(CONIR.full.p_ds, key="unif", adjustment="cos", mono="strict", transect="point", convert.units=0.001)
CONIR.full.hazard1.p_ds <- ds(CONIR.full.p_ds, key="hr", adjustment="cos", transect="point", convert.units=0.001)
CONIR.full.hazard2.p_ds <- ds(CONIR.full.p_ds, key="hr", adjustment="poly", transect="point", convert.units=0.001)
###best model: halfnormal cosine adjustment


####test log(DISTANCE) models
null.CONIR.full <- glmer(dist.c ~ level.st + (1|transect), data = coni2, family = gaussian(link=log))
m1.CONIR.full <- glmer(dist.c ~ level.st + score.st + (1|transect), data = coni2, family = gaussian(link=log))
m2.CONIR.full <- glmer(dist.c ~ level.st + vertical.st + (1|transect), data = coni2, family = gaussian(link=log))
m3.CONIR.full <- glmer(dist.c ~ level.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m4.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + (1|transect), data = coni2, family = gaussian(link=log))
m5.CONIR.full <- glmer(dist.c ~ level.st + score.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m6.CONIR.full <- glmer(dist.c ~ level.st + vertical.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m7.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + SunsetT.st + (1|transect), data = coni2, family = gaussian(link=log))
m8.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Wind.st + (1|transect), data = coni2, family = gaussian(link=log))
m9.CONIR.full <- lmer(log(dist.c) ~ level.st + score.st + vertical.st + Temp.st + (1|transect), data = coni2)
m10.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m11.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Wind.st + Temp.st + (1|transect), data = coni2, family = gaussian(link=log))
m12.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Wind.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m13.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))
m14.CONIR.full <- glmer(dist.c ~ level.st + score.st + vertical.st + Wind.st + Temp.st + Humidity.st + (1|transect), data = coni2, family = gaussian(link=log))

aic.CONIR.full <- model.sel(null.CONIR.full, m1.CONIR.full, m2.CONIR.full, m3.CONIR.full, m4.CONIR.full, m5.CONIR.full,
                            m6.CONIR.full, m7.CONIR.full, m8.CONIR.full, m9.CONIR.full, m10.CONIR.full, m11.CONIR.full,
                            m12.CONIR.full, m13.CONIR.full, m14.CONIR.full)

# > summary(m9.CONIR.full)
# Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
# Family: gaussian  ( log )
# Formula: dist.c ~ level.st + score.st + vertical.st + Temp.st + (1 | transect)
# Data: coni2
# 
# AIC      BIC   logLik deviance df.resid 
# 13370.4  13406.3  -6678.2  13356.4     1231 
# 
# Scaled residuals: 
#   Min      1Q  Median      3Q     Max 
# -4.4312 -0.6161 -0.1614  0.3790  4.9127 
# 
# Random effects:
#   Groups   Name        Variance Std.Dev.
# transect (Intercept)  329.8   18.16   
# Residual             2591.9   50.91   
# Number of obs: 1238, groups:  transect, 8
# 
# Fixed effects:
#   Estimate Std. Error t value Pr(>|z|)    
# (Intercept)  4.659404   0.126779   36.75  < 2e-16 ***
#   level.st    -0.658952   0.011423  -57.68  < 2e-16 ***
#   score.st     0.019898   0.008064    2.47   0.0136 *  
#   vertical.st  0.044219   0.010660    4.15 3.35e-05 ***
#   Temp.st      0.085930   0.045070    1.91   0.0566 .  
# ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Correlation of Fixed Effects:
#   (Intr) lvl.st scr.st vrtcl.
# level.st     0.078                     
# score.st    -0.004 -0.451              
# vertical.st -0.033 -0.160  0.048       
# Temp.st      0.049 -0.064  0.013 -0.197

####predict distances####
coni2.predicted <- cbind(coni2, fitted = predict(m9.CONIR.full, newdata=coni2, type="response", se.fit=FALSE, allow.new.levels=TRUE))
colnames(coni2.predicted)[colnames(coni2.predicted)=="fitted"] <- "Predicted.dist"
coni2.predicted$Predicted.dist1 <- exp(coni2.predicted$Predicted.dist)
coni2.predicted$Error <- abs(coni2.predicted$dist.c - coni2.predicted$Predicted.dist1)
coni2.predicted$Error2 <- (coni2.predicted$dist.c - coni2.predicted$Predicted.dist1)

describe(coni2.predicted$Error) #descriptive stats for distance error
describe(coni2.predicted$Error2)

coni2.predicted$AbsErrorPercentage <- (abs(coni2.predicted$Predicted.dist - coni2.predicted$dist.c)/coni2.predicted$dist.c)
coni2.predicted$ErrorPercentage <- ((coni2.predicted$Predicted.dist - coni2.predicted$dist.c)/coni2.predicted$dist.c)
describe(coni2.predicted$AbsErrorPercentage)
describe(coni2.predicted$ErrorPercentage)

colnames(coni2.predicted)[colnames(coni2.predicted)=="distance"] <- "d1"
colnames(coni2.predicted)[colnames(coni2.predicted)=="dist.c"] <- "distance"
colnames(coni2.predicted)[colnames(coni2.predicted)=="distance.c"] <- "d2"

####bootstrap time
conir.boot <- coni2
sample.error <- data.frame()
for(j in 1:1000){
  s <- conir.boot[sample(nrow(conir.boot), 
                         size = floor(runif(1, 10, nrow(conir.boot))), 
                         replace = FALSE), ]
  s <- s[sample(nrow(s)),]
  folds <- cut(seq(1,nrow(s)),breaks=10,labels=FALSE)
  s.r2 <- data.frame()
  s.predicted <- data.frame()
  for(i in 1:10){
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- s[testIndexes, ]
    trainData <- s[-testIndexes, ]
    fitted_model <- glm(dist.c ~ level.st + I(level.st^2) + vertical.st, 
                        data = trainData, family = gaussian(link=log))
    predicted_testData <- cbind(testData, 
                                fitted = predict(fitted_model, 
                                                 newdata=testData, type="response", 
                                                 se.fit=FALSE, allow.new.levels=TRUE))
    colnames(predicted_testData)[colnames(predicted_testData)=="fitted"] <- "Predicted.dist"
    predicted_testData$Error <- abs(predicted_testData$Predicted.dist - predicted_testData$dist.c)
    predicted_testData$Error2 <- (predicted_testData$Predicted.dist - predicted_testData$dist.c)
    s.predicted <- rbind(s.predicted, predicted_testData)
    validate_model <- lm(dist.c ~ Predicted.dist, data = predicted_testData)
    adj.r2 <- summary(validate_model)$adj.r.squared
    s.r2 <- rbind(s.r2, adj.r2)
  }
  temp.error <- data.frame(samples = nrow(s), abs.error = mean(s.predicted$Error),
                           bias = mean(s.predicted$Error2), fit = mean(s.r2$X0))
  sample.error <- rbind(sample.error, temp.error)
  
}
setwd("C:/Users/Daniel/Desktop/Analysis/Distance Estimation/CONI/CNN/")
sample.error <- read.csv(file = "conir_asym.csv")
sample.error.cv <- data.frame(Truncate = numeric(),
                              CV = numeric(),
                              SD = numeric(),
                              BP = numeric())
m <- bptest(lm(abs.error ~ samples, data = sample.error))
a <- max(sample.error$samples)
b <- cv(sample.error$abs.error)
c <- sd(sample.error$abs.error)
d <- as.numeric(m$p.value)
sample.error.cv <- rbind(sample.error.cv, data.frame(Truncate = a, CV = b, SD = c, BP = d))
while(a > 20){
  a = a - 1
  data <- sample.error[sample.error$samples < a, ]
  m1 <- bptest(lm(abs.error ~ samples, data = data))
  b1 <- cv(data$abs.error)
  c1 <- sd(data$abs.error)
  d1 <- as.numeric(m1$p.value)
  sample.error.cv <- rbind(sample.error.cv, data.frame(Truncate = a, CV = b1, SD = c1, BP = d1))
}
