library(data.table)
library(ggplot2)
library(plyr)
library(segmented)
library(MuMIn)
library(lme4)
library(gridExtra)

setwd("C:/Users/Daniel/Desktop/Analysis/Road vs Interior")
library(MuMIn)
RVF <- read.csv(file="RoadsVsForest.csv", header=TRUE)
str(RVF)
RVF$Observer<- as.factor(RVF$Observer)

#loops model selection via AICc/deltaAIC for all species, requires manual model selection by parsimony for deltaAIC < 2,
#saves AIC output to list "model.rank"
Sp_List<-unique(RVF$Sound)
model.rank=list()
top.model=list()
edr=list()
confidence=list()

###Start Loop
for(i in 1:length(Sp_List)){
  data<-RVF[RVF$Sound==Sp_List[i], ]
  data$x <- -data$Distance^2 # transformed distance
  
  ##global model
  global<-glm(Detected~x+ Habitat:x + Wind:x +Humidity:x -1, data=data, family=binomial(link=cloglog))
  global.observer<-glm(Detected~x+ Habitat:x + Wind:x + Humidity:x + Observer:x -1, data=data, family=binomial(link=cloglog))
  wind<-glm(Detected~x+Habitat:x + Wind:x -1, data=data, family=binomial(link=cloglog))
  wind.observer<-glm(Detected~x+Habitat:x + Wind:x + Observer:x -1, data=data, family=binomial(link=cloglog))
  noweather<-glm(Detected~x+Habitat:x -1, data=data, family=binomial(link=cloglog))
  noweather.observer<-glm(Detected~x+Habitat:x + Observer:x -1, data=data, family=binomial(link=cloglog))
  null<-glm(Detected~x -1, data=data, family=binomial(link=cloglog))
  
  
  #Rank models by AIC
  aic<-model.sel(global, global.observer, wind, wind.observer, noweather, noweather.observer, null)
  model.rank[[as.character(paste(Sp_List[i] , sep=""))]]=aic

}


#Display result of AIC weighting for manual model selection
model.rank["BEKI"]
model.rank["PISI"]
model.rank["BOOW"]
model.rank["NSWO"]
model.rank["BBWA"]
model.rank["TEWA"]
model.rank["LEOW"]
model.rank["BLWA"]
model.rank["BAOW"]
model.rank["WETO"]
model.rank["CORA"]
model.rank["RBNU"]
model.rank["LISP"]
model.rank["CCSP"]
model.rank["GGOW"]
model.rank["YERA"]
model.rank["WTSP"]
model.rank["BAWW"]
model.rank["OSFL"]
model.rank["RBGR"]
model.rank["WAVI"]
model.rank["CATO"]
model.rank["OVEN"]
model.rank["DEJU"]
model.rank["BHCO"]
model.rank["11313Hz"]
model.rank["2000Hz"]
model.rank["4000Hz"]
model.rank["2828Hz"]
model.rank["5656Hz"]
model.rank["8000Hz"]
model.rank["1000Hz"]
model.rank["1414Hz"]

#Display top selected model for each species
top.model["BEKI"]
top.model["PISI"]
top.model["BOOW"]
top.model["NSWO"]
top.model["BBWA"]
top.model["TEWA"]
top.model["LEOW"]
top.model["BLWA"]
top.model["BAOW"]
top.model["WETO"]
top.model["CORA"]
top.model["RBNU"]
top.model["LISP"]
top.model["CCSP"]
top.model["GGOW"]
top.model["YERA"]
top.model["WTSP"]
top.model["BAWW"]
top.model["OSFL"]
top.model["RBGR"]
top.model["WAVI"]
top.model["CATO"]
top.model["OVEN"]
top.model["DEJU"]
top.model["BHCO"]
top.model["T11313Hz"]
top.model["T2000Hz"]
top.model["T4000Hz"]
top.model["T2828Hz"]
top.model["T5656Hz"]
top.model["T8000Hz"]
top.model["T1000Hz"]
top.model["T1414Hz"]

######################
######################
######################
#NO WEATHER MODELS
######################
######################
######################
#BBWA - TOP MODEL = NO WEATHER
BBWA<-RVF[RVF$Sound=="BBWA", ]
BBWA$x <--BBWA$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=BBWA, family=binomial("cloglog"))
summary(m)
top.model[["BBWA"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BBWA.ROAD_C"]]<-edrROAD_C
edr[["BBWA.ROAD_D"]]<-edrROAD_D
edr[["BBWA.EDGE_C"]]<-edrEDGE_C
edr[["BBWA.EDGE_D"]]<-edrEDGE_D
edr[["BBWA.FOREST_C"]]<-edrFOREST_C
edr[["BBWA.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BBWA$Habitat, x=BBWA$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BBWA.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_C"]]
confidence[["BBWA.ROAD_D"]]
confidence[["BBWA.EDGE_C"]]
confidence[["BBWA.EDGE_D"]]
confidence[["BBWA.FOREST_C"]]
confidence[["BBWA.FOREST_D"]]


######################
######################
######################
#BEKI - TOP MODEL = NO WEATHER
BEKI<-RVF[RVF$Sound=="BEKI", ]
BEKI$x <--BEKI$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=BEKI, family=binomial("cloglog"))
summary(m)
top.model[["BEKI"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BEKI.ROAD_C"]]<-edrROAD_C
edr[["BEKI.ROAD_D"]]<-edrROAD_D
edr[["BEKI.EDGE_C"]]<-edrEDGE_C
edr[["BEKI.EDGE_D"]]<-edrEDGE_D
edr[["BEKI.FOREST_C"]]<-edrFOREST_C
edr[["BEKI.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BEKI$Habitat, x=BEKI$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BEKI.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_C"]]
confidence[["BEKI.ROAD_D"]]
confidence[["BEKI.EDGE_C"]]
confidence[["BEKI.EDGE_D"]]
confidence[["BEKI.FOREST_C"]]
confidence[["BEKI.FOREST_D"]]


######################
######################
######################
#BLWA - TOP MODEL = NO WEATHER
BLWA<-RVF[RVF$Sound=="BLWA", ]
BLWA$x <--BLWA$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=BLWA, family=binomial("cloglog"))
summary(m)
top.model[["BLWA"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BLWA.ROAD_C"]]<-edrROAD_C
edr[["BLWA.ROAD_D"]]<-edrROAD_D
edr[["BLWA.EDGE_C"]]<-edrEDGE_C
edr[["BLWA.EDGE_D"]]<-edrEDGE_D
edr[["BLWA.FOREST_C"]]<-edrFOREST_C
edr[["BLWA.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BLWA$Habitat, x=BLWA$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BLWA.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_C"]]
confidence[["BLWA.ROAD_D"]]
confidence[["BLWA.EDGE_C"]]
confidence[["BLWA.EDGE_D"]]
confidence[["BLWA.FOREST_C"]]
confidence[["BLWA.FOREST_D"]]


######################
######################
######################
#DEJU - TOP MODEL = NO WEATHER
DEJU<-RVF[RVF$Sound=="DEJU", ]
DEJU$x <--DEJU$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=DEJU, family=binomial("cloglog"))
summary(m)
top.model[["DEJU"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["DEJU.ROAD_C"]]<-edrROAD_C
edr[["DEJU.ROAD_D"]]<-edrROAD_D
edr[["DEJU.EDGE_C"]]<-edrEDGE_C
edr[["DEJU.EDGE_D"]]<-edrEDGE_D
edr[["DEJU.FOREST_C"]]<-edrFOREST_C
edr[["DEJU.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=DEJU$Habitat, x=DEJU$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["DEJU.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_C"]]
confidence[["DEJU.ROAD_D"]]
confidence[["DEJU.EDGE_C"]]
confidence[["DEJU.EDGE_D"]]
confidence[["DEJU.FOREST_C"]]
confidence[["DEJU.FOREST_D"]]


######################
######################
######################
#LISP - TOP MODEL = NO WEATHER
LISP<-RVF[RVF$Sound=="LISP", ]
LISP$x <--LISP$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=LISP, family=binomial("cloglog"))
summary(m)
top.model[["LISP"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["LISP.ROAD_C"]]<-edrROAD_C
edr[["LISP.ROAD_D"]]<-edrROAD_D
edr[["LISP.EDGE_C"]]<-edrEDGE_C
edr[["LISP.EDGE_D"]]<-edrEDGE_D
edr[["LISP.FOREST_C"]]<-edrFOREST_C
edr[["LISP.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=LISP$Habitat, x=LISP$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["LISP.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_C"]]
confidence[["LISP.ROAD_D"]]
confidence[["LISP.EDGE_C"]]
confidence[["LISP.EDGE_D"]]
confidence[["LISP.FOREST_C"]]
confidence[["LISP.FOREST_D"]]


######################
######################
######################
#OSFL - TOP MODEL = NO WEATHER
OSFL<-RVF[RVF$Sound=="OSFL", ]
OSFL$x <--OSFL$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=OSFL, family=binomial("cloglog"))
summary(m)
top.model[["OSFL"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["OSFL.ROAD_C"]]<-edrROAD_C
edr[["OSFL.ROAD_D"]]<-edrROAD_D
edr[["OSFL.EDGE_C"]]<-edrEDGE_C
edr[["OSFL.EDGE_D"]]<-edrEDGE_D
edr[["OSFL.FOREST_C"]]<-edrFOREST_C
edr[["OSFL.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=OSFL$Habitat, x=OSFL$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["OSFL.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_C"]]
confidence[["OSFL.ROAD_D"]]
confidence[["OSFL.EDGE_C"]]
confidence[["OSFL.EDGE_D"]]
confidence[["OSFL.FOREST_C"]]
confidence[["OSFL.FOREST_D"]]


######################
######################
######################
#PISI - TOP MODEL = NO WEATHER
PISI<-RVF[RVF$Sound=="PISI", ]
PISI$x <--PISI$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=PISI, family=binomial("cloglog"))
summary(m)
top.model[["PISI"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["PISI.ROAD_C"]]<-edrROAD_C
edr[["PISI.ROAD_D"]]<-edrROAD_D
edr[["PISI.EDGE_C"]]<-edrEDGE_C
edr[["PISI.EDGE_D"]]<-edrEDGE_D
edr[["PISI.FOREST_C"]]<-edrFOREST_C
edr[["PISI.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=PISI$Habitat, x=PISI$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["PISI.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_C"]]
confidence[["PISI.ROAD_D"]]
confidence[["PISI.EDGE_C"]]
confidence[["PISI.EDGE_D"]]
confidence[["PISI.FOREST_C"]]
confidence[["PISI.FOREST_D"]]


######################
######################
######################
#RBGR - TOP MODEL = NO WEATHER
RBGR<-RVF[RVF$Sound=="RBGR", ]
RBGR$x <--RBGR$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=RBGR, family=binomial("cloglog"))
summary(m)
top.model[["RBGR"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["RBGR.ROAD_C"]]<-edrROAD_C
edr[["RBGR.ROAD_D"]]<-edrROAD_D
edr[["RBGR.EDGE_C"]]<-edrEDGE_C
edr[["RBGR.EDGE_D"]]<-edrEDGE_D
edr[["RBGR.FOREST_C"]]<-edrFOREST_C
edr[["RBGR.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=RBGR$Habitat, x=RBGR$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["RBGR.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_C"]]
confidence[["RBGR.ROAD_D"]]
confidence[["RBGR.EDGE_C"]]
confidence[["RBGR.EDGE_D"]]
confidence[["RBGR.FOREST_C"]]
confidence[["RBGR.FOREST_D"]]


######################
######################
######################
#TEWA - TOP MODEL = NO WEATHER
TEWA<-RVF[RVF$Sound=="TEWA", ]
TEWA$x <--TEWA$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=TEWA, family=binomial("cloglog"))
summary(m)
top.model[["TEWA"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["TEWA.ROAD_C"]]<-edrROAD_C
edr[["TEWA.ROAD_D"]]<-edrROAD_D
edr[["TEWA.EDGE_C"]]<-edrEDGE_C
edr[["TEWA.EDGE_D"]]<-edrEDGE_D
edr[["TEWA.FOREST_C"]]<-edrFOREST_C
edr[["TEWA.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=TEWA$Habitat, x=TEWA$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["TEWA.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_C"]]
confidence[["TEWA.ROAD_D"]]
confidence[["TEWA.EDGE_C"]]
confidence[["TEWA.EDGE_D"]]
confidence[["TEWA.FOREST_C"]]
confidence[["TEWA.FOREST_D"]]


######################
######################
######################
#WAVI - TOP MODEL = NO WEATHER
WAVI<-RVF[RVF$Sound=="WAVI", ]
WAVI$x <--WAVI$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=WAVI, family=binomial("cloglog"))
summary(m)
top.model[["WAVI"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["WAVI.ROAD_C"]]<-edrROAD_C
edr[["WAVI.ROAD_D"]]<-edrROAD_D
edr[["WAVI.EDGE_C"]]<-edrEDGE_C
edr[["WAVI.EDGE_D"]]<-edrEDGE_D
edr[["WAVI.FOREST_C"]]<-edrFOREST_C
edr[["WAVI.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WAVI$Habitat, x=WAVI$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["WAVI.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_C"]]
confidence[["WAVI.ROAD_D"]]
confidence[["WAVI.EDGE_C"]]
confidence[["WAVI.EDGE_D"]]
confidence[["WAVI.FOREST_C"]]
confidence[["WAVI.FOREST_D"]]


######################
######################
######################
#WTSP - TOP MODEL = NO WEATHER
WTSP<-RVF[RVF$Sound=="WTSP", ]
WTSP$x <--WTSP$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=WTSP, family=binomial("cloglog"))
summary(m)
top.model[["WTSP"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["WTSP.ROAD_C"]]<-edrROAD_C
edr[["WTSP.ROAD_D"]]<-edrROAD_D
edr[["WTSP.EDGE_C"]]<-edrEDGE_C
edr[["WTSP.EDGE_D"]]<-edrEDGE_D
edr[["WTSP.FOREST_C"]]<-edrFOREST_C
edr[["WTSP.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WTSP$Habitat, x=WTSP$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["WTSP.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_C"]]
confidence[["WTSP.ROAD_D"]]
confidence[["WTSP.EDGE_C"]]
confidence[["WTSP.EDGE_D"]]
confidence[["WTSP.FOREST_C"]]
confidence[["WTSP.FOREST_D"]]


######################
######################
######################
#YERA - TOP MODEL = NO WEATHER
YERA<-RVF[RVF$Sound=="YERA", ]
YERA$x <--YERA$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=YERA, family=binomial("cloglog"))
summary(m)
top.model[["YERA"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[5]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
#cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
#(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["YERA.ROAD_C"]]<-edrROAD_C
edr[["YERA.ROAD_D"]]<-edrROAD_D
edr[["YERA.EDGE_C"]]<-edrEDGE_C
edr[["YERA.EDGE_D"]]<-edrEDGE_D
edr[["YERA.FOREST_C"]]<-edrFOREST_C
#edr[["YERA.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
#FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=YERA$Habitat, x=YERA$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[5]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  #cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  #(edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  #FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
#FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["YERA.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
#confidence[["YERA.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_C"]]
confidence[["YERA.ROAD_D"]]
confidence[["YERA.EDGE_C"]]
confidence[["YERA.EDGE_D"]]
confidence[["YERA.FOREST_C"]]
#confidence[["YERA.FOREST_D"]]


######################
######################
######################
#T11313Hz - TOP MODEL = NO WEATHER
T11313Hz<-RVF[RVF$Sound=="11313Hz", ]
T11313Hz$x <--T11313Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T11313Hz, family=binomial("cloglog"))
summary(m)
top.model[["T11313Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T11313Hz.ROAD_C"]]<-edrROAD_C
edr[["T11313Hz.ROAD_D"]]<-edrROAD_D
edr[["T11313Hz.EDGE_C"]]<-edrEDGE_C
edr[["T11313Hz.EDGE_D"]]<-edrEDGE_D
edr[["T11313Hz.FOREST_C"]]<-edrFOREST_C
edr[["T11313Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T11313Hz$Habitat, x=T11313Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T11313Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_C"]]
confidence[["T11313Hz.ROAD_D"]]
confidence[["T11313Hz.EDGE_C"]]
confidence[["T11313Hz.EDGE_D"]]
confidence[["T11313Hz.FOREST_C"]]
confidence[["T11313Hz.FOREST_D"]]


######################
######################
######################
#T8000Hz - TOP MODEL = NO WEATHER
T8000Hz<-RVF[RVF$Sound=="8000Hz", ]
T8000Hz$x <--T8000Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T8000Hz, family=binomial("cloglog"))
summary(m)
top.model[["T8000Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T8000Hz.ROAD_C"]]<-edrROAD_C
edr[["T8000Hz.ROAD_D"]]<-edrROAD_D
edr[["T8000Hz.EDGE_C"]]<-edrEDGE_C
edr[["T8000Hz.EDGE_D"]]<-edrEDGE_D
edr[["T8000Hz.FOREST_C"]]<-edrFOREST_C
edr[["T8000Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T8000Hz$Habitat, x=T8000Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T8000Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_C"]]
confidence[["T8000Hz.ROAD_D"]]
confidence[["T8000Hz.EDGE_C"]]
confidence[["T8000Hz.EDGE_D"]]
confidence[["T8000Hz.FOREST_C"]]
confidence[["T8000Hz.FOREST_D"]]


######################
######################
######################
#T4000Hz - TOP MODEL = NO WEATHER
T4000Hz<-RVF[RVF$Sound=="4000Hz", ]
T4000Hz$x <--T4000Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T4000Hz, family=binomial("cloglog"))
summary(m)
top.model[["T4000Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T4000Hz.ROAD_C"]]<-edrROAD_C
edr[["T4000Hz.ROAD_D"]]<-edrROAD_D
edr[["T4000Hz.EDGE_C"]]<-edrEDGE_C
edr[["T4000Hz.EDGE_D"]]<-edrEDGE_D
edr[["T4000Hz.FOREST_C"]]<-edrFOREST_C
edr[["T4000Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T4000Hz$Habitat, x=T4000Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T4000Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_C"]]
confidence[["T4000Hz.ROAD_D"]]
confidence[["T4000Hz.EDGE_C"]]
confidence[["T4000Hz.EDGE_D"]]
confidence[["T4000Hz.FOREST_C"]]
confidence[["T4000Hz.FOREST_D"]]


######################
######################
######################
#T5656Hz - TOP MODEL = NO WEATHER
T5656Hz<-RVF[RVF$Sound=="5656Hz", ]
T5656Hz$x <--T5656Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T5656Hz, family=binomial("cloglog"))
summary(m)
top.model[["T5656Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T5656Hz.ROAD_C"]]<-edrROAD_C
edr[["T5656Hz.ROAD_D"]]<-edrROAD_D
edr[["T5656Hz.EDGE_C"]]<-edrEDGE_C
edr[["T5656Hz.EDGE_D"]]<-edrEDGE_D
edr[["T5656Hz.FOREST_C"]]<-edrFOREST_C
edr[["T5656Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T5656Hz$Habitat, x=T5656Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T5656Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_C"]]
confidence[["T5656Hz.ROAD_D"]]
confidence[["T5656Hz.EDGE_C"]]
confidence[["T5656Hz.EDGE_D"]]
confidence[["T5656Hz.FOREST_C"]]
confidence[["T5656Hz.FOREST_D"]]


######################
######################
######################
#T1000Hz - TOP MODEL = NO WEATHER
T1000Hz<-RVF[RVF$Sound=="1000Hz", ]
T1000Hz$x <--T1000Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T1000Hz, family=binomial("cloglog"))
summary(m)
top.model[["T1000Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T1000Hz.ROAD_C"]]<-edrROAD_C
edr[["T1000Hz.ROAD_D"]]<-edrROAD_D
edr[["T1000Hz.EDGE_C"]]<-edrEDGE_C
edr[["T1000Hz.EDGE_D"]]<-edrEDGE_D
edr[["T1000Hz.FOREST_C"]]<-edrFOREST_C
edr[["T1000Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T1000Hz$Habitat, x=T1000Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T1000Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_C"]]
confidence[["T1000Hz.ROAD_D"]]
confidence[["T1000Hz.EDGE_C"]]
confidence[["T1000Hz.EDGE_D"]]
confidence[["T1000Hz.FOREST_C"]]
confidence[["T1000Hz.FOREST_D"]]


######################
######################
######################
#T1414Hz - TOP MODEL = NO WEATHER
T1414Hz<-RVF[RVF$Sound=="1414Hz", ]
T1414Hz$x <--T1414Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat -1, data=T1414Hz, family=binomial("cloglog"))
summary(m)
top.model[["T1414Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]
cROAD_D <- cf[1]
cEDGE_C <- cf[1]+cf[2]
cEDGE_D <- cf[1]+cf[3]
cFOREST_C <- cf[1]+cf[4]
cFOREST_D <- cf[1]+cf[5]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T1414Hz.ROAD_C"]]<-edrROAD_C
edr[["T1414Hz.ROAD_D"]]<-edrROAD_D
edr[["T1414Hz.EDGE_C"]]<-edrEDGE_C
edr[["T1414Hz.EDGE_D"]]<-edrEDGE_D
edr[["T1414Hz.FOREST_C"]]<-edrFOREST_C
edr[["T1414Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T1414Hz$Habitat, x=T1414Hz$x)
  m_star <- glm(y ~ x + x:Habitat -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]
  cROAD_D_star <- cf_star[1]
  cEDGE_C_star <- cf_star[1]+cf_star[2]
  cEDGE_D_star <- cf_star[1]+cf_star[3]
  cFOREST_C_star <- cf_star[1]+cf_star[4]
  cFOREST_D_star <- cf_star[1]+cf_star[5]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T1414Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95),na.rm=TRUE)
confidence[["T1414Hz.ROAD_C"]]
confidence[["T1414Hz.ROAD_D"]]
confidence[["T1414Hz.EDGE_C"]]
confidence[["T1414Hz.EDGE_D"]]
confidence[["T1414Hz.FOREST_C"]]
confidence[["T1414Hz.FOREST_D"]]

######################
######################
######################
#GLOBAL MODELS
######################
######################
######################
#BAOW - TOP MODEL = GLOBAL
BAOW<-RVF[RVF$Sound=="BAOW", ]
BAOW$x <--BAOW$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=BAOW, family=binomial("cloglog"))
summary(m)
top.model[["BAOW"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BAOW.ROAD_C"]]<-edrROAD_C
edr[["BAOW.ROAD_D"]]<-edrROAD_D
edr[["BAOW.EDGE_C"]]<-edrEDGE_C
edr[["BAOW.EDGE_D"]]<-edrEDGE_D
edr[["BAOW.FOREST_C"]]<-edrFOREST_C
edr[["BAOW.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BAOW$Habitat, x=BAOW$x, Wind=BAOW$Wind, Humidity=BAOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BAOW.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAOW.ROAD_C"]]
confidence[["BAOW.ROAD_D"]]
confidence[["BAOW.EDGE_C"]]
confidence[["BAOW.EDGE_D"]]
confidence[["BAOW.FOREST_C"]]
confidence[["BAOW.FOREST_D"]]


######################
######################
######################
#BAWW - TOP MODEL = GLOBAL
BAWW<-RVF[RVF$Sound=="BAWW", ]
BAWW$x <--BAWW$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=BAWW, family=binomial("cloglog"))
summary(m)
top.model[["BAWW"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BAWW.ROAD_C"]]<-edrROAD_C
edr[["BAWW.ROAD_D"]]<-edrROAD_D
edr[["BAWW.EDGE_C"]]<-edrEDGE_C
edr[["BAWW.EDGE_D"]]<-edrEDGE_D
edr[["BAWW.FOREST_C"]]<-edrFOREST_C
edr[["BAWW.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BAWW$Habitat, x=BAWW$x, Wind=BAWW$Wind, Humidity=BAWW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BAWW.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_C"]]
confidence[["BAWW.ROAD_D"]]
confidence[["BAWW.EDGE_C"]]
confidence[["BAWW.EDGE_D"]]
confidence[["BAWW.FOREST_C"]]
confidence[["BAWW.FOREST_D"]]


######################
######################
######################
#BHCO - TOP MODEL = GLOBAL
BHCO<-RVF[RVF$Sound=="BHCO", ]
BHCO$x <--BHCO$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=BHCO, family=binomial("cloglog"))
summary(m)
top.model[["BHCO"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BHCO.ROAD_C"]]<-edrROAD_C
edr[["BHCO.ROAD_D"]]<-edrROAD_D
edr[["BHCO.EDGE_C"]]<-edrEDGE_C
edr[["BHCO.EDGE_D"]]<-edrEDGE_D
edr[["BHCO.FOREST_C"]]<-edrFOREST_C
edr[["BHCO.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BHCO$Habitat, x=BHCO$x, Wind=BHCO$Wind, Humidity=BHCO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BHCO.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_C"]]
confidence[["BHCO.ROAD_D"]]
confidence[["BHCO.EDGE_C"]]
confidence[["BHCO.EDGE_D"]]
confidence[["BHCO.FOREST_C"]]
confidence[["BHCO.FOREST_D"]]


######################
######################
######################
#BOOW - TOP MODEL = GLOBAL
BOOW<-RVF[RVF$Sound=="BOOW", ]
BOOW$x <--BOOW$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=BOOW, family=binomial("cloglog"))
summary(m)
top.model[["BOOW"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["BOOW.ROAD_C"]]<-edrROAD_C
edr[["BOOW.ROAD_D"]]<-edrROAD_D
edr[["BOOW.EDGE_C"]]<-edrEDGE_C
edr[["BOOW.EDGE_D"]]<-edrEDGE_D
edr[["BOOW.FOREST_C"]]<-edrFOREST_C
edr[["BOOW.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BOOW$Habitat, x=BOOW$x, Wind=BOOW$Wind, Humidity=BOOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["BOOW.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_C"]]
confidence[["BOOW.ROAD_D"]]
confidence[["BOOW.EDGE_C"]]
confidence[["BOOW.EDGE_D"]]
confidence[["BOOW.FOREST_C"]]
confidence[["BOOW.FOREST_D"]]


######################
######################
######################
#CCSP - TOP MODEL = GLOBAL
CCSP<-RVF[RVF$Sound=="CCSP", ]
CCSP$x <--CCSP$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=CCSP, family=binomial("cloglog"))
summary(m)
top.model[["CCSP"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["CCSP.ROAD_C"]]<-edrROAD_C
edr[["CCSP.ROAD_D"]]<-edrROAD_D
edr[["CCSP.EDGE_C"]]<-edrEDGE_C
edr[["CCSP.EDGE_D"]]<-edrEDGE_D
edr[["CCSP.FOREST_C"]]<-edrFOREST_C
edr[["CCSP.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CCSP$Habitat, x=CCSP$x, Wind=CCSP$Wind, Humidity=CCSP$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["CCSP.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_C"]]
confidence[["CCSP.ROAD_D"]]
confidence[["CCSP.EDGE_C"]]
confidence[["CCSP.EDGE_D"]]
confidence[["CCSP.FOREST_C"]]
confidence[["CCSP.FOREST_D"]]


######################
######################
######################
#CORA - TOP MODEL = GLOBAL
CORA<-RVF[RVF$Sound=="CORA", ]
CORA$x <--CORA$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=CORA, family=binomial("cloglog"))
summary(m)
top.model[["CORA"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["CORA.ROAD_C"]]<-edrROAD_C
edr[["CORA.ROAD_D"]]<-edrROAD_D
edr[["CORA.EDGE_C"]]<-edrEDGE_C
edr[["CORA.EDGE_D"]]<-edrEDGE_D
edr[["CORA.FOREST_C"]]<-edrFOREST_C
edr[["CORA.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CORA$Habitat, x=CORA$x, Wind=CORA$Wind, Humidity=CORA$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["CORA.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_C"]]
confidence[["CORA.ROAD_D"]]
confidence[["CORA.EDGE_C"]]
confidence[["CORA.EDGE_D"]]
confidence[["CORA.FOREST_C"]]
confidence[["CORA.FOREST_D"]]


######################
######################
######################
#GGOW - TOP MODEL = GLOBAL
GGOW<-RVF[RVF$Sound=="GGOW", ]
GGOW$x <--GGOW$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=GGOW, family=binomial("cloglog"))
summary(m)
top.model[["GGOW"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["GGOW.ROAD_C"]]<-edrROAD_C
edr[["GGOW.ROAD_D"]]<-edrROAD_D
edr[["GGOW.EDGE_C"]]<-edrEDGE_C
edr[["GGOW.EDGE_D"]]<-edrEDGE_D
edr[["GGOW.FOREST_C"]]<-edrFOREST_C
edr[["GGOW.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=GGOW$Habitat, x=GGOW$x, Wind=GGOW$Wind, Humidity=GGOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["GGOW.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_C"]]
confidence[["GGOW.ROAD_D"]]
confidence[["GGOW.EDGE_C"]]
confidence[["GGOW.EDGE_D"]]
confidence[["GGOW.FOREST_C"]]
confidence[["GGOW.FOREST_D"]]


######################
######################
######################
#LEOW - TOP MODEL = GLOBAL
LEOW<-RVF[RVF$Sound=="LEOW", ]
LEOW$x <--LEOW$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=LEOW, family=binomial("cloglog"))
summary(m)
top.model[["LEOW"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["LEOW.ROAD_C"]]<-edrROAD_C
edr[["LEOW.ROAD_D"]]<-edrROAD_D
edr[["LEOW.EDGE_C"]]<-edrEDGE_C
edr[["LEOW.EDGE_D"]]<-edrEDGE_D
edr[["LEOW.FOREST_C"]]<-edrFOREST_C
edr[["LEOW.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=LEOW$Habitat, x=LEOW$x, Wind=LEOW$Wind, Humidity=LEOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["LEOW.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_C"]]
confidence[["LEOW.ROAD_D"]]
confidence[["LEOW.EDGE_C"]]
confidence[["LEOW.EDGE_D"]]
confidence[["LEOW.FOREST_C"]]
confidence[["LEOW.FOREST_D"]]


######################
######################
######################
#NSWO - TOP MODEL = GLOBAL
NSWO<-RVF[RVF$Sound=="NSWO", ]
NSWO$x <--NSWO$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=NSWO, family=binomial("cloglog"))
summary(m)
top.model[["NSWO"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["NSWO.ROAD_C"]]<-edrROAD_C
edr[["NSWO.ROAD_D"]]<-edrROAD_D
edr[["NSWO.EDGE_C"]]<-edrEDGE_C
edr[["NSWO.EDGE_D"]]<-edrEDGE_D
edr[["NSWO.FOREST_C"]]<-edrFOREST_C
edr[["NSWO.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=NSWO$Habitat, x=NSWO$x, Wind=NSWO$Wind, Humidity=NSWO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["NSWO.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_C"]]
confidence[["NSWO.ROAD_D"]]
confidence[["NSWO.EDGE_C"]]
confidence[["NSWO.EDGE_D"]]
confidence[["NSWO.FOREST_C"]]
confidence[["NSWO.FOREST_D"]]


######################
######################
######################
#OVEN - TOP MODEL = GLOBAL
OVEN<-RVF[RVF$Sound=="OVEN", ]
OVEN$x <--OVEN$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=OVEN, family=binomial("cloglog"))
summary(m)
top.model[["OVEN"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["OVEN.ROAD_C"]]<-edrROAD_C
edr[["OVEN.ROAD_D"]]<-edrROAD_D
edr[["OVEN.EDGE_C"]]<-edrEDGE_C
edr[["OVEN.EDGE_D"]]<-edrEDGE_D
edr[["OVEN.FOREST_C"]]<-edrFOREST_C
edr[["OVEN.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=OVEN$Habitat, x=OVEN$x, Wind=OVEN$Wind, Humidity=OVEN$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["OVEN.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_C"]]
confidence[["OVEN.ROAD_D"]]
confidence[["OVEN.EDGE_C"]]
confidence[["OVEN.EDGE_D"]]
confidence[["OVEN.FOREST_C"]]
confidence[["OVEN.FOREST_D"]]


######################
######################
######################
#RBNU - TOP MODEL = GLOBAL
RBNU<-RVF[RVF$Sound=="RBNU", ]
RBNU$x <--RBNU$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=RBNU, family=binomial("cloglog"))
summary(m)
top.model[["RBNU"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["RBNU.ROAD_C"]]<-edrROAD_C
edr[["RBNU.ROAD_D"]]<-edrROAD_D
edr[["RBNU.EDGE_C"]]<-edrEDGE_C
edr[["RBNU.EDGE_D"]]<-edrEDGE_D
edr[["RBNU.FOREST_C"]]<-edrFOREST_C
edr[["RBNU.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=RBNU$Habitat, x=RBNU$x, Wind=RBNU$Wind, Humidity=RBNU$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["RBNU.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_C"]]
confidence[["RBNU.ROAD_D"]]
confidence[["RBNU.EDGE_C"]]
confidence[["RBNU.EDGE_D"]]
confidence[["RBNU.FOREST_C"]]
confidence[["RBNU.FOREST_D"]]


######################
######################
######################
#CATO - TOP MODEL = GLOBAL
CATO<-RVF[RVF$Sound=="CATO", ]
CATO$x <--CATO$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=CATO, family=binomial("cloglog"))
summary(m)
top.model[["CATO"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[5]+cf[7]+cf[8]
cROAD_D <- cf[1]+cf[7]+cf[8]
cEDGE_C <- cf[1]+cf[2]+cf[7]+cf[8]
cEDGE_D <- cf[1]+cf[3]+cf[7]+cf[8]
cFOREST_C <- cf[1]+cf[4]+cf[7]+cf[8]
#cFOREST_D <- cf[1]+cf[5]+cf[7]+cf[8]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
#(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["CATO.ROAD_C"]]<-edrROAD_C
edr[["CATO.ROAD_D"]]<-edrROAD_D
edr[["CATO.EDGE_C"]]<-edrEDGE_C
edr[["CATO.EDGE_D"]]<-edrEDGE_D
edr[["CATO.FOREST_C"]]<-edrFOREST_C
#edr[["CATO.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
#FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CATO$Habitat, x=CATO$x, Wind=CATO$Wind, Humidity=CATO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[5]+cf_star[7]+cf_star[8]
  cROAD_D_star <- cf_star[1]+cf_star[7]+cf_star[8]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[8]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[8]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[7]+cf_star[8]
  #cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  #(edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  #FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
#FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["CATO.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
#confidence[["CATO.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_C"]]
confidence[["CATO.ROAD_D"]]
confidence[["CATO.EDGE_C"]]
confidence[["CATO.EDGE_D"]]
confidence[["CATO.FOREST_C"]]
#confidence[["CATO.FOREST_D"]]


######################
######################
######################
#WETO - TOP MODEL = GLOBAL
WETO<-RVF[RVF$Sound=="WETO", ]
WETO$x <--WETO$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=WETO, family=binomial("cloglog"))
summary(m)
top.model[["WETO"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[5]+cf[7]+cf[8]
cROAD_D <- cf[1]+cf[7]+cf[8]
cEDGE_C <- cf[1]+cf[2]+cf[7]+cf[8]
cEDGE_D <- cf[1]+cf[3]+cf[7]+cf[8]
cFOREST_C <- cf[1]+cf[4]+cf[7]+cf[8]
#cFOREST_D <- cf[1]+cf[5]+cf[7]+cf[8]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
#(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["WETO.ROAD_C"]]<-edrROAD_C
edr[["WETO.ROAD_D"]]<-edrROAD_D
edr[["WETO.EDGE_C"]]<-edrEDGE_C
edr[["WETO.EDGE_D"]]<-edrEDGE_D
edr[["WETO.FOREST_C"]]<-edrFOREST_C
#edr[["WETO.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
#FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WETO$Habitat, x=WETO$x, Wind=WETO$Wind, Humidity=WETO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[5]+cf_star[7]+cf_star[8]
  cROAD_D_star <- cf_star[1]+cf_star[7]+cf_star[8]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[8]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[8]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[7]+cf_star[8]
  #cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  #(edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  #FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
#FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["WETO.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
#confidence[["WETO.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_C"]]
confidence[["WETO.ROAD_D"]]
confidence[["WETO.EDGE_C"]]
confidence[["WETO.EDGE_D"]]
confidence[["WETO.FOREST_C"]]
#confidence[["WETO.FOREST_D"]]


######################
######################
######################
#T2000Hz - TOP MODEL = GLOBAL
T2000Hz<-RVF[RVF$Sound=="2000Hz", ]
T2000Hz$x <--T2000Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=T2000Hz, family=binomial("cloglog"))
summary(m)
top.model[["T2000Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T2000Hz.ROAD_C"]]<-edrROAD_C
edr[["T2000Hz.ROAD_D"]]<-edrROAD_D
edr[["T2000Hz.EDGE_C"]]<-edrEDGE_C
edr[["T2000Hz.EDGE_D"]]<-edrEDGE_D
edr[["T2000Hz.FOREST_C"]]<-edrFOREST_C
edr[["T2000Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T2000Hz$Habitat, x=T2000Hz$x, Wind=T2000Hz$Wind, Humidity=T2000Hz$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T2000Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_C"]]
confidence[["T2000Hz.ROAD_D"]]
confidence[["T2000Hz.EDGE_C"]]
confidence[["T2000Hz.EDGE_D"]]
confidence[["T2000Hz.FOREST_C"]]
confidence[["T2000Hz.FOREST_D"]]


######################
######################
######################
#T2828Hz - TOP MODEL = GLOBAL
T2828Hz<-RVF[RVF$Sound=="2828Hz", ]
T2828Hz$x <--T2828Hz$Distance^2
m <- glm(Detected ~ x + x:Habitat + x:Wind + x:Humidity -1, data=T2828Hz, family=binomial("cloglog"))
summary(m)
top.model[["T2828Hz"]]<-m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_C <- cf[1]+cf[6]+cf[8]+cf[9]
cROAD_D <- cf[1]+cf[8]+cf[9]
cEDGE_C <- cf[1]+cf[2]+cf[8]+cf[9]
cEDGE_D <- cf[1]+cf[3]+cf[8]+cf[9]
cFOREST_C <- cf[1]+cf[4]+cf[8]+cf[9]
cFOREST_D <- cf[1]+cf[5]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_C <- sqrt(1/cROAD_C))
(edrROAD_D <- sqrt(1/cROAD_D))
(edrEDGE_C <- sqrt(1/cEDGE_C))
(edrEDGE_D <- sqrt(1/cEDGE_D))
(edrFOREST_C <- sqrt(1/cFOREST_C))
(edrFOREST_D <- sqrt(1/cFOREST_D))
edr[["T2828Hz.ROAD_C"]]<-edrROAD_C
edr[["T2828Hz.ROAD_D"]]<-edrROAD_D
edr[["T2828Hz.EDGE_C"]]<-edrEDGE_C
edr[["T2828Hz.EDGE_D"]]<-edrEDGE_D
edr[["T2828Hz.FOREST_C"]]<-edrFOREST_C
edr[["T2828Hz.FOREST_D"]]<-edrFOREST_D


##Simulate from the fitted model, bootstrapping from predicted values for n=1000
f <- fitted(m)

ROAD_C <- list()
ROAD_D <- list()
EDGE_C <- list()
EDGE_D <- list()
FOREST_C <- list()
FOREST_D <- list()
for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T2828Hz$Habitat, x=T2828Hz$x, Wind=T2828Hz$Wind, Humidity=T2828Hz$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_C_star <- cf_star[1]+cf_star[6]+cf_star[8]+cf_star[9]
  cROAD_D_star <- cf_star[1]+cf_star[8]+cf_star[9]
  cEDGE_C_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cEDGE_D_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFOREST_C_star <- cf_star[1]+cf_star[4]+cf_star[8]+cf_star[9]
  cFOREST_D_star <- cf_star[1]+cf_star[5]+cf_star[8]+cf_star[9]
  (edrROAD_C_star <- sqrt(1/cROAD_C_star))
  (edrROAD_D_star <- sqrt(1/cROAD_D_star))
  (edrEDGE_C_star <- sqrt(1/cEDGE_C_star))
  (edrEDGE_D_star <- sqrt(1/cEDGE_D_star))
  (edrFOREST_C_star <- sqrt(1/cFOREST_C_star))
  (edrFOREST_D_star <- sqrt(1/cFOREST_D_star))
  ROAD_C[[i]] <- edrROAD_C_star
  ROAD_D[[i]] <- edrROAD_D_star
  EDGE_C[[i]] <- edrEDGE_C_star
  EDGE_D[[i]] <- edrEDGE_D_star
  FOREST_C[[i]] <- edrFOREST_C_star
  FOREST_D[[i]] <- edrFOREST_D_star
}

#Calculating 90% confidence intervals on bootstrapped values
ROAD_C_mat <- do.call(rbind, ROAD_C)
ROAD_D_mat <- do.call(rbind, ROAD_D)
EDGE_C_mat <- do.call(rbind, EDGE_C)
EDGE_D_mat <- do.call(rbind, EDGE_D)
FOREST_C_mat <- do.call(rbind, FOREST_C)
FOREST_D_mat <- do.call(rbind, FOREST_D)
confidence[["T2828Hz.ROAD_C"]]<-apply(ROAD_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_D"]]<-apply(ROAD_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.EDGE_C"]]<-apply(EDGE_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.EDGE_D"]]<-apply(EDGE_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FOREST_C"]]<-apply(FOREST_C_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FOREST_D"]]<-apply(FOREST_D_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_C"]]
confidence[["T2828Hz.ROAD_D"]]
confidence[["T2828Hz.EDGE_C"]]
confidence[["T2828Hz.EDGE_D"]]
confidence[["T2828Hz.FOREST_C"]]
confidence[["T2828Hz.FOREST_D"]]



#Export lists of EDR and 90% confidence intervals as spreadsheets
edr.values <- data.frame(unlist(edr))
write.csv(edr.values, file = "edr_values.csv")

confidence.values <- data.frame(unlist(confidence))
write.csv(confidence.values, file = "confidence.csv")

