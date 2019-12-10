setwd("C:/Users/Daniel/Desktop/Analysis/Human vs ARU/GLM")
ARU <- read.csv(file="Human_ARU.csv", header=TRUE)
#reads in CSV file with road and forest comparison data
str(ARU)
library(MuMIn)
library(pROC)
library(psych)
library(ggplot2)
ARU$Habitat<- as.factor(ARU$Habitat)
ARU$Forest<- as.factor(ARU$Forest)
ARU$Habitat <- relevel(ARU$Habitat, ref = "Roadside")
ARU$Forest <- relevel(ARU$Forest, ref = "Closed")
ARU$x <- -ARU$Distance^2

describe(ARU)

####Is Forest type important or Open vs Closed?
F1 <- glm(Detected~x+Habitat, data=ARU, family=binomial)
F2 <- glm(Detected~x+Forest, data=ARU, family=binomial)
F3 <- glm(Detected~x*Habitat, data=ARU, family=binomial)
F4 <- glm(Detected~x*Forest, data=ARU, family=binomial)
F5 <- glm(Detected~x+Forest+Habitat, data=ARU, family=binomial)
F6 <- glm(Detected~x, data=ARU, family=binomial)
pAIC <- model.sel(F1,F2,F3,F4,F5,F6)
pAIC
####Forest type top model (F3)
####Check for correlation between weather variables
cor1 <- cor(ARU$Wind, ARU$Temp)
cor2 <- cor(ARU$Wind, ARU$Humidity)
cor3 <- cor(ARU$Temp, ARU$Humidity)
####r=-0.67 Temp, Humidity

####Model selection####
model.rank=list()
Sp_List<-unique(ARU$Sound)

for(i in 1:length(Sp_List)){
  data<-ARU[ARU$Sound==Sp_List[i], ]
  data$x <- -data$Distance^2 # transformed distance
  #data$Habitat <- relevel(data$Habitat, ref = "Road")
  null <- glm(Detected~x, data=data, family=binomial)
  m1 <- glm(Detected~x + Habitat + Type, data=data, family=binomial)
  m2 <- glm(Detected~x * Habitat + Type, data=data, family=binomial)
  m3 <- glm(Detected~x + Habitat * Type, data=data, family=binomial)
  m4 <- glm(Detected~x * Habitat + Type + Wind, data=data, family=binomial)
  m5 <- glm(Detected~x * Habitat + Type + Wind + Humidity, data=data, family=binomial)
  m6 <- glm(Detected~x * Habitat + Type + Humidity, data=data, family=binomial)
  m7 <- glm(Detected~x * Habitat + Habitat * Type, data=data, family=binomial)
  m8 <- glm(Detected~x * Habitat + Habitat * Type + Wind, data=data, family=binomial)
  m9 <- glm(Detected~x * Habitat + Habitat * Type + Wind + Humidity, data=data, family=binomial)
  m10 <- glm(Detected~x * Habitat + Habitat * Type + Humidity, data=data, family=binomial)
  
  aic<-model.sel(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,null)
  model.rank[[as.character(paste(Sp_List[i] , sep=""))]]=aic
  
}


model.rank["BEKI"]
model.rank["PISI"]
model.rank["BOOW"]
model.rank["NSWO"]
model.rank["BBWA"]
model.rank["TEWA"]
model.rank["LEOW"]
model.rank["BLWA"]
model.rank["BADO"]
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

#output model selection tables
for(i in 1:length(model.rank)){
  model.selection.table <- data.frame(model.rank[i])
  write.csv(model.selection.table, file=paste0(names(model.rank[i]),".csv"))
}

####top models for each species
top.model=list()
ROC=list()
m2list=list("8000Hz","BLWA","BEKI","PISI","CCSP","OVEN","OVEN","WAVI","DEJU","LISP","BHCO","5656Hz","4000Hz","2828Hz")
m4list=list("YERA")
m5list=list("WETO","CATO")
m6list=list("CORA","RBGR","OSFL","RBNU","WTSP","NSWO","2000Hz","1414Hz","1000Hz")
m7list=list("BAWW","TEWA")
m10list=list("BBWA","LEOW","GGOW","BADO","BOOW")
for(i in 1:length(m2list)){
  data<-ARU[ARU$Sound==m2list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Type, data=train.data, family=binomial)
  top.model[[as.character(paste(m2list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m2list[i] , sep=""))]]=ROC.data
}

for(i in 1:length(m4list)){
  data<-ARU[ARU$Sound==m4list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Type + Wind, data=train.data, family=binomial)
  top.model[[as.character(paste(m4list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m4list[i] , sep=""))]]=ROC.data
}

for(i in 1:length(m5list)){
  data<-ARU[ARU$Sound==m5list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.data, family=binomial)
  top.model[[as.character(paste(m5list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m5list[i] , sep=""))]]=ROC.data
}

for(i in 1:length(m6list)){
  data<-ARU[ARU$Sound==m6list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Type + Humidity, data=train.data, family=binomial)
  top.model[[as.character(paste(m6list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m6list[i] , sep=""))]]=ROC.data
}

for(i in 1:length(m7list)){
  data<-ARU[ARU$Sound==m7list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Habitat*Type, data=train.data, family=binomial)
  top.model[[as.character(paste(m7list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m7list[i] , sep=""))]]=ROC.data
}

for(i in 1:length(m10list)){
  data<-ARU[ARU$Sound==m10list[i], ]
  data$x <- -data$Distance^2 # transformed distance
  indexes=sample(1:nrow(data), size=0.3*nrow(data))
  test.data=data[indexes,]
  train.data=data[-indexes,]
  nrow(test.data)    #255 obs
  nrow(train.data)    #595 obs
  
  
  model<-glm(Detected~x*Habitat + Habitat*Type + Humidity, data=train.data, family=binomial)
  top.model[[as.character(paste(m10list[i] , sep=""))]]=model
  
  test.data.predicted<-cbind(test.data, fitted = predict(model, newdata = test.data, type = "response", se.fit = TRUE))
  colnames(test.data.predicted)[colnames(test.data.predicted)=="fitted.fit"] <- "Predicted"
  ROC.data<-roc(test.data.predicted$Detected, test.data.predicted$Predicted) 
  ROC[[as.character(paste(m10list[i] , sep=""))]]=ROC.data
}

model.coef=list()
for(i in 1:length(top.model)){
  estimate.sum <- data.frame("Variable" = variable.names(top.model[[i]]))
  estimate.sum <- cbind(estimate.sum, "Estimate" = coef(top.model[[i]]))
  estimate.sum <- cbind(estimate.sum, "StdError"=coef(summary(top.model[[i]]))[,2])
  rownames(estimate.sum) <- c()
  estimate.sum$Estimate <- signif(estimate.sum$Estimate, digits = 4)
  estimate.sum$StdError <- signif(estimate.sum$StdError, digits = 4)
  write.csv(estimate.sum, paste0("coef_", names(top.model[i]), ".csv"))
}

ARUcoef <- read.csv(file="model_coef_wSTDERROR.csv", header=TRUE)
ARUcoef <- t(ARUcoef)
ARUcoef <- data.frame(ARUcoef)
write.csv(ARUcoef, "SPP_ESTIMATES_WSTDERROR.csv")


####mean betas
mean.beta <- data.frame("Sound" = NA, "Variable" = NA, "Estimate" = NA)
for(i in 1:length(top.model)){
  mean.beta1 <- data.frame("Sound" = names(top.model[i]))
  mean.beta1 <- cbind(mean.beta1, "Variable" = variable.names(top.model[[i]]))
  mean.beta1 <- cbind(mean.beta1, "Estimate" = coef(top.model[[i]]))
  mean.beta <- rbind(mean.beta1, mean.beta)
}
rownames(mean.beta) <- c()

#A1 - distance
A1d<-mean.beta[mean.beta$Variable=="x", ]
A1d <- A1d[complete.cases(A1d),]
mean(A1d[["Estimate"]])
sd(A1d[["Estimate"]])

#A1 - veg
A1conifer <- mean.beta[mean.beta$Variable=="HabitatConiferous",]
A1conifer <- A1conifer[complete.cases(A1conifer),]
mean(A1conifer[["Estimate"]])
sd(A1conifer[["Estimate"]])

A1dec <- mean.beta[mean.beta$Variable=="HabitatDeciduous",]
A1dec <- A1dec[complete.cases(A1dec),]
mean(A1dec[["Estimate"]])
sd(A1dec[["Estimate"]])

#A1 - ARU
A1s2 <- mean.beta[mean.beta$Variable=="TypeSM2",]
A1s2 <- A1s2[complete.cases(A1s2),]
mean(A1s2[["Estimate"]])
sd(A1s2[["Estimate"]])

A1s3 <- mean.beta[mean.beta$Variable=="TypeSM3",]
A1s3 <- A1s3[complete.cases(A1s3),]
mean(A1s3[["Estimate"]])
sd(A1s3[["Estimate"]])

A1RF <- mean.beta[mean.beta$Variable=="TypeRiverForks",]
A1RF <- A1RF[complete.cases(A1RF),]
mean(A1RF[["Estimate"]])
sd(A1RF[["Estimate"]])

A1z <- mean.beta[mean.beta$Variable=="TypeZoom",]
A1z <- A1z[complete.cases(A1z),]
mean(A1z[["Estimate"]])
sd(A1z[["Estimate"]])

#A1 - weather
A1h <- mean.beta[mean.beta$Variable=="Humidity",]
A1h <- A1h[complete.cases(A1h),]
mean(A1h[["Estimate"]])
sd(A1h[["Estimate"]])

A1w <- mean.beta[mean.beta$Variable=="Wind",]
A1w <- A1w[complete.cases(A1w),]
mean(A1w[["Estimate"]])
sd(A1w[["Estimate"]])


####EDR summary stats
EDRsum <- read.csv(file="EDRsummaries.csv", header=TRUE)

mean(EDRsum[["Human"]], na.rm=TRUE)
sd(EDRsum[["Human"]], na.rm=TRUE) 

mean(EDRsum[["SM2"]], na.rm=TRUE)
sd(EDRsum[["SM2"]], na.rm=TRUE) 

mean(EDRsum[["SM3"]], na.rm=TRUE)
sd(EDRsum[["SM3"]], na.rm=TRUE) 

mean(EDRsum[["RiverForks"]], na.rm=TRUE)
sd(EDRsum[["RiverForks"]], na.rm=TRUE) 

mean(EDRsum[["Zoom"]], na.rm=TRUE)
sd(EDRsum[["Zoom"]], na.rm=TRUE)

EDRroad <- EDRsum[EDRsum$Habitat=="Road",]
EDRroad <- EDRroad[complete.cases(EDRroad),]
EDRroad <- EDRroad[-c(1:2)]
EDRroad <- unlist(EDRroad)
mean(EDRroad)
sd(EDRroad)

EDRconi <- EDRsum[EDRsum$Habitat=="Conifer",]
EDRconi <- EDRconi[complete.cases(EDRconi),]
EDRconi <- EDRconi[-c(1:2)]
EDRconi <- unlist(EDRconi)
mean(EDRconi)
sd(EDRconi)

EDRdec <- EDRsum[EDRsum$Habitat=="Deciduous",]
EDRdec <- EDRdec[complete.cases(EDRdec),]
EDRdec <- EDRdec[-c(1:2)]
EDRdec <- unlist(EDRdec)
mean(EDRdec)
sd(EDRdec)

####MDD summary stats
MDDsum <- read.csv(file="MDDsummaries.csv", header=TRUE)

mean(MDDsum[["Human"]], na.rm=TRUE)
sd(MDDsum[["Human"]], na.rm=TRUE) 

mean(MDDsum[["SM2"]], na.rm=TRUE)
sd(MDDsum[["SM2"]], na.rm=TRUE) 

mean(MDDsum[["SM3"]], na.rm=TRUE)
sd(MDDsum[["SM3"]], na.rm=TRUE) 

mean(MDDsum[["RiverForks"]], na.rm=TRUE)
sd(MDDsum[["RiverForks"]], na.rm=TRUE) 

mean(MDDsum[["Zoom"]], na.rm=TRUE)
sd(MDDsum[["Zoom"]], na.rm=TRUE) 

MDDroad <- MDDsum[MDDsum$Habitat=="Road",]
MDDroad <- MDDroad[complete.cases(MDDroad),]
MDDroad <- MDDroad[-c(1:2)]
MDDroad <- unlist(MDDroad)
mean(MDDroad)
sd(MDDroad)

MDDconi <- MDDsum[MDDsum$Habitat=="Conifer",]
MDDconi <- MDDconi[complete.cases(MDDconi),]
MDDconi <- MDDconi[-c(1:2)]
MDDconi <- unlist(MDDconi)
mean(MDDconi)
sd(MDDconi)

MDDdec <- MDDsum[MDDsum$Habitat=="Deciduous",]
MDDdec <- MDDdec[complete.cases(MDDdec),]
MDDdec <- MDDdec[-c(1:2)]
MDDdec <- unlist(MDDdec)
mean(MDDdec)
sd(MDDdec)

####EDR correction stats
EDRcorr <- read.csv(file="EDRcorrect.csv", header=TRUE)

mean(EDRcorr[["Human"]], na.rm=TRUE)
sd(EDRcorr[["Human"]], na.rm=TRUE) 

mean(EDRcorr[["SM2"]], na.rm=TRUE)
sd(EDRcorr[["SM2"]], na.rm=TRUE) 

mean(EDRcorr[["SM3"]], na.rm=TRUE)
sd(EDRcorr[["SM3"]], na.rm=TRUE) 

mean(EDRcorr[["RiverForks"]], na.rm=TRUE)
sd(EDRcorr[["RiverForks"]], na.rm=TRUE) 

mean(EDRcorr[["Zoom"]], na.rm=TRUE)
sd(EDRcorr[["Zoom"]], na.rm=TRUE)

####MDD correction stats
MDDcorr <- read.csv(file="MDDcorrect.csv", header=TRUE)

mean(MDDcorr[["Human"]], na.rm=TRUE)
sd(MDDcorr[["Human"]], na.rm=TRUE) 

mean(MDDcorr[["SM2"]], na.rm=TRUE)
sd(MDDcorr[["SM2"]], na.rm=TRUE) 

mean(MDDcorr[["SM3"]], na.rm=TRUE)
sd(MDDcorr[["SM3"]], na.rm=TRUE) 

mean(MDDcorr[["RiverForks"]], na.rm=TRUE)
sd(MDDcorr[["RiverForks"]], na.rm=TRUE) 

mean(MDDcorr[["Zoom"]], na.rm=TRUE)
sd(MDDcorr[["Zoom"]], na.rm=TRUE)

####Look at interactions
BADO_INT <- ARU[ARU$Sound=="BADO",]
indexes=sample(1:nrow(BADO_INT), size=0.3*nrow(BADO_INT))
test.BADO_INT=BADO_INT[indexes,]
train.BADO_INT=BADO_INT[-indexes,]
nrow(test.BADO_INT)    #255 obs
nrow(train.BADO_INT)    #595 obs


model<-glm(Detected~x*Habitat + Habitat*Type -1, data=train.BADO_INT, family=binomial(link=cloglog))

test.BADO_INT.predicted<-cbind(test.BADO_INT, fitted = predict(model, newdata = test.BADO_INT, type = "response", se.fit = TRUE))
colnames(test.BADO_INT.predicted)[colnames(test.BADO_INT.predicted)=="fitted.fit"] <- "Predicted"
p.tmp<-ggplot(test.BADO_INT.predicted, aes(Distance, Predicted, col=factor(Habitat))) +
  geom_point(data=test.BADO_INT.predicted,aes(x=Distance,y=Predicted), pch = 21) +
  geom_line(aes(group=interaction(Type,Habitat), col=factor(Habitat), lty=factor(Type)))+ 
  ylab("Probability of detecting CCSP") + 
  xlab("Distance from Point Count (m)") + 
  theme(panel.background = element_rect(fill = 'white'))+ 
  theme(axis.line = element_line(colour = "black", size = 1))
p.tmp

BAWW_INT <- ARU[ARU$Sound=="BAWW",]
indexes=sample(1:nrow(BAWW_INT), size=0.3*nrow(BAWW_INT))
test.BAWW_INT=BAWW_INT[indexes,]
train.BAWW_INT=BAWW_INT[-indexes,]
nrow(test.BAWW_INT)    #255 obs
nrow(train.BAWW_INT)    #595 obs


model<-glm(Detected~x*Habitat + Habitat*Type -1, data=train.BAWW_INT, family=binomial(link=cloglog))

test.BAWW_INT.predicted<-cbind(test.BAWW_INT, fitted = predict(model, newdata = test.BAWW_INT, type = "response", se.fit = TRUE))
colnames(test.BAWW_INT.predicted)[colnames(test.BAWW_INT.predicted)=="fitted.fit"] <- "Predicted"
p.tmp<-ggplot(test.BAWW_INT.predicted, aes(Distance, Predicted, col=factor(Habitat))) +
  geom_point(data=test.BAWW_INT.predicted,aes(x=Distance,y=Predicted), pch = 21) +
  geom_line(aes(group=interaction(Type,Habitat), col=factor(Habitat), lty=factor(Type)))+ 
  ylab("Probability of detecting CCSP") + 
  xlab("Distance from Point Count (m)") + 
  theme(panel.background = element_rect(fill = 'white'))+ 
  theme(axis.line = element_line(colour = "black", size = 1))
p.tmp

TEWA_INT <- ARU[ARU$Sound=="TEWA",]
indexes=sample(1:nrow(TEWA_INT), size=0.3*nrow(TEWA_INT))
test.TEWA_INT=TEWA_INT[indexes,]
train.TEWA_INT=TEWA_INT[-indexes,]
nrow(test.TEWA_INT)    #255 obs
nrow(train.TEWA_INT)    #595 obs


model<-glm(Detected~x*Habitat + Habitat*Type -1, data=train.TEWA_INT, family=binomial(link=cloglog))

test.TEWA_INT.predicted<-cbind(test.TEWA_INT, fitted = predict(model, newdata = test.TEWA_INT, type = "response", se.fit = TRUE))
colnames(test.TEWA_INT.predicted)[colnames(test.TEWA_INT.predicted)=="fitted.fit"] <- "Predicted"
p.tmp<-ggplot(test.TEWA_INT.predicted, aes(Distance, Predicted, col=factor(Habitat))) +
  geom_point(data=test.TEWA_INT.predicted,aes(x=Distance,y=Predicted), pch = 21) +
  geom_line(aes(group=interaction(Type,Habitat), col=factor(Habitat), lty=factor(Type)))+ 
  ylab("Probability of detecting CCSP") + 
  xlab("Distance from Point Count (m)") + 
  theme(panel.background = element_rect(fill = 'white'))+ 
  theme(axis.line = element_line(colour = "black", size = 1))
p.tmp

GGOW_INT <- ARU[ARU$Sound=="GGOW",]
indexes=sample(1:nrow(GGOW_INT), size=0.3*nrow(GGOW_INT))
test.GGOW_INT=GGOW_INT[indexes,]
train.GGOW_INT=GGOW_INT[-indexes,]
nrow(test.GGOW_INT)    #255 obs
nrow(train.GGOW_INT)    #595 obs


model<-glm(Detected~x*Habitat + Habitat*Type -1, data=train.GGOW_INT, family=binomial(link=cloglog))

test.GGOW_INT.predicted<-cbind(test.GGOW_INT, fitted = predict(model, newdata = test.GGOW_INT, type = "response", se.fit = TRUE))
colnames(test.GGOW_INT.predicted)[colnames(test.GGOW_INT.predicted)=="fitted.fit"] <- "Predicted"
p.tmp<-ggplot(test.GGOW_INT.predicted, aes(Distance, Predicted, col=factor(Habitat))) +
  geom_point(data=test.GGOW_INT.predicted,aes(x=Distance,y=Predicted), pch = 21) +
  geom_line(aes(group=interaction(Type,Habitat), col=factor(Habitat), lty=factor(Type)))+ 
  ylab("Probability of detecting CCSP") + 
  xlab("Distance from Point Count (m)") + 
  theme(panel.background = element_rect(fill = 'white'))+ 
  theme(axis.line = element_line(colour = "black", size = 1))
p.tmp


#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
####x*Habitat + Type - no weather variables####
####CCSP
CCSP<-ARU[ARU$Sound=="CCSP", ]


indexes=sample(1:nrow(CCSP), size=0.3*nrow(CCSP))
test.CCSP=CCSP[indexes,]
train.CCSP=CCSP[-indexes,]
nrow(test.CCSP)    #255 obs
nrow(train.CCSP)    #595 obs


CCSP.model<-glm(Detected~x*Habitat + Type, data=train.CCSP, family=binomial)
top.model[["CCSP"]]=CCSP.model

test.CCSP.predicted<-cbind(test.CCSP, fitted = predict(CCSP.model, newdata = test.CCSP, type = "response", se.fit = TRUE))
colnames(test.CCSP.predicted)[colnames(test.CCSP.predicted)=="fitted.fit"] <- "Predicted"
ROC.CCSP<-roc(test.CCSP.predicted$Detected, test.CCSP.predicted$Predicted) 
ROC.CCSP
ROCplot.CCSP<-plot(ROC.CCSP)

####BEKI
BEKI<-ARU[ARU$Sound=="BEKI", ]


indexes=sample(1:nrow(BEKI), size=0.3*nrow(BEKI))
test.BEKI=BEKI[indexes,]
train.BEKI=BEKI[-indexes,]
nrow(test.BEKI)    #255 obs
nrow(train.BEKI)    #595 obs


BEKI.model<-glm(Detected~x*Habitat + Type, data=train.BEKI, family=binomial)
top.model[["BEKI"]]=BEKI.model


test.BEKI.predicted<-cbind(test.BEKI, fitted = predict(BEKI.model, newdata = test.BEKI, type = "response", se.fit = TRUE))
colnames(test.BEKI.predicted)[colnames(test.BEKI.predicted)=="fitted.fit"] <- "Predicted"
ROC.BEKI<-roc(test.BEKI.predicted$Detected, test.BEKI.predicted$Predicted) 
ROC.BEKI
ROCplot.BEKI<-plot(ROC.BEKI)

####PISI
PISI<-ARU[ARU$Sound=="PISI", ]


indexes=sample(1:nrow(PISI), size=0.3*nrow(PISI))
test.PISI=PISI[indexes,]
train.PISI=PISI[-indexes,]
nrow(test.PISI)    #255 obs
nrow(train.PISI)    #595 obs


PISI.model<-glm(Detected~x*Habitat + Type, data=train.PISI, family=binomial)
top.model[["PISI"]]=PISI.model


test.PISI.predicted<-cbind(test.PISI, fitted = predict(PISI.model, newdata = test.PISI, type = "response", se.fit = TRUE))
colnames(test.PISI.predicted)[colnames(test.PISI.predicted)=="fitted.fit"] <- "Predicted"
ROC.PISI<-roc(test.PISI.predicted$Detected, test.PISI.predicted$Predicted) 
ROC.PISI
ROCplot.PISI<-plot(ROC.PISI)

####TEWA
TEWA<-ARU[ARU$Sound=="TEWA", ]


indexes=sample(1:nrow(TEWA), size=0.3*nrow(TEWA))
test.TEWA=TEWA[indexes,]
train.TEWA=TEWA[-indexes,]
nrow(test.TEWA)    #255 obs
nrow(train.TEWA)    #595 obs


TEWA.model<-glm(Detected~x*Habitat + Type, data=train.TEWA, family=binomial)
top.model[["TEWA"]]=TEWA.model

test.TEWA.predicted<-cbind(test.TEWA, fitted = predict(TEWA.model, newdata = test.TEWA, type = "response", se.fit = TRUE))
colnames(test.TEWA.predicted)[colnames(test.TEWA.predicted)=="fitted.fit"] <- "Predicted"
ROC.TEWA<-roc(test.TEWA.predicted$Detected, test.TEWA.predicted$Predicted) 
ROC.TEWA
ROCplot.TEWA<-plot(ROC.TEWA)

####BLWA
BLWA<-ARU[ARU$Sound=="BLWA", ]


indexes=sample(1:nrow(BLWA), size=0.3*nrow(BLWA))
test.BLWA=BLWA[indexes,]
train.BLWA=BLWA[-indexes,]
nrow(test.BLWA)    #255 obs
nrow(train.BLWA)    #595 obs


BLWA.model<-glm(Detected~x*Habitat + Type, data=train.BLWA, family=binomial)
top.model[["BLWA"]]=BLWA.model

test.BLWA.predicted<-cbind(test.BLWA, fitted = predict(BLWA.model, newdata = test.BLWA, type = "response", se.fit = TRUE))
colnames(test.BLWA.predicted)[colnames(test.BLWA.predicted)=="fitted.fit"] <- "Predicted"
ROC.BLWA<-roc(test.BLWA.predicted$Detected, test.BLWA.predicted$Predicted) 
ROC.BLWA
ROCplot.BLWA<-plot(ROC.BLWA)

####LISP
LISP<-ARU[ARU$Sound=="LISP", ]


indexes=sample(1:nrow(LISP), size=0.3*nrow(LISP))
test.LISP=LISP[indexes,]
train.LISP=LISP[-indexes,]
nrow(test.LISP)    #255 obs
nrow(train.LISP)    #595 obs


LISP.model<-glm(Detected~x*Habitat + Type, data=train.LISP, family=binomial)
top.model[["LISP"]]=LISP.model

test.LISP.predicted<-cbind(test.LISP, fitted = predict(LISP.model, newdata = test.LISP, type = "response", se.fit = TRUE))
colnames(test.LISP.predicted)[colnames(test.LISP.predicted)=="fitted.fit"] <- "Predicted"
ROC.LISP<-roc(test.LISP.predicted$Detected, test.LISP.predicted$Predicted) 
ROC.LISP
ROCplot.LISP<-plot(ROC.LISP)

####BAWW
BAWW<-ARU[ARU$Sound=="BAWW", ]


indexes=sample(1:nrow(BAWW), size=0.3*nrow(BAWW))
test.BAWW=BAWW[indexes,]
train.BAWW=BAWW[-indexes,]
nrow(test.BAWW)    #255 obs
nrow(train.BAWW)    #595 obs


BAWW.model<-glm(Detected~x*Habitat + Type, data=train.BAWW, family=binomial)
top.model[["BAWW"]]=BAWW.model

test.BAWW.predicted<-cbind(test.BAWW, fitted = predict(BAWW.model, newdata = test.BAWW, type = "response", se.fit = TRUE))
colnames(test.BAWW.predicted)[colnames(test.BAWW.predicted)=="fitted.fit"] <- "Predicted"
ROC.BAWW<-roc(test.BAWW.predicted$Detected, test.BAWW.predicted$Predicted) 
ROC.BAWW
ROCplot.BAWW<-plot(ROC.BAWW)

####WAVI
WAVI<-ARU[ARU$Sound=="WAVI", ]


indexes=sample(1:nrow(WAVI), size=0.3*nrow(WAVI))
test.WAVI=WAVI[indexes,]
train.WAVI=WAVI[-indexes,]
nrow(test.WAVI)    #255 obs
nrow(train.WAVI)    #595 obs


WAVI.model<-glm(Detected~x*Habitat + Type, data=train.WAVI, family=binomial)
top.model[["WAVI"]]=WAVI.model

test.WAVI.predicted<-cbind(test.WAVI, fitted = predict(WAVI.model, newdata = test.WAVI, type = "response", se.fit = TRUE))
colnames(test.WAVI.predicted)[colnames(test.WAVI.predicted)=="fitted.fit"] <- "Predicted"
ROC.WAVI<-roc(test.WAVI.predicted$Detected, test.WAVI.predicted$Predicted) 
ROC.WAVI
ROCplot.WAVI<-plot(ROC.WAVI)

####OVEN
OVEN<-ARU[ARU$Sound=="OVEN", ]


indexes=sample(1:nrow(OVEN), size=0.3*nrow(OVEN))
test.OVEN=OVEN[indexes,]
train.OVEN=OVEN[-indexes,]
nrow(test.OVEN)    #255 obs
nrow(train.OVEN)    #595 obs


OVEN.model<-glm(Detected~x*Habitat + Type, data=train.OVEN, family=binomial)
top.model[["OVEN"]]=OVEN.model

test.OVEN.predicted<-cbind(test.OVEN, fitted = predict(OVEN.model, newdata = test.OVEN, type = "response", se.fit = TRUE))
colnames(test.OVEN.predicted)[colnames(test.OVEN.predicted)=="fitted.fit"] <- "Predicted"
ROC.OVEN<-roc(test.OVEN.predicted$Detected, test.OVEN.predicted$Predicted) 
ROC.OVEN
ROCplot.OVEN<-plot(ROC.OVEN)

####DEJU
DEJU<-ARU[ARU$Sound=="DEJU", ]


indexes=sample(1:nrow(DEJU), size=0.3*nrow(DEJU))
test.DEJU=DEJU[indexes,]
train.DEJU=DEJU[-indexes,]
nrow(test.DEJU)    #255 obs
nrow(train.DEJU)    #595 obs


DEJU.model<-glm(Detected~x*Habitat + Type, data=train.DEJU, family=binomial)
top.model[["DEJU"]]=DEJU.model

test.DEJU.predicted<-cbind(test.DEJU, fitted = predict(DEJU.model, newdata = test.DEJU, type = "response", se.fit = TRUE))
colnames(test.DEJU.predicted)[colnames(test.DEJU.predicted)=="fitted.fit"] <- "Predicted"
ROC.DEJU<-roc(test.DEJU.predicted$Detected, test.DEJU.predicted$Predicted) 
ROC.DEJU
ROCplot.DEJU<-plot(ROC.DEJU)

####BHCO
BHCO<-ARU[ARU$Sound=="BHCO", ]


indexes=sample(1:nrow(BHCO), size=0.3*nrow(BHCO))
test.BHCO=BHCO[indexes,]
train.BHCO=BHCO[-indexes,]
nrow(test.BHCO)    #255 obs
nrow(train.BHCO)    #595 obs


BHCO.model<-glm(Detected~x*Habitat + Type, data=train.BHCO, family=binomial)
top.model[["BHCO"]]=BHCO.model

test.BHCO.predicted<-cbind(test.BHCO, fitted = predict(BHCO.model, newdata = test.BHCO, type = "response", se.fit = TRUE))
colnames(test.BHCO.predicted)[colnames(test.BHCO.predicted)=="fitted.fit"] <- "Predicted"
ROC.BHCO<-roc(test.BHCO.predicted$Detected, test.BHCO.predicted$Predicted) 
ROC.BHCO
ROCplot.BHCO<-plot(ROC.BHCO)

####T4000Hz
T4000Hz<-ARU[ARU$Sound=="4000Hz", ]


indexes=sample(1:nrow(T4000Hz), size=0.3*nrow(T4000Hz))
test.T4000Hz=T4000Hz[indexes,]
train.T4000Hz=T4000Hz[-indexes,]
nrow(test.T4000Hz)    #255 obs
nrow(train.T4000Hz)    #595 obs


T4000Hz.model<-glm(Detected~x*Habitat + Type, data=train.T4000Hz, family=binomial)
top.model[["4000Hz"]]=T4000Hz.model

test.T4000Hz.predicted<-cbind(test.T4000Hz, fitted = predict(T4000Hz.model, newdata = test.T4000Hz, type = "response", se.fit = TRUE))
colnames(test.T4000Hz.predicted)[colnames(test.T4000Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T4000Hz<-roc(test.T4000Hz.predicted$Detected, test.T4000Hz.predicted$Predicted) 
ROC.T4000Hz
ROCplot.T4000Hz<-plot(ROC.T4000Hz)

####T2828Hz
T2828Hz<-ARU[ARU$Sound=="2828Hz", ]


indexes=sample(1:nrow(T2828Hz), size=0.3*nrow(T2828Hz))
test.T2828Hz=T2828Hz[indexes,]
train.T2828Hz=T2828Hz[-indexes,]
nrow(test.T2828Hz)    #255 obs
nrow(train.T2828Hz)    #595 obs


T2828Hz.model<-glm(Detected~x*Habitat + Type, data=train.T2828Hz, family=binomial)
top.model[["2828Hz"]]=T2828Hz.model

test.T2828Hz.predicted<-cbind(test.T2828Hz, fitted = predict(T2828Hz.model, newdata = test.T2828Hz, type = "response", se.fit = TRUE))
colnames(test.T2828Hz.predicted)[colnames(test.T2828Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T2828Hz<-roc(test.T2828Hz.predicted$Detected, test.T2828Hz.predicted$Predicted) 
ROC.T2828Hz
ROCplot.T2828Hz<-plot(ROC.T2828Hz)

####T5656Hz
T5656Hz<-ARU[ARU$Sound=="5656Hz", ]


indexes=sample(1:nrow(T5656Hz), size=0.3*nrow(T5656Hz))
test.T5656Hz=T5656Hz[indexes,]
train.T5656Hz=T5656Hz[-indexes,]
nrow(test.T5656Hz)    #255 obs
nrow(train.T5656Hz)    #595 obs


T5656Hz.model<-glm(Detected~x*Habitat + Type, data=train.T5656Hz, family=binomial)
top.model[["5656Hz"]]=T5656Hz.model

test.T5656Hz.predicted<-cbind(test.T5656Hz, fitted = predict(T5656Hz.model, newdata = test.T5656Hz, type = "response", se.fit = TRUE))
colnames(test.T5656Hz.predicted)[colnames(test.T5656Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T5656Hz<-roc(test.T5656Hz.predicted$Detected, test.T5656Hz.predicted$Predicted) 
ROC.T5656Hz
ROCplot.T5656Hz<-plot(ROC.T5656Hz)

####T8000Hz
T8000Hz<-ARU[ARU$Sound=="8000Hz", ]


indexes=sample(1:nrow(T8000Hz), size=0.3*nrow(T8000Hz))
test.T8000Hz=T8000Hz[indexes,]
train.T8000Hz=T8000Hz[-indexes,]
nrow(test.T8000Hz)    #255 obs
nrow(train.T8000Hz)    #595 obs


T8000Hz.model<-glm(Detected~x*Habitat + Type, data=train.T8000Hz, family=binomial)
top.model[["8000Hz"]]=T8000Hz.model

test.T8000Hz.predicted<-cbind(test.T8000Hz, fitted = predict(T8000Hz.model, newdata = test.T8000Hz, type = "response", se.fit = TRUE))
colnames(test.T8000Hz.predicted)[colnames(test.T8000Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T8000Hz<-roc(test.T8000Hz.predicted$Detected, test.T8000Hz.predicted$Predicted) 
ROC.T8000Hz
ROCplot.T8000Hz<-plot(ROC.T8000Hz)

####x*Habitat + Type + Wind + Humidity - all weather variables####
####BOOW
BOOW<-ARU[ARU$Sound=="BOOW", ]


indexes=sample(1:nrow(BOOW), size=0.3*nrow(BOOW))
test.BOOW=BOOW[indexes,]
train.BOOW=BOOW[-indexes,]
nrow(test.BOOW)    #255 obs
nrow(train.BOOW)    #595 obs


BOOW.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.BOOW, family=binomial)
top.model[["BOOW"]]=BOOW.model

test.BOOW.predicted<-cbind(test.BOOW, fitted = predict(BOOW.model, newdata = test.BOOW, type = "response", se.fit = TRUE))
colnames(test.BOOW.predicted)[colnames(test.BOOW.predicted)=="fitted.fit"] <- "Predicted"
ROC.BOOW<-roc(test.BOOW.predicted$Detected, test.BOOW.predicted$Predicted) 
ROC.BOOW
ROCplot.BOOW<-plot(ROC.BOOW)

####NSWO
NSWO<-ARU[ARU$Sound=="NSWO", ]


indexes=sample(1:nrow(NSWO), size=0.3*nrow(NSWO))
test.NSWO=NSWO[indexes,]
train.NSWO=NSWO[-indexes,]
nrow(test.NSWO)    #255 obs
nrow(train.NSWO)    #595 obs


NSWO.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.NSWO, family=binomial)
top.model[["NSWO"]]=NSWO.model

test.NSWO.predicted<-cbind(test.NSWO, fitted = predict(NSWO.model, newdata = test.NSWO, type = "response", se.fit = TRUE))
colnames(test.NSWO.predicted)[colnames(test.NSWO.predicted)=="fitted.fit"] <- "Predicted"
ROC.NSWO<-roc(test.NSWO.predicted$Detected, test.NSWO.predicted$Predicted) 
ROC.NSWO
ROCplot.NSWO<-plot(ROC.NSWO)

####BBWA
BBWA<-ARU[ARU$Sound=="BBWA", ]


indexes=sample(1:nrow(BBWA), size=0.3*nrow(BBWA))
test.BBWA=BBWA[indexes,]
train.BBWA=BBWA[-indexes,]
nrow(test.BBWA)    #255 obs
nrow(train.BBWA)    #595 obs


BBWA.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.BBWA, family=binomial)
top.model[["BBWA"]]=BBWA.model

test.BBWA.predicted<-cbind(test.BBWA, fitted = predict(BBWA.model, newdata = test.BBWA, type = "response", se.fit = TRUE))
colnames(test.BBWA.predicted)[colnames(test.BBWA.predicted)=="fitted.fit"] <- "Predicted"
ROC.BBWA<-roc(test.BBWA.predicted$Detected, test.BBWA.predicted$Predicted) 
ROC.BBWA
ROCplot.BBWA<-plot(ROC.BBWA)

####LEOW
LEOW<-ARU[ARU$Sound=="LEOW", ]


indexes=sample(1:nrow(LEOW), size=0.3*nrow(LEOW))
test.LEOW=LEOW[indexes,]
train.LEOW=LEOW[-indexes,]
nrow(test.LEOW)    #255 obs
nrow(train.LEOW)    #595 obs


LEOW.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.LEOW, family=binomial)
top.model[["LEOW"]]=LEOW.model

test.LEOW.predicted<-cbind(test.LEOW, fitted = predict(LEOW.model, newdata = test.LEOW, type = "response", se.fit = TRUE))
colnames(test.LEOW.predicted)[colnames(test.LEOW.predicted)=="fitted.fit"] <- "Predicted"
ROC.LEOW<-roc(test.LEOW.predicted$Detected, test.LEOW.predicted$Predicted) 
ROC.LEOW
ROCplot.LEOW<-plot(ROC.LEOW)

####BADO
BADO<-ARU[ARU$Sound=="BADO", ]


indexes=sample(1:nrow(BADO), size=0.3*nrow(BADO))
test.BADO=BADO[indexes,]
train.BADO=BADO[-indexes,]
nrow(test.BADO)    #255 obs
nrow(train.BADO)    #595 obs


BADO.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.BADO, family=binomial)
top.model[["BADO"]]=BADO.model

test.BADO.predicted<-cbind(test.BADO, fitted = predict(BADO.model, newdata = test.BADO, type = "response", se.fit = TRUE))
colnames(test.BADO.predicted)[colnames(test.BADO.predicted)=="fitted.fit"] <- "Predicted"
ROC.BADO<-roc(test.BADO.predicted$Detected, test.BADO.predicted$Predicted) 
ROC.BADO
ROCplot.BADO<-plot(ROC.BADO)

####WETO
WETO<-ARU[ARU$Sound=="WETO", ]


indexes=sample(1:nrow(WETO), size=0.3*nrow(WETO))
test.WETO=WETO[indexes,]
train.WETO=WETO[-indexes,]
nrow(test.WETO)    #255 obs
nrow(train.WETO)    #595 obs


WETO.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.WETO, family=binomial)
top.model[["WETO"]]=WETO.model

test.WETO.predicted<-cbind(test.WETO, fitted = predict(WETO.model, newdata = test.WETO, type = "response", se.fit = TRUE))
colnames(test.WETO.predicted)[colnames(test.WETO.predicted)=="fitted.fit"] <- "Predicted"
ROC.WETO<-roc(test.WETO.predicted$Detected, test.WETO.predicted$Predicted) 
ROC.WETO
ROCplot.WETO<-plot(ROC.WETO)

####CORA
CORA<-ARU[ARU$Sound=="CORA", ]


indexes=sample(1:nrow(CORA), size=0.3*nrow(CORA))
test.CORA=CORA[indexes,]
train.CORA=CORA[-indexes,]
nrow(test.CORA)    #255 obs
nrow(train.CORA)    #595 obs


CORA.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.CORA, family=binomial)
top.model[["CORA"]]=CORA.model

test.CORA.predicted<-cbind(test.CORA, fitted = predict(CORA.model, newdata = test.CORA, type = "response", se.fit = TRUE))
colnames(test.CORA.predicted)[colnames(test.CORA.predicted)=="fitted.fit"] <- "Predicted"
ROC.CORA<-roc(test.CORA.predicted$Detected, test.CORA.predicted$Predicted) 
ROC.CORA
ROCplot.CORA<-plot(ROC.CORA)

####RBNU
RBNU<-ARU[ARU$Sound=="RBNU", ]


indexes=sample(1:nrow(RBNU), size=0.3*nrow(RBNU))
test.RBNU=RBNU[indexes,]
train.RBNU=RBNU[-indexes,]
nrow(test.RBNU)    #255 obs
nrow(train.RBNU)    #595 obs


RBNU.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.RBNU, family=binomial)
top.model[["RBNU"]]=RBNU.model

test.RBNU.predicted<-cbind(test.RBNU, fitted = predict(RBNU.model, newdata = test.RBNU, type = "response", se.fit = TRUE))
colnames(test.RBNU.predicted)[colnames(test.RBNU.predicted)=="fitted.fit"] <- "Predicted"
ROC.RBNU<-roc(test.RBNU.predicted$Detected, test.RBNU.predicted$Predicted) 
ROC.RBNU
ROCplot.RBNU<-plot(ROC.RBNU)

####GGOW
GGOW<-ARU[ARU$Sound=="GGOW", ]


indexes=sample(1:nrow(GGOW), size=0.3*nrow(GGOW))
test.GGOW=GGOW[indexes,]
train.GGOW=GGOW[-indexes,]
nrow(test.GGOW)    #255 obs
nrow(train.GGOW)    #595 obs


GGOW.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.GGOW, family=binomial)
top.model[["GGOW"]]=GGOW.model

test.GGOW.predicted<-cbind(test.GGOW, fitted = predict(GGOW.model, newdata = test.GGOW, type = "response", se.fit = TRUE))
colnames(test.GGOW.predicted)[colnames(test.GGOW.predicted)=="fitted.fit"] <- "Predicted"
ROC.GGOW<-roc(test.GGOW.predicted$Detected, test.GGOW.predicted$Predicted) 
ROC.GGOW
ROCplot.GGOW<-plot(ROC.GGOW)

####WTSP
WTSP<-ARU[ARU$Sound=="WTSP", ]


indexes=sample(1:nrow(WTSP), size=0.3*nrow(WTSP))
test.WTSP=WTSP[indexes,]
train.WTSP=WTSP[-indexes,]
nrow(test.WTSP)    #255 obs
nrow(train.WTSP)    #595 obs


WTSP.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.WTSP, family=binomial)
top.model[["WTSP"]]=WTSP.model

test.WTSP.predicted<-cbind(test.WTSP, fitted = predict(WTSP.model, newdata = test.WTSP, type = "response", se.fit = TRUE))
colnames(test.WTSP.predicted)[colnames(test.WTSP.predicted)=="fitted.fit"] <- "Predicted"
ROC.WTSP<-roc(test.WTSP.predicted$Detected, test.WTSP.predicted$Predicted) 
ROC.WTSP
ROCplot.WTSP<-plot(ROC.WTSP)

####OSFL
OSFL<-ARU[ARU$Sound=="OSFL", ]


indexes=sample(1:nrow(OSFL), size=0.3*nrow(OSFL))
test.OSFL=OSFL[indexes,]
train.OSFL=OSFL[-indexes,]
nrow(test.OSFL)    #255 obs
nrow(train.OSFL)    #595 obs


OSFL.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.OSFL, family=binomial)
top.model[["OSFL"]]=OSFL.model

test.OSFL.predicted<-cbind(test.OSFL, fitted = predict(OSFL.model, newdata = test.OSFL, type = "response", se.fit = TRUE))
colnames(test.OSFL.predicted)[colnames(test.OSFL.predicted)=="fitted.fit"] <- "Predicted"
ROC.OSFL<-roc(test.OSFL.predicted$Detected, test.OSFL.predicted$Predicted) 
ROC.OSFL
ROCplot.OSFL<-plot(ROC.OSFL)

####RBGR
RBGR<-ARU[ARU$Sound=="RBGR", ]


indexes=sample(1:nrow(RBGR), size=0.3*nrow(RBGR))
test.RBGR=RBGR[indexes,]
train.RBGR=RBGR[-indexes,]
nrow(test.RBGR)    #255 obs
nrow(train.RBGR)    #595 obs


RBGR.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.RBGR, family=binomial)
top.model[["RBGR"]]=RBGR.model

test.RBGR.predicted<-cbind(test.RBGR, fitted = predict(RBGR.model, newdata = test.RBGR, type = "response", se.fit = TRUE))
colnames(test.RBGR.predicted)[colnames(test.RBGR.predicted)=="fitted.fit"] <- "Predicted"
ROC.RBGR<-roc(test.RBGR.predicted$Detected, test.RBGR.predicted$Predicted) 
ROC.RBGR
ROCplot.RBGR<-plot(ROC.RBGR)

####CATO
CATO<-ARU[ARU$Sound=="CATO", ]


indexes=sample(1:nrow(CATO), size=0.3*nrow(CATO))
test.CATO=CATO[indexes,]
train.CATO=CATO[-indexes,]
nrow(test.CATO)    #255 obs
nrow(train.CATO)    #595 obs


CATO.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.CATO, family=binomial)
top.model[["CATO"]]=CATO.model

test.CATO.predicted<-cbind(test.CATO, fitted = predict(CATO.model, newdata = test.CATO, type = "response", se.fit = TRUE))
colnames(test.CATO.predicted)[colnames(test.CATO.predicted)=="fitted.fit"] <- "Predicted"
ROC.CATO<-roc(test.CATO.predicted$Detected, test.CATO.predicted$Predicted) 
ROC.CATO
ROCplot.CATO<-plot(ROC.CATO)

####T2000Hz
T2000Hz<-ARU[ARU$Sound=="2000Hz", ]


indexes=sample(1:nrow(T2000Hz), size=0.3*nrow(T2000Hz))
test.T2000Hz=T2000Hz[indexes,]
train.T2000Hz=T2000Hz[-indexes,]
nrow(test.T2000Hz)    #255 obs
nrow(train.T2000Hz)    #595 obs


T2000Hz.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.T2000Hz, family=binomial)
top.model[["2000Hz"]]=T2000Hz.model

test.T2000Hz.predicted<-cbind(test.T2000Hz, fitted = predict(T2000Hz.model, newdata = test.T2000Hz, type = "response", se.fit = TRUE))
colnames(test.T2000Hz.predicted)[colnames(test.T2000Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T2000Hz<-roc(test.T2000Hz.predicted$Detected, test.T2000Hz.predicted$Predicted) 
ROC.T2000Hz
ROCplot.T2000Hz<-plot(ROC.T2000Hz)

####T1000Hz
T1000Hz<-ARU[ARU$Sound=="1000Hz", ]


indexes=sample(1:nrow(T1000Hz), size=0.3*nrow(T1000Hz))
test.T1000Hz=T1000Hz[indexes,]
train.T1000Hz=T1000Hz[-indexes,]
nrow(test.T1000Hz)    #255 obs
nrow(train.T1000Hz)    #595 obs


T1000Hz.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.T1000Hz, family=binomial)
top.model[["1000Hz"]]=T1000Hz.model

test.T1000Hz.predicted<-cbind(test.T1000Hz, fitted = predict(T1000Hz.model, newdata = test.T1000Hz, type = "response", se.fit = TRUE))
colnames(test.T1000Hz.predicted)[colnames(test.T1000Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T1000Hz<-roc(test.T1000Hz.predicted$Detected, test.T1000Hz.predicted$Predicted) 
ROC.T1000Hz
ROCplot.T1000Hz<-plot(ROC.T1000Hz)

####T1414Hz
T1414Hz<-ARU[ARU$Sound=="1414Hz", ]


indexes=sample(1:nrow(T1414Hz), size=0.3*nrow(T1414Hz))
test.T1414Hz=T1414Hz[indexes,]
train.T1414Hz=T1414Hz[-indexes,]
nrow(test.T1414Hz)    #255 obs
nrow(train.T1414Hz)    #595 obs


T1414Hz.model<-glm(Detected~x*Habitat + Type + Wind + Humidity, data=train.T1414Hz, family=binomial)
top.model[["1414Hz"]]=T1414Hz.model

test.T1414Hz.predicted<-cbind(test.T1414Hz, fitted = predict(T1414Hz.model, newdata = test.T1414Hz, type = "response", se.fit = TRUE))
colnames(test.T1414Hz.predicted)[colnames(test.T1414Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T1414Hz<-roc(test.T1414Hz.predicted$Detected, test.T1414Hz.predicted$Predicted) 
ROC.T1414Hz
ROCplot.T1414Hz<-plot(ROC.T1414Hz)

####x*Habitat + Type + Wind - wind only####
####YERA
YERA<-ARU[ARU$Sound=="YERA", ]


indexes=sample(1:nrow(YERA), size=0.3*nrow(YERA))
test.YERA=YERA[indexes,]
train.YERA=YERA[-indexes,]
nrow(test.YERA)    #255 obs
nrow(train.YERA)    #595 obs


YERA.model<-glm(Detected~x*Habitat + Type + Wind, data=train.YERA, family=binomial)
top.model[["YERA"]]=YERA.model

test.YERA.predicted<-cbind(test.YERA, fitted = predict(YERA.model, newdata = test.YERA, type = "response", se.fit = TRUE))
colnames(test.YERA.predicted)[colnames(test.YERA.predicted)=="fitted.fit"] <- "Predicted"
ROC.YERA<-roc(test.YERA.predicted$Detected, test.YERA.predicted$Predicted) 
ROC.YERA
ROCplot.YERA<-plot(ROC.YERA)

####T11313Hz
T11313Hz<-ARU[ARU$Sound=="11313Hz", ]


indexes=sample(1:nrow(T11313Hz), size=0.3*nrow(T11313Hz))
test.T11313Hz=T11313Hz[indexes,]
train.T11313Hz=T11313Hz[-indexes,]
nrow(test.T11313Hz)    #255 obs
nrow(train.T11313Hz)    #595 obs


T11313Hz.model<-glm(Detected~x*Habitat + Type + Wind, data=train.T11313Hz, family=binomial)
top.model[["11313Hz"]]=T11313Hz.model

test.T11313Hz.predicted<-cbind(test.T11313Hz, fitted = predict(T11313Hz.model, newdata = test.T11313Hz, type = "response", se.fit = TRUE))
colnames(test.T11313Hz.predicted)[colnames(test.T11313Hz.predicted)=="fitted.fit"] <- "Predicted"
ROC.T11313Hz<-roc(test.T11313Hz.predicted$Detected, test.T11313Hz.predicted$Predicted) 
ROC.T11313Hz
ROCplot.T11313Hz<-plot(ROC.T11313Hz)

####Generate model coefficient table####

CFTABLE <- data.frame(Sound=character(),
                   Intercept=numeric(),
                   Distance=numeric(),
                   Conifer=numeric(),
                   Deciduous=numeric(),
                   InteractConifer=numeric(),
                   InteractDeciduous=numeric(),
                   RiverForks=numeric(),
                   SM2=numeric(),
                   SM3=numeric(),
                   Zoom=numeric(),
                   Wind=numeric(),
                   Humidity=numeric())

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BEKI",
                                     Intercept=coef(summary(top.model[["BEKI"]]))[1],
                                     Distance=coef(summary(top.model[["BEKI"]]))[2],
                                     Conifer=coef(summary(top.model[["BEKI"]]))[3],
                                     Deciduous=coef(summary(top.model[["BEKI"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BEKI"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["BEKI"]]))[10],
                                     RiverForks=coef(summary(top.model[["BEKI"]]))[5],
                                     SM2=coef(summary(top.model[["BEKI"]]))[6],
                                     SM3=coef(summary(top.model[["BEKI"]]))[7],
                                     Zoom=coef(summary(top.model[["BEKI"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="PISI",
                                     Intercept=coef(summary(top.model[["PISI"]]))[1],
                                     Distance=coef(summary(top.model[["PISI"]]))[2],
                                     Conifer=coef(summary(top.model[["PISI"]]))[3],
                                     Deciduous=coef(summary(top.model[["PISI"]]))[4],
                                     InteractConifer=coef(summary(top.model[["PISI"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["PISI"]]))[10],
                                     RiverForks=coef(summary(top.model[["PISI"]]))[5],
                                     SM2=coef(summary(top.model[["PISI"]]))[6],
                                     SM3=coef(summary(top.model[["PISI"]]))[7],
                                     Zoom=coef(summary(top.model[["PISI"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="TEWA",
                                     Intercept=coef(summary(top.model[["TEWA"]]))[1],
                                     Distance=coef(summary(top.model[["TEWA"]]))[2],
                                     Conifer=coef(summary(top.model[["TEWA"]]))[3],
                                     Deciduous=coef(summary(top.model[["TEWA"]]))[4],
                                     InteractConifer=coef(summary(top.model[["TEWA"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["TEWA"]]))[10],
                                     RiverForks=coef(summary(top.model[["TEWA"]]))[5],
                                     SM2=coef(summary(top.model[["TEWA"]]))[6],
                                     SM3=coef(summary(top.model[["TEWA"]]))[7],
                                     Zoom=coef(summary(top.model[["TEWA"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BLWA",
                                     Intercept=coef(summary(top.model[["BLWA"]]))[1],
                                     Distance=coef(summary(top.model[["BLWA"]]))[2],
                                     Conifer=coef(summary(top.model[["BLWA"]]))[3],
                                     Deciduous=coef(summary(top.model[["BLWA"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BLWA"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["BLWA"]]))[10],
                                     RiverForks=coef(summary(top.model[["BLWA"]]))[5],
                                     SM2=coef(summary(top.model[["BLWA"]]))[6],
                                     SM3=coef(summary(top.model[["BLWA"]]))[7],
                                     Zoom=coef(summary(top.model[["BLWA"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="LISP",
                                     Intercept=coef(summary(top.model[["LISP"]]))[1],
                                     Distance=coef(summary(top.model[["LISP"]]))[2],
                                     Conifer=coef(summary(top.model[["LISP"]]))[3],
                                     Deciduous=coef(summary(top.model[["LISP"]]))[4],
                                     InteractConifer=coef(summary(top.model[["LISP"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["LISP"]]))[10],
                                     RiverForks=coef(summary(top.model[["LISP"]]))[5],
                                     SM2=coef(summary(top.model[["LISP"]]))[6],
                                     SM3=coef(summary(top.model[["LISP"]]))[7],
                                     Zoom=coef(summary(top.model[["LISP"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="CCSP",
                                     Intercept=coef(summary(top.model[["CCSP"]]))[1],
                                     Distance=coef(summary(top.model[["CCSP"]]))[2],
                                     Conifer=coef(summary(top.model[["CCSP"]]))[3],
                                     Deciduous=coef(summary(top.model[["CCSP"]]))[4],
                                     InteractConifer=coef(summary(top.model[["CCSP"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["CCSP"]]))[10],
                                     RiverForks=coef(summary(top.model[["CCSP"]]))[5],
                                     SM2=coef(summary(top.model[["CCSP"]]))[6],
                                     SM3=coef(summary(top.model[["CCSP"]]))[7],
                                     Zoom=coef(summary(top.model[["CCSP"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BAWW",
                                     Intercept=coef(summary(top.model[["BAWW"]]))[1],
                                     Distance=coef(summary(top.model[["BAWW"]]))[2],
                                     Conifer=coef(summary(top.model[["BAWW"]]))[3],
                                     Deciduous=coef(summary(top.model[["BAWW"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BAWW"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["BAWW"]]))[10],
                                     RiverForks=coef(summary(top.model[["BAWW"]]))[5],
                                     SM2=coef(summary(top.model[["BAWW"]]))[6],
                                     SM3=coef(summary(top.model[["BAWW"]]))[7],
                                     Zoom=coef(summary(top.model[["BAWW"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="WAVI",
                                     Intercept=coef(summary(top.model[["WAVI"]]))[1],
                                     Distance=coef(summary(top.model[["WAVI"]]))[2],
                                     Conifer=coef(summary(top.model[["WAVI"]]))[3],
                                     Deciduous=coef(summary(top.model[["WAVI"]]))[4],
                                     InteractConifer=coef(summary(top.model[["WAVI"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["WAVI"]]))[10],
                                     RiverForks=coef(summary(top.model[["WAVI"]]))[5],
                                     SM2=coef(summary(top.model[["WAVI"]]))[6],
                                     SM3=coef(summary(top.model[["WAVI"]]))[7],
                                     Zoom=coef(summary(top.model[["WAVI"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="OVEN",
                                     Intercept=coef(summary(top.model[["OVEN"]]))[1],
                                     Distance=coef(summary(top.model[["OVEN"]]))[2],
                                     Conifer=coef(summary(top.model[["OVEN"]]))[3],
                                     Deciduous=coef(summary(top.model[["OVEN"]]))[4],
                                     InteractConifer=coef(summary(top.model[["OVEN"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["OVEN"]]))[10],
                                     RiverForks=coef(summary(top.model[["OVEN"]]))[5],
                                     SM2=coef(summary(top.model[["OVEN"]]))[6],
                                     SM3=coef(summary(top.model[["OVEN"]]))[7],
                                     Zoom=coef(summary(top.model[["OVEN"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="DEJU",
                                     Intercept=coef(summary(top.model[["DEJU"]]))[1],
                                     Distance=coef(summary(top.model[["DEJU"]]))[2],
                                     Conifer=coef(summary(top.model[["DEJU"]]))[3],
                                     Deciduous=coef(summary(top.model[["DEJU"]]))[4],
                                     InteractConifer=coef(summary(top.model[["DEJU"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["DEJU"]]))[10],
                                     RiverForks=coef(summary(top.model[["DEJU"]]))[5],
                                     SM2=coef(summary(top.model[["DEJU"]]))[6],
                                     SM3=coef(summary(top.model[["DEJU"]]))[7],
                                     Zoom=coef(summary(top.model[["DEJU"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BHCO",
                                     Intercept=coef(summary(top.model[["BHCO"]]))[1],
                                     Distance=coef(summary(top.model[["BHCO"]]))[2],
                                     Conifer=coef(summary(top.model[["BHCO"]]))[3],
                                     Deciduous=coef(summary(top.model[["BHCO"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BHCO"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["BHCO"]]))[10],
                                     RiverForks=coef(summary(top.model[["BHCO"]]))[5],
                                     SM2=coef(summary(top.model[["BHCO"]]))[6],
                                     SM3=coef(summary(top.model[["BHCO"]]))[7],
                                     Zoom=coef(summary(top.model[["BHCO"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="4000Hz",
                                     Intercept=coef(summary(top.model[["4000Hz"]]))[1],
                                     Distance=coef(summary(top.model[["4000Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["4000Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["4000Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["4000Hz"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["4000Hz"]]))[10],
                                     RiverForks=coef(summary(top.model[["4000Hz"]]))[5],
                                     SM2=coef(summary(top.model[["4000Hz"]]))[6],
                                     SM3=coef(summary(top.model[["4000Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["4000Hz"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="2828Hz",
                                     Intercept=coef(summary(top.model[["2828Hz"]]))[1],
                                     Distance=coef(summary(top.model[["2828Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["2828Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["2828Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["2828Hz"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["2828Hz"]]))[10],
                                     RiverForks=coef(summary(top.model[["2828Hz"]]))[5],
                                     SM2=coef(summary(top.model[["2828Hz"]]))[6],
                                     SM3=coef(summary(top.model[["2828Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["2828Hz"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="5656Hz",
                                     Intercept=coef(summary(top.model[["5656Hz"]]))[1],
                                     Distance=coef(summary(top.model[["5656Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["5656Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["5656Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["5656Hz"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["5656Hz"]]))[10],
                                     RiverForks=coef(summary(top.model[["5656Hz"]]))[5],
                                     SM2=coef(summary(top.model[["5656Hz"]]))[6],
                                     SM3=coef(summary(top.model[["5656Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["5656Hz"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="8000Hz",
                                     Intercept=coef(summary(top.model[["8000Hz"]]))[1],
                                     Distance=coef(summary(top.model[["8000Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["8000Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["8000Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["8000Hz"]]))[9],
                                     InteractDeciduous=coef(summary(top.model[["8000Hz"]]))[10],
                                     RiverForks=coef(summary(top.model[["8000Hz"]]))[5],
                                     SM2=coef(summary(top.model[["8000Hz"]]))[6],
                                     SM3=coef(summary(top.model[["8000Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["8000Hz"]]))[8],
                                     Wind=NA,
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BOOW",
                                     Intercept=coef(summary(top.model[["BOOW"]]))[1],
                                     Distance=coef(summary(top.model[["BOOW"]]))[2],
                                     Conifer=coef(summary(top.model[["BOOW"]]))[3],
                                     Deciduous=coef(summary(top.model[["BOOW"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BOOW"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["BOOW"]]))[12],
                                     RiverForks=coef(summary(top.model[["BOOW"]]))[5],
                                     SM2=coef(summary(top.model[["BOOW"]]))[6],
                                     SM3=coef(summary(top.model[["BOOW"]]))[7],
                                     Zoom=coef(summary(top.model[["BOOW"]]))[8],
                                     Wind=coef(summary(top.model[["BOOW"]]))[9],
                                     Humidity=coef(summary(top.model[["BOOW"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="NSWO",
                                     Intercept=coef(summary(top.model[["NSWO"]]))[1],
                                     Distance=coef(summary(top.model[["NSWO"]]))[2],
                                     Conifer=coef(summary(top.model[["NSWO"]]))[3],
                                     Deciduous=coef(summary(top.model[["NSWO"]]))[4],
                                     InteractConifer=coef(summary(top.model[["NSWO"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["NSWO"]]))[12],
                                     RiverForks=coef(summary(top.model[["NSWO"]]))[5],
                                     SM2=coef(summary(top.model[["NSWO"]]))[6],
                                     SM3=coef(summary(top.model[["NSWO"]]))[7],
                                     Zoom=coef(summary(top.model[["NSWO"]]))[8],
                                     Wind=coef(summary(top.model[["NSWO"]]))[9],
                                     Humidity=coef(summary(top.model[["NSWO"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BBWA",
                                     Intercept=coef(summary(top.model[["BBWA"]]))[1],
                                     Distance=coef(summary(top.model[["BBWA"]]))[2],
                                     Conifer=coef(summary(top.model[["BBWA"]]))[3],
                                     Deciduous=coef(summary(top.model[["BBWA"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BBWA"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["BBWA"]]))[12],
                                     RiverForks=coef(summary(top.model[["BBWA"]]))[5],
                                     SM2=coef(summary(top.model[["BBWA"]]))[6],
                                     SM3=coef(summary(top.model[["BBWA"]]))[7],
                                     Zoom=coef(summary(top.model[["BBWA"]]))[8],
                                     Wind=coef(summary(top.model[["BBWA"]]))[9],
                                     Humidity=coef(summary(top.model[["BBWA"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="LEOW",
                                     Intercept=coef(summary(top.model[["LEOW"]]))[1],
                                     Distance=coef(summary(top.model[["LEOW"]]))[2],
                                     Conifer=coef(summary(top.model[["LEOW"]]))[3],
                                     Deciduous=coef(summary(top.model[["LEOW"]]))[4],
                                     InteractConifer=coef(summary(top.model[["LEOW"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["LEOW"]]))[12],
                                     RiverForks=coef(summary(top.model[["LEOW"]]))[5],
                                     SM2=coef(summary(top.model[["LEOW"]]))[6],
                                     SM3=coef(summary(top.model[["LEOW"]]))[7],
                                     Zoom=coef(summary(top.model[["LEOW"]]))[8],
                                     Wind=coef(summary(top.model[["LEOW"]]))[9],
                                     Humidity=coef(summary(top.model[["LEOW"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="BADO",
                                     Intercept=coef(summary(top.model[["BADO"]]))[1],
                                     Distance=coef(summary(top.model[["BADO"]]))[2],
                                     Conifer=coef(summary(top.model[["BADO"]]))[3],
                                     Deciduous=coef(summary(top.model[["BADO"]]))[4],
                                     InteractConifer=coef(summary(top.model[["BADO"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["BADO"]]))[12],
                                     RiverForks=coef(summary(top.model[["BADO"]]))[5],
                                     SM2=coef(summary(top.model[["BADO"]]))[6],
                                     SM3=coef(summary(top.model[["BADO"]]))[7],
                                     Zoom=coef(summary(top.model[["BADO"]]))[8],
                                     Wind=coef(summary(top.model[["BADO"]]))[9],
                                     Humidity=coef(summary(top.model[["BADO"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="WETO",
                                     Intercept=coef(summary(top.model[["WETO"]]))[1],
                                     Distance=coef(summary(top.model[["WETO"]]))[2],
                                     Conifer=coef(summary(top.model[["WETO"]]))[3],
                                     Deciduous=coef(summary(top.model[["WETO"]]))[4],
                                     InteractConifer=coef(summary(top.model[["WETO"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["WETO"]]))[12],
                                     RiverForks=coef(summary(top.model[["WETO"]]))[5],
                                     SM2=coef(summary(top.model[["WETO"]]))[6],
                                     SM3=coef(summary(top.model[["WETO"]]))[7],
                                     Zoom=coef(summary(top.model[["WETO"]]))[8],
                                     Wind=coef(summary(top.model[["WETO"]]))[9],
                                     Humidity=coef(summary(top.model[["WETO"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="CORA",
                                     Intercept=coef(summary(top.model[["CORA"]]))[1],
                                     Distance=coef(summary(top.model[["CORA"]]))[2],
                                     Conifer=coef(summary(top.model[["CORA"]]))[3],
                                     Deciduous=coef(summary(top.model[["CORA"]]))[4],
                                     InteractConifer=coef(summary(top.model[["CORA"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["CORA"]]))[12],
                                     RiverForks=coef(summary(top.model[["CORA"]]))[5],
                                     SM2=coef(summary(top.model[["CORA"]]))[6],
                                     SM3=coef(summary(top.model[["CORA"]]))[7],
                                     Zoom=coef(summary(top.model[["CORA"]]))[8],
                                     Wind=coef(summary(top.model[["CORA"]]))[9],
                                     Humidity=coef(summary(top.model[["CORA"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="RBNU",
                                     Intercept=coef(summary(top.model[["RBNU"]]))[1],
                                     Distance=coef(summary(top.model[["RBNU"]]))[2],
                                     Conifer=coef(summary(top.model[["RBNU"]]))[3],
                                     Deciduous=coef(summary(top.model[["RBNU"]]))[4],
                                     InteractConifer=coef(summary(top.model[["RBNU"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["RBNU"]]))[12],
                                     RiverForks=coef(summary(top.model[["RBNU"]]))[5],
                                     SM2=coef(summary(top.model[["RBNU"]]))[6],
                                     SM3=coef(summary(top.model[["RBNU"]]))[7],
                                     Zoom=coef(summary(top.model[["RBNU"]]))[8],
                                     Wind=coef(summary(top.model[["RBNU"]]))[9],
                                     Humidity=coef(summary(top.model[["RBNU"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="GGOW",
                                     Intercept=coef(summary(top.model[["GGOW"]]))[1],
                                     Distance=coef(summary(top.model[["GGOW"]]))[2],
                                     Conifer=coef(summary(top.model[["GGOW"]]))[3],
                                     Deciduous=coef(summary(top.model[["GGOW"]]))[4],
                                     InteractConifer=coef(summary(top.model[["GGOW"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["GGOW"]]))[12],
                                     RiverForks=coef(summary(top.model[["GGOW"]]))[5],
                                     SM2=coef(summary(top.model[["GGOW"]]))[6],
                                     SM3=coef(summary(top.model[["GGOW"]]))[7],
                                     Zoom=coef(summary(top.model[["GGOW"]]))[8],
                                     Wind=coef(summary(top.model[["GGOW"]]))[9],
                                     Humidity=coef(summary(top.model[["GGOW"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="WTSP",
                                     Intercept=coef(summary(top.model[["WTSP"]]))[1],
                                     Distance=coef(summary(top.model[["WTSP"]]))[2],
                                     Conifer=coef(summary(top.model[["WTSP"]]))[3],
                                     Deciduous=coef(summary(top.model[["WTSP"]]))[4],
                                     InteractConifer=coef(summary(top.model[["WTSP"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["WTSP"]]))[12],
                                     RiverForks=coef(summary(top.model[["WTSP"]]))[5],
                                     SM2=coef(summary(top.model[["WTSP"]]))[6],
                                     SM3=coef(summary(top.model[["WTSP"]]))[7],
                                     Zoom=coef(summary(top.model[["WTSP"]]))[8],
                                     Wind=coef(summary(top.model[["WTSP"]]))[9],
                                     Humidity=coef(summary(top.model[["WTSP"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="OSFL",
                                     Intercept=coef(summary(top.model[["OSFL"]]))[1],
                                     Distance=coef(summary(top.model[["OSFL"]]))[2],
                                     Conifer=coef(summary(top.model[["OSFL"]]))[3],
                                     Deciduous=coef(summary(top.model[["OSFL"]]))[4],
                                     InteractConifer=coef(summary(top.model[["OSFL"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["OSFL"]]))[12],
                                     RiverForks=coef(summary(top.model[["OSFL"]]))[5],
                                     SM2=coef(summary(top.model[["OSFL"]]))[6],
                                     SM3=coef(summary(top.model[["OSFL"]]))[7],
                                     Zoom=coef(summary(top.model[["OSFL"]]))[8],
                                     Wind=coef(summary(top.model[["OSFL"]]))[9],
                                     Humidity=coef(summary(top.model[["OSFL"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="RBGR",
                                     Intercept=coef(summary(top.model[["RBGR"]]))[1],
                                     Distance=coef(summary(top.model[["RBGR"]]))[2],
                                     Conifer=coef(summary(top.model[["RBGR"]]))[3],
                                     Deciduous=coef(summary(top.model[["RBGR"]]))[4],
                                     InteractConifer=coef(summary(top.model[["RBGR"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["RBGR"]]))[12],
                                     RiverForks=coef(summary(top.model[["RBGR"]]))[5],
                                     SM2=coef(summary(top.model[["RBGR"]]))[6],
                                     SM3=coef(summary(top.model[["RBGR"]]))[7],
                                     Zoom=coef(summary(top.model[["RBGR"]]))[8],
                                     Wind=coef(summary(top.model[["RBGR"]]))[9],
                                     Humidity=coef(summary(top.model[["RBGR"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="CATO",
                                     Intercept=coef(summary(top.model[["CATO"]]))[1],
                                     Distance=coef(summary(top.model[["CATO"]]))[2],
                                     Conifer=coef(summary(top.model[["CATO"]]))[3],
                                     Deciduous=coef(summary(top.model[["CATO"]]))[4],
                                     InteractConifer=coef(summary(top.model[["CATO"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["CATO"]]))[12],
                                     RiverForks=coef(summary(top.model[["CATO"]]))[5],
                                     SM2=coef(summary(top.model[["CATO"]]))[6],
                                     SM3=coef(summary(top.model[["CATO"]]))[7],
                                     Zoom=coef(summary(top.model[["CATO"]]))[8],
                                     Wind=coef(summary(top.model[["CATO"]]))[9],
                                     Humidity=coef(summary(top.model[["CATO"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="2000Hz",
                                     Intercept=coef(summary(top.model[["2000Hz"]]))[1],
                                     Distance=coef(summary(top.model[["2000Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["2000Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["2000Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["2000Hz"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["2000Hz"]]))[12],
                                     RiverForks=coef(summary(top.model[["2000Hz"]]))[5],
                                     SM2=coef(summary(top.model[["2000Hz"]]))[6],
                                     SM3=coef(summary(top.model[["2000Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["2000Hz"]]))[8],
                                     Wind=coef(summary(top.model[["2000Hz"]]))[9],
                                     Humidity=coef(summary(top.model[["2000Hz"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="1000Hz",
                                     Intercept=coef(summary(top.model[["1000Hz"]]))[1],
                                     Distance=coef(summary(top.model[["1000Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["1000Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["1000Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["1000Hz"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["1000Hz"]]))[12],
                                     RiverForks=coef(summary(top.model[["1000Hz"]]))[5],
                                     SM2=coef(summary(top.model[["1000Hz"]]))[6],
                                     SM3=coef(summary(top.model[["1000Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["1000Hz"]]))[8],
                                     Wind=coef(summary(top.model[["1000Hz"]]))[9],
                                     Humidity=coef(summary(top.model[["1000Hz"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="1414Hz",
                                     Intercept=coef(summary(top.model[["1414Hz"]]))[1],
                                     Distance=coef(summary(top.model[["1414Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["1414Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["1414Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["1414Hz"]]))[11],
                                     InteractDeciduous=coef(summary(top.model[["1414Hz"]]))[12],
                                     RiverForks=coef(summary(top.model[["1414Hz"]]))[5],
                                     SM2=coef(summary(top.model[["1414Hz"]]))[6],
                                     SM3=coef(summary(top.model[["1414Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["1414Hz"]]))[8],
                                     Wind=coef(summary(top.model[["1414Hz"]]))[9],
                                     Humidity=coef(summary(top.model[["1414Hz"]]))[10]))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="YERA",
                                     Intercept=coef(summary(top.model[["YERA"]]))[1],
                                     Distance=coef(summary(top.model[["YERA"]]))[2],
                                     Conifer=coef(summary(top.model[["YERA"]]))[3],
                                     Deciduous=coef(summary(top.model[["YERA"]]))[4],
                                     InteractConifer=coef(summary(top.model[["YERA"]]))[10],
                                     InteractDeciduous=coef(summary(top.model[["YERA"]]))[11],
                                     RiverForks=coef(summary(top.model[["YERA"]]))[5],
                                     SM2=coef(summary(top.model[["YERA"]]))[6],
                                     SM3=coef(summary(top.model[["YERA"]]))[7],
                                     Zoom=coef(summary(top.model[["YERA"]]))[8],
                                     Wind=coef(summary(top.model[["YERA"]]))[9],
                                     Humidity=NA))

CFTABLE <- rbind(CFTABLE, data.frame(Sound="11313Hz",
                                     Intercept=coef(summary(top.model[["11313Hz"]]))[1],
                                     Distance=coef(summary(top.model[["11313Hz"]]))[2],
                                     Conifer=coef(summary(top.model[["11313Hz"]]))[3],
                                     Deciduous=coef(summary(top.model[["11313Hz"]]))[4],
                                     InteractConifer=coef(summary(top.model[["11313Hz"]]))[10],
                                     InteractDeciduous=coef(summary(top.model[["11313Hz"]]))[11],
                                     RiverForks=coef(summary(top.model[["11313Hz"]]))[5],
                                     SM2=coef(summary(top.model[["11313Hz"]]))[6],
                                     SM3=coef(summary(top.model[["11313Hz"]]))[7],
                                     Zoom=coef(summary(top.model[["11313Hz"]]))[8],
                                     Wind=coef(summary(top.model[["11313Hz"]]))[9],
                                     Humidity=NA))
               
write.csv(CFTABLE, file="SPECIES_model_coef.csv")

top.model1=list()
####
####EDR CALC####
edr=list()
confidence=list()
df1 <- data.frame(Sound=character(), Habitat=character(), Type=character(), EDR=numeric(), stringsAsFactors = FALSE)

####NO WEATHER EDR
#BEKI
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=BEKI, family=binomial("cloglog"))
summary(m)
top.model1[["BEKI"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BEKI.ROAD_HUM"]]<-edrROAD_HUM
edr[["BEKI.ROAD_RF"]]<-edrROAD_RF
edr[["BEKI.ROAD_SM2"]]<-edrROAD_SM2
edr[["BEKI.ROAD_SM3"]]<-edrROAD_SM3
edr[["BEKI.ROAD_ZM"]]<-edrROAD_ZM
edr[["BEKI.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BEKI.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BEKI.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BEKI.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BEKI.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BEKI.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BEKI.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BEKI.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BEKI.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BEKI.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BEKI$Habitat, x=BEKI$x, Type=BEKI$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BEKI.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BEKI.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#PISI
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=PISI, family=binomial("cloglog"))
summary(m)
top.model1[["PISI"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["PISI.ROAD_HUM"]]<-edrROAD_HUM
edr[["PISI.ROAD_RF"]]<-edrROAD_RF
edr[["PISI.ROAD_SM2"]]<-edrROAD_SM2
edr[["PISI.ROAD_SM3"]]<-edrROAD_SM3
edr[["PISI.ROAD_ZM"]]<-edrROAD_ZM
edr[["PISI.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["PISI.FORESTC_RF"]]<-edrFORESTC_RF
edr[["PISI.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["PISI.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["PISI.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["PISI.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["PISI.FORESTD_RF"]]<-edrFORESTD_RF
edr[["PISI.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["PISI.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["PISI.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=PISI$Habitat, x=PISI$x, Type=PISI$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["PISI.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["PISI.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#TEWA
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=TEWA, family=binomial("cloglog"))
summary(m)
top.model1[["TEWA"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["TEWA.ROAD_HUM"]]<-edrROAD_HUM
edr[["TEWA.ROAD_RF"]]<-edrROAD_RF
edr[["TEWA.ROAD_SM2"]]<-edrROAD_SM2
edr[["TEWA.ROAD_SM3"]]<-edrROAD_SM3
edr[["TEWA.ROAD_ZM"]]<-edrROAD_ZM
edr[["TEWA.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["TEWA.FORESTC_RF"]]<-edrFORESTC_RF
edr[["TEWA.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["TEWA.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["TEWA.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["TEWA.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["TEWA.FORESTD_RF"]]<-edrFORESTD_RF
edr[["TEWA.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["TEWA.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["TEWA.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=TEWA$Habitat, x=TEWA$x, Type=TEWA$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["TEWA.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["TEWA.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#BLWA
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=BLWA, family=binomial("cloglog"))
summary(m)
top.model1[["BLWA"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BLWA.ROAD_HUM"]]<-edrROAD_HUM
edr[["BLWA.ROAD_RF"]]<-edrROAD_RF
edr[["BLWA.ROAD_SM2"]]<-edrROAD_SM2
edr[["BLWA.ROAD_SM3"]]<-edrROAD_SM3
edr[["BLWA.ROAD_ZM"]]<-edrROAD_ZM
edr[["BLWA.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BLWA.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BLWA.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BLWA.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BLWA.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BLWA.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BLWA.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BLWA.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BLWA.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BLWA.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BLWA$Habitat, x=BLWA$x, Type=BLWA$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BLWA.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BLWA.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#LISP
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=LISP, family=binomial("cloglog"))
summary(m)
top.model1[["LISP"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["LISP.ROAD_HUM"]]<-edrROAD_HUM
edr[["LISP.ROAD_RF"]]<-edrROAD_RF
edr[["LISP.ROAD_SM2"]]<-edrROAD_SM2
edr[["LISP.ROAD_SM3"]]<-edrROAD_SM3
edr[["LISP.ROAD_ZM"]]<-edrROAD_ZM
edr[["LISP.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["LISP.FORESTC_RF"]]<-edrFORESTC_RF
edr[["LISP.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["LISP.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["LISP.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["LISP.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["LISP.FORESTD_RF"]]<-edrFORESTD_RF
edr[["LISP.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["LISP.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["LISP.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=LISP$Habitat, x=LISP$x, Type=LISP$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["LISP.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LISP.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#CCSP
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=CCSP, family=binomial("cloglog"))
summary(m)
top.model1[["CCSP"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["CCSP.ROAD_HUM"]]<-edrROAD_HUM
edr[["CCSP.ROAD_RF"]]<-edrROAD_RF
edr[["CCSP.ROAD_SM2"]]<-edrROAD_SM2
edr[["CCSP.ROAD_SM3"]]<-edrROAD_SM3
edr[["CCSP.ROAD_ZM"]]<-edrROAD_ZM
edr[["CCSP.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["CCSP.FORESTC_RF"]]<-edrFORESTC_RF
edr[["CCSP.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["CCSP.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["CCSP.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["CCSP.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["CCSP.FORESTD_RF"]]<-edrFORESTD_RF
edr[["CCSP.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["CCSP.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["CCSP.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CCSP$Habitat, x=CCSP$x, Type=CCSP$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["CCSP.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CCSP.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#BAWW
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=BAWW, family=binomial("cloglog"))
summary(m)
top.model1[["BAWW"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BAWW.ROAD_HUM"]]<-edrROAD_HUM
edr[["BAWW.ROAD_RF"]]<-edrROAD_RF
edr[["BAWW.ROAD_SM2"]]<-edrROAD_SM2
edr[["BAWW.ROAD_SM3"]]<-edrROAD_SM3
edr[["BAWW.ROAD_ZM"]]<-edrROAD_ZM
edr[["BAWW.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BAWW.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BAWW.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BAWW.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BAWW.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BAWW.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BAWW.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BAWW.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BAWW.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BAWW.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BAWW$Habitat, x=BAWW$x, Type=BAWW$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BAWW.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BAWW.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#WAVI
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=WAVI, family=binomial("cloglog"))
summary(m)
top.model1[["WAVI"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["WAVI.ROAD_HUM"]]<-edrROAD_HUM
edr[["WAVI.ROAD_RF"]]<-edrROAD_RF
edr[["WAVI.ROAD_SM2"]]<-edrROAD_SM2
edr[["WAVI.ROAD_SM3"]]<-edrROAD_SM3
edr[["WAVI.ROAD_ZM"]]<-edrROAD_ZM
edr[["WAVI.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["WAVI.FORESTC_RF"]]<-edrFORESTC_RF
edr[["WAVI.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["WAVI.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["WAVI.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["WAVI.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["WAVI.FORESTD_RF"]]<-edrFORESTD_RF
edr[["WAVI.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["WAVI.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["WAVI.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WAVI$Habitat, x=WAVI$x, Type=WAVI$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["WAVI.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WAVI.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#OVEN
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=OVEN, family=binomial("cloglog"))
summary(m)
top.model1[["OVEN"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["OVEN.ROAD_HUM"]]<-edrROAD_HUM
edr[["OVEN.ROAD_RF"]]<-edrROAD_RF
edr[["OVEN.ROAD_SM2"]]<-edrROAD_SM2
edr[["OVEN.ROAD_SM3"]]<-edrROAD_SM3
edr[["OVEN.ROAD_ZM"]]<-edrROAD_ZM
edr[["OVEN.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["OVEN.FORESTC_RF"]]<-edrFORESTC_RF
edr[["OVEN.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["OVEN.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["OVEN.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["OVEN.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["OVEN.FORESTD_RF"]]<-edrFORESTD_RF
edr[["OVEN.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["OVEN.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["OVEN.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=OVEN$Habitat, x=OVEN$x, Type=OVEN$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["OVEN.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OVEN.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#DEJU
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=DEJU, family=binomial("cloglog"))
summary(m)
top.model1[["DEJU"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["DEJU.ROAD_HUM"]]<-edrROAD_HUM
edr[["DEJU.ROAD_RF"]]<-edrROAD_RF
edr[["DEJU.ROAD_SM2"]]<-edrROAD_SM2
edr[["DEJU.ROAD_SM3"]]<-edrROAD_SM3
edr[["DEJU.ROAD_ZM"]]<-edrROAD_ZM
edr[["DEJU.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["DEJU.FORESTC_RF"]]<-edrFORESTC_RF
edr[["DEJU.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["DEJU.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["DEJU.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["DEJU.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["DEJU.FORESTD_RF"]]<-edrFORESTD_RF
edr[["DEJU.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["DEJU.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["DEJU.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=DEJU$Habitat, x=DEJU$x, Type=DEJU$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["DEJU.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["DEJU.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#BHCO
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=BHCO, family=binomial("cloglog"))
summary(m)
top.model1[["BHCO"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BHCO.ROAD_HUM"]]<-edrROAD_HUM
edr[["BHCO.ROAD_RF"]]<-edrROAD_RF
edr[["BHCO.ROAD_SM2"]]<-edrROAD_SM2
edr[["BHCO.ROAD_SM3"]]<-edrROAD_SM3
edr[["BHCO.ROAD_ZM"]]<-edrROAD_ZM
edr[["BHCO.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BHCO.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BHCO.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BHCO.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BHCO.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BHCO.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BHCO.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BHCO.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BHCO.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BHCO.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BHCO$Habitat, x=BHCO$x, Type=BHCO$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BHCO.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BHCO.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#4000Hz
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=T4000Hz, family=binomial("cloglog"))
summary(m)
top.model1[["4000Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["4000Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["4000Hz.ROAD_RF"]]<-edrROAD_RF
edr[["4000Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["4000Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["4000Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["4000Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["4000Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["4000Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["4000Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["4000Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["4000Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["4000Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["4000Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["4000Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["4000Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T4000Hz$Habitat, x=T4000Hz$x, Type=T4000Hz$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T4000Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T4000Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#2828Hz
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=T2828Hz, family=binomial("cloglog"))
summary(m)
top.model1[["2828Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["2828Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["2828Hz.ROAD_RF"]]<-edrROAD_RF
edr[["2828Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["2828Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["2828Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["2828Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["2828Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["2828Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["2828Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["2828Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["2828Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["2828Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["2828Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["2828Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["2828Hz.FORESTD_ZM"]]<-edrFORESTD_ZM
#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T2828Hz$Habitat, x=T2828Hz$x, Type=T2828Hz$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T2828Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2828Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#5656Hz
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=T5656Hz, family=binomial("cloglog"))
summary(m)
top.model1[["5656Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["5656Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["5656Hz.ROAD_RF"]]<-edrROAD_RF
edr[["5656Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["5656Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["5656Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["5656Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["5656Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["5656Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["5656Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["5656Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["5656Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["5656Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["5656Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["5656Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["5656Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T5656Hz$Habitat, x=T5656Hz$x, Type=T5656Hz$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T5656Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T5656Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#8000Hz
m <- glm(Detected ~ x + x:Habitat + x:Type -1, data=T8000Hz, family=binomial("cloglog"))
summary(m)
top.model1[["8000Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]
cROAD_RF <- cf[1]+cf[2]+cf[5]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]
cROAD_ZM <- cf[1]+cf[2]+cf[8]
cFORESTC_HUM <- cf[1]+cf[3]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]
cFORESTD_HUM <- cf[1]
cFORESTD_RF <- cf[1]+cf[5]
cFORESTD_SM2 <- cf[1]+cf[6]
cFORESTD_SM3 <- cf[1]+cf[7]
cFORESTD_ZM <- cf[1]+cf[8]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["8000Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["8000Hz.ROAD_RF"]]<-edrROAD_RF
edr[["8000Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["8000Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["8000Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["8000Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["8000Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["8000Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["8000Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["8000Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["8000Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["8000Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["8000Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["8000Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["8000Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T8000Hz$Habitat, x=T8000Hz$x, Type=T8000Hz$Type)
  m_star <- glm(y ~ x + x:Habitat + x:Type -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]
  cFORESTD_HUM_star <- cf_star[1]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T8000Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T8000Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

####Weather EDR####

#BOOW
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=BOOW, family=binomial("cloglog"))
summary(m)
top.model1[["BOOW"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BOOW.ROAD_HUM"]]<-edrROAD_HUM
edr[["BOOW.ROAD_RF"]]<-edrROAD_RF
edr[["BOOW.ROAD_SM2"]]<-edrROAD_SM2
edr[["BOOW.ROAD_SM3"]]<-edrROAD_SM3
edr[["BOOW.ROAD_ZM"]]<-edrROAD_ZM
edr[["BOOW.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BOOW.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BOOW.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BOOW.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BOOW.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BOOW.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BOOW.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BOOW.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BOOW.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BOOW.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BOOW$Habitat, x=BOOW$x, Type=BOOW$Type, Wind=BOOW$Wind, Humidity=BOOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BOOW.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BOOW.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#NSWO
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=NSWO, family=binomial("cloglog"))
summary(m)
top.model1[["NSWO"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["NSWO.ROAD_HUM"]]<-edrROAD_HUM
edr[["NSWO.ROAD_RF"]]<-edrROAD_RF
edr[["NSWO.ROAD_SM2"]]<-edrROAD_SM2
edr[["NSWO.ROAD_SM3"]]<-edrROAD_SM3
edr[["NSWO.ROAD_ZM"]]<-edrROAD_ZM
edr[["NSWO.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["NSWO.FORESTC_RF"]]<-edrFORESTC_RF
edr[["NSWO.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["NSWO.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["NSWO.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["NSWO.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["NSWO.FORESTD_RF"]]<-edrFORESTD_RF
edr[["NSWO.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["NSWO.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["NSWO.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=NSWO$Habitat, x=NSWO$x, Type=NSWO$Type, Wind=NSWO$Wind, Humidity=NSWO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["NSWO.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["NSWO.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#BBWA
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=BBWA, family=binomial("cloglog"))
summary(m)
top.model1[["BBWA"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BBWA.ROAD_HUM"]]<-edrROAD_HUM
edr[["BBWA.ROAD_RF"]]<-edrROAD_RF
edr[["BBWA.ROAD_SM2"]]<-edrROAD_SM2
edr[["BBWA.ROAD_SM3"]]<-edrROAD_SM3
edr[["BBWA.ROAD_ZM"]]<-edrROAD_ZM
edr[["BBWA.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BBWA.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BBWA.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BBWA.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BBWA.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BBWA.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BBWA.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BBWA.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BBWA.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BBWA.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BBWA$Habitat, x=BBWA$x, Type=BBWA$Type, Wind=BBWA$Wind, Humidity=BBWA$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BBWA.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BBWA.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#LEOW
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=LEOW, family=binomial("cloglog"))
summary(m)
top.model1[["LEOW"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["LEOW.ROAD_HUM"]]<-edrROAD_HUM
edr[["LEOW.ROAD_RF"]]<-edrROAD_RF
edr[["LEOW.ROAD_SM2"]]<-edrROAD_SM2
edr[["LEOW.ROAD_SM3"]]<-edrROAD_SM3
edr[["LEOW.ROAD_ZM"]]<-edrROAD_ZM
edr[["LEOW.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["LEOW.FORESTC_RF"]]<-edrFORESTC_RF
edr[["LEOW.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["LEOW.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["LEOW.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["LEOW.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["LEOW.FORESTD_RF"]]<-edrFORESTD_RF
edr[["LEOW.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["LEOW.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["LEOW.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=LEOW$Habitat, x=LEOW$x, Type=LEOW$Type, Wind=LEOW$Wind, Humidity=LEOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["LEOW.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["LEOW.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#BADO
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=BADO, family=binomial("cloglog"))
summary(m)
top.model1[["BADO"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["BADO.ROAD_HUM"]]<-edrROAD_HUM
edr[["BADO.ROAD_RF"]]<-edrROAD_RF
edr[["BADO.ROAD_SM2"]]<-edrROAD_SM2
edr[["BADO.ROAD_SM3"]]<-edrROAD_SM3
edr[["BADO.ROAD_ZM"]]<-edrROAD_ZM
edr[["BADO.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["BADO.FORESTC_RF"]]<-edrFORESTC_RF
edr[["BADO.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["BADO.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["BADO.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["BADO.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["BADO.FORESTD_RF"]]<-edrFORESTD_RF
edr[["BADO.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["BADO.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["BADO.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=BADO$Habitat, x=BADO$x, Type=BADO$Type, Wind=BADO$Wind, Humidity=BADO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["BADO.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["BADO.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#WETO
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=WETO, family=binomial("cloglog"))
summary(m)
top.model1[["WETO"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["WETO.ROAD_HUM"]]<-edrROAD_HUM
edr[["WETO.ROAD_RF"]]<-edrROAD_RF
edr[["WETO.ROAD_SM2"]]<-edrROAD_SM2
edr[["WETO.ROAD_SM3"]]<-edrROAD_SM3
edr[["WETO.ROAD_ZM"]]<-edrROAD_ZM
edr[["WETO.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["WETO.FORESTC_RF"]]<-edrFORESTC_RF
edr[["WETO.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["WETO.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["WETO.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["WETO.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["WETO.FORESTD_RF"]]<-edrFORESTD_RF
edr[["WETO.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["WETO.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["WETO.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WETO$Habitat, x=WETO$x, Type=WETO$Type, Wind=WETO$Wind, Humidity=WETO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["WETO.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WETO.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#CORA
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=CORA, family=binomial("cloglog"))
summary(m)
top.model1[["CORA"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["CORA.ROAD_HUM"]]<-edrROAD_HUM
edr[["CORA.ROAD_RF"]]<-edrROAD_RF
edr[["CORA.ROAD_SM2"]]<-edrROAD_SM2
edr[["CORA.ROAD_SM3"]]<-edrROAD_SM3
edr[["CORA.ROAD_ZM"]]<-edrROAD_ZM
edr[["CORA.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["CORA.FORESTC_RF"]]<-edrFORESTC_RF
edr[["CORA.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["CORA.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["CORA.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["CORA.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["CORA.FORESTD_RF"]]<-edrFORESTD_RF
edr[["CORA.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["CORA.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["CORA.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CORA$Habitat, x=CORA$x, Type=CORA$Type, Wind=CORA$Wind, Humidity=CORA$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["CORA.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CORA.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#RBNU
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=RBNU, family=binomial("cloglog"))
summary(m)
top.model1[["RBNU"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["RBNU.ROAD_HUM"]]<-edrROAD_HUM
edr[["RBNU.ROAD_RF"]]<-edrROAD_RF
edr[["RBNU.ROAD_SM2"]]<-edrROAD_SM2
edr[["RBNU.ROAD_SM3"]]<-edrROAD_SM3
edr[["RBNU.ROAD_ZM"]]<-edrROAD_ZM
edr[["RBNU.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["RBNU.FORESTC_RF"]]<-edrFORESTC_RF
edr[["RBNU.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["RBNU.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["RBNU.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["RBNU.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["RBNU.FORESTD_RF"]]<-edrFORESTD_RF
edr[["RBNU.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["RBNU.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["RBNU.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=RBNU$Habitat, x=RBNU$x, Type=RBNU$Type, Wind=RBNU$Wind, Humidity=RBNU$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["RBNU.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBNU.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#GGOW
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=GGOW, family=binomial("cloglog"))
summary(m)
top.model1[["GGOW"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["GGOW.ROAD_HUM"]]<-edrROAD_HUM
edr[["GGOW.ROAD_RF"]]<-edrROAD_RF
edr[["GGOW.ROAD_SM2"]]<-edrROAD_SM2
edr[["GGOW.ROAD_SM3"]]<-edrROAD_SM3
edr[["GGOW.ROAD_ZM"]]<-edrROAD_ZM
edr[["GGOW.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["GGOW.FORESTC_RF"]]<-edrFORESTC_RF
edr[["GGOW.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["GGOW.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["GGOW.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["GGOW.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["GGOW.FORESTD_RF"]]<-edrFORESTD_RF
edr[["GGOW.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["GGOW.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["GGOW.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=GGOW$Habitat, x=GGOW$x, Type=GGOW$Type, Wind=GGOW$Wind, Humidity=GGOW$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["GGOW.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["GGOW.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#WTSP
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=WTSP, family=binomial("cloglog"))
summary(m)
top.model1[["WTSP"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["WTSP.ROAD_HUM"]]<-edrROAD_HUM
edr[["WTSP.ROAD_RF"]]<-edrROAD_RF
edr[["WTSP.ROAD_SM2"]]<-edrROAD_SM2
edr[["WTSP.ROAD_SM3"]]<-edrROAD_SM3
edr[["WTSP.ROAD_ZM"]]<-edrROAD_ZM
edr[["WTSP.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["WTSP.FORESTC_RF"]]<-edrFORESTC_RF
edr[["WTSP.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["WTSP.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["WTSP.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["WTSP.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["WTSP.FORESTD_RF"]]<-edrFORESTD_RF
edr[["WTSP.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["WTSP.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["WTSP.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=WTSP$Habitat, x=WTSP$x, Type=WTSP$Type, Wind=WTSP$Wind, Humidity=WTSP$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["WTSP.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["WTSP.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#OSFL
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=OSFL, family=binomial("cloglog"))
summary(m)
top.model1[["OSFL"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["OSFL.ROAD_HUM"]]<-edrROAD_HUM
edr[["OSFL.ROAD_RF"]]<-edrROAD_RF
edr[["OSFL.ROAD_SM2"]]<-edrROAD_SM2
edr[["OSFL.ROAD_SM3"]]<-edrROAD_SM3
edr[["OSFL.ROAD_ZM"]]<-edrROAD_ZM
edr[["OSFL.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["OSFL.FORESTC_RF"]]<-edrFORESTC_RF
edr[["OSFL.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["OSFL.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["OSFL.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["OSFL.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["OSFL.FORESTD_RF"]]<-edrFORESTD_RF
edr[["OSFL.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["OSFL.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["OSFL.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=OSFL$Habitat, x=OSFL$x, Type=OSFL$Type, Wind=OSFL$Wind, Humidity=OSFL$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["OSFL.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["OSFL.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#RBGR
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=RBGR, family=binomial("cloglog"))
summary(m)
top.model1[["RBGR"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["RBGR.ROAD_HUM"]]<-edrROAD_HUM
edr[["RBGR.ROAD_RF"]]<-edrROAD_RF
edr[["RBGR.ROAD_SM2"]]<-edrROAD_SM2
edr[["RBGR.ROAD_SM3"]]<-edrROAD_SM3
edr[["RBGR.ROAD_ZM"]]<-edrROAD_ZM
edr[["RBGR.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["RBGR.FORESTC_RF"]]<-edrFORESTC_RF
edr[["RBGR.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["RBGR.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["RBGR.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["RBGR.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["RBGR.FORESTD_RF"]]<-edrFORESTD_RF
edr[["RBGR.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["RBGR.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["RBGR.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=RBGR$Habitat, x=RBGR$x, Type=RBGR$Type, Wind=RBGR$Wind, Humidity=RBGR$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["RBGR.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["RBGR.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#CATO
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=CATO, family=binomial("cloglog"))
summary(m)
top.model1[["CATO"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["CATO.ROAD_HUM"]]<-edrROAD_HUM
edr[["CATO.ROAD_RF"]]<-edrROAD_RF
edr[["CATO.ROAD_SM2"]]<-edrROAD_SM2
edr[["CATO.ROAD_SM3"]]<-edrROAD_SM3
edr[["CATO.ROAD_ZM"]]<-edrROAD_ZM
edr[["CATO.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["CATO.FORESTC_RF"]]<-edrFORESTC_RF
edr[["CATO.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["CATO.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["CATO.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["CATO.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["CATO.FORESTD_RF"]]<-edrFORESTD_RF
edr[["CATO.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["CATO.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["CATO.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=CATO$Habitat, x=CATO$x, Type=CATO$Type, Wind=CATO$Wind, Humidity=CATO$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["CATO.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["CATO.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#2000Hz
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=T2000Hz, family=binomial("cloglog"))
summary(m)
top.model1[["2000Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["2000Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["2000Hz.ROAD_RF"]]<-edrROAD_RF
edr[["2000Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["2000Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["2000Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["2000Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["2000Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["2000Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["2000Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["2000Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["2000Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["2000Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["2000Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["2000Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["2000Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T2000Hz$Habitat, x=T2000Hz$x, Type=T2000Hz$Type, Wind=T2000Hz$Wind, Humidity=T2000Hz$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T2000Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T2000Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#1000Hz
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=T1000Hz, family=binomial("cloglog"))
summary(m)
top.model1[["1000Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["1000Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["1000Hz.ROAD_RF"]]<-edrROAD_RF
edr[["1000Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["1000Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["1000Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["1000Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["1000Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["1000Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["1000Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["1000Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["1000Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["1000Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["1000Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["1000Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["1000Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T1000Hz$Habitat, x=T1000Hz$x, Type=T1000Hz$Type, Wind=T1000Hz$Wind, Humidity=T1000Hz$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T1000Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1000Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#1414Hz
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=T1414Hz, family=binomial("cloglog"))
summary(m)
top.model1[["1414Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]+cf[10]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]+cf[10]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]+cf[10]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]+cf[10]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]+cf[10]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]+cf[10]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]+cf[10]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]+cf[10]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]+cf[10]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]+cf[10]
cFORESTD_HUM <- cf[1]+cf[9]+cf[10]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]+cf[10]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]+cf[10]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]+cf[10]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]+cf[10]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["1414Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["1414Hz.ROAD_RF"]]<-edrROAD_RF
edr[["1414Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["1414Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["1414Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["1414Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["1414Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["1414Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["1414Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["1414Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["1414Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["1414Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["1414Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["1414Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["1414Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T1414Hz$Habitat, x=T1414Hz$x, Type=T1414Hz$Type, Wind=T1414Hz$Wind, Humidity=T1414Hz$Humidity)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind + x:Humidity -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]+cf_star[10]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]+cf_star[10]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]+cf_star[10]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]+cf_star[10]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]+cf_star[10]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]+cf_star[10]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]+cf_star[10]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]+cf_star[10]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]+cf_star[10]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]+cf_star[10]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]+cf_star[10]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T1414Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T1414Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

####Wind only EDR####

#YERA
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind -1, data=YERA, family=binomial("cloglog"))
summary(m)
top.model1[["YERA"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]
cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]
cFORESTD_HUM <- cf[1]+cf[9]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["YERA.ROAD_HUM"]]<-edrROAD_HUM
edr[["YERA.ROAD_RF"]]<-edrROAD_RF
edr[["YERA.ROAD_SM2"]]<-edrROAD_SM2
edr[["YERA.ROAD_SM3"]]<-edrROAD_SM3
edr[["YERA.ROAD_ZM"]]<-edrROAD_ZM
edr[["YERA.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["YERA.FORESTC_RF"]]<-edrFORESTC_RF
edr[["YERA.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["YERA.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["YERA.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["YERA.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["YERA.FORESTD_RF"]]<-edrFORESTD_RF
edr[["YERA.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["YERA.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["YERA.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=YERA$Habitat, x=YERA$x, Type=YERA$Type, Wind=YERA$Wind)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["YERA.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["YERA.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

#11313Hz
m <- glm(Detected ~ x + x:Habitat + x:Type + x:Wind -1, data=T11313Hz, family=binomial("cloglog"))
summary(m)
top.model1[["11313Hz"]]=m

#calculating coefficients for each habitat
cf <- coef(m)
cROAD_HUM <- cf[1]+cf[2]+cf[9]
cROAD_RF <- cf[1]+cf[2]+cf[5]+cf[9]
cROAD_SM2 <- cf[1]+cf[2]+cf[6]+cf[9]
cROAD_SM3 <- cf[1]+cf[2]+cf[7]+cf[9]
cROAD_ZM <- cf[1]+cf[2]+cf[8]+cf[9]
cFORESTC_HUM <- cf[1]+cf[3]+cf[9]
cFORESTC_RF <- cf[1]+cf[3]+cf[5]+cf[9]
  cFORESTC_SM2 <- cf[1]+cf[3]+cf[6]+cf[9]
cFORESTC_SM3 <- cf[1]+cf[3]+cf[7]+cf[9]
cFORESTC_ZM <- cf[1]+cf[3]+cf[8]+cf[9]
cFORESTD_HUM <- cf[1]+cf[9]
cFORESTD_RF <- cf[1]+cf[5]+cf[9]
cFORESTD_SM2 <- cf[1]+cf[6]+cf[9]
cFORESTD_SM3 <- cf[1]+cf[7]+cf[9]
cFORESTD_ZM <- cf[1]+cf[8]+cf[9]

#calculating EDR values
(edrROAD_HUM <- sqrt(1/cROAD_HUM))
(edrROAD_RF <- sqrt(1/cROAD_RF))
(edrROAD_SM2 <- sqrt(1/cROAD_SM2))
(edrROAD_SM3 <- sqrt(1/cROAD_SM3))
(edrROAD_ZM <- sqrt(1/cROAD_ZM))
(edrFORESTC_HUM <- sqrt(1/cFORESTC_HUM))
(edrFORESTC_RF <- sqrt(1/cFORESTC_RF))
(edrFORESTC_SM2 <- sqrt(1/cFORESTC_SM2))
(edrFORESTC_SM3 <- sqrt(1/cFORESTC_SM3))
(edrFORESTC_ZM <- sqrt(1/cFORESTC_ZM))
(edrFORESTD_HUM <- sqrt(1/cFORESTD_HUM))
(edrFORESTD_RF <- sqrt(1/cFORESTD_RF))
(edrFORESTD_SM2 <- sqrt(1/cFORESTD_SM2))
(edrFORESTD_SM3 <- sqrt(1/cFORESTD_SM3))
(edrFORESTD_ZM <- sqrt(1/cFORESTD_ZM))
edr[["11313Hz.ROAD_HUM"]]<-edrROAD_HUM
edr[["11313Hz.ROAD_RF"]]<-edrROAD_RF
edr[["11313Hz.ROAD_SM2"]]<-edrROAD_SM2
edr[["11313Hz.ROAD_SM3"]]<-edrROAD_SM3
edr[["11313Hz.ROAD_ZM"]]<-edrROAD_ZM
edr[["11313Hz.FORESTC_HUM"]]<-edrFORESTC_HUM
edr[["11313Hz.FORESTC_RF"]]<-edrFORESTC_RF
edr[["11313Hz.FORESTC_SM2"]]<-edrFORESTC_SM2
edr[["11313Hz.FORESTC_SM3"]]<-edrFORESTC_SM3
edr[["11313Hz.FORESTC_ZM"]]<-edrFORESTC_ZM
edr[["11313Hz.FORESTD_HUM"]]<-edrFORESTD_HUM
edr[["11313Hz.FORESTD_RF"]]<-edrFORESTD_RF
edr[["11313Hz.FORESTD_SM2"]]<-edrFORESTD_SM2
edr[["11313Hz.FORESTD_SM3"]]<-edrFORESTD_SM3
edr[["11313Hz.FORESTD_ZM"]]<-edrFORESTD_ZM

#Calculating 90% CIs
f <- fitted(m)
ROAD_HUM <- list()
ROAD_RF <- list()
ROAD_SM2 <- list()
ROAD_SM3 <- list()
ROAD_ZM <- list()
FORESTC_HUM <- list()
FORESTC_RF <- list()
FORESTC_SM2 <- list()
FORESTC_SM3 <- list()
FORESTC_ZM <- list()
FORESTD_HUM <- list()
FORESTD_RF <- list()
FORESTD_SM2 <- list()
FORESTD_SM3 <- list()
FORESTD_ZM <- list()

for (i in 1:1000) {
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y=y_star, Habitat=T11313Hz$Habitat, x=T11313Hz$x, Type=T11313Hz$Type, Wind=T11313Hz$Wind)
  m_star <- glm(y ~ x + x:Habitat + x:Type + x:Wind -1, data=df_star, family=binomial("cloglog"))
  
  cf_star <- coef(m_star)
  cROAD_HUM_star <- cf_star[1]+cf_star[2]+cf_star[9]
  cROAD_RF_star <- cf_star[1]+cf_star[2]+cf_star[5]+cf_star[9]
  cROAD_SM2_star <- cf_star[1]+cf_star[2]+cf_star[6]+cf_star[9]
  cROAD_SM3_star <- cf_star[1]+cf_star[2]+cf_star[7]+cf_star[9]
  cROAD_ZM_star <- cf_star[1]+cf_star[2]+cf_star[8]+cf_star[9]
  cFORESTC_HUM_star <- cf_star[1]+cf_star[3]+cf_star[9]
  cFORESTC_RF_star <- cf_star[1]+cf_star[3]+cf_star[5]+cf_star[9]
  cFORESTC_SM2_star <- cf_star[1]+cf_star[3]+cf_star[6]+cf_star[9]
  cFORESTC_SM3_star <- cf_star[1]+cf_star[3]+cf_star[7]+cf_star[9]
  cFORESTC_ZM_star <- cf_star[1]+cf_star[3]+cf_star[8]+cf_star[9]
  cFORESTD_HUM_star <- cf_star[1]+cf_star[9]
  cFORESTD_RF_star <- cf_star[1]+cf_star[5]+cf_star[9]
  cFORESTD_SM2_star <- cf_star[1]+cf_star[6]+cf_star[9]
  cFORESTD_SM3_star <- cf_star[1]+cf_star[7]+cf_star[9]
  cFORESTD_ZM_star <- cf_star[1]+cf_star[8]+cf_star[9]
  
  (edrROAD_HUM_star <- sqrt(1/cROAD_HUM_star))
  (edrROAD_RF_star <- sqrt(1/cROAD_RF_star))
  (edrROAD_SM2_star <- sqrt(1/cROAD_SM2_star))
  (edrROAD_SM3_star <- sqrt(1/cROAD_SM3_star))
  (edrROAD_ZM_star <- sqrt(1/cROAD_ZM_star))
  (edrFORESTC_HUM_star <- sqrt(1/cFORESTC_HUM_star))
  (edrFORESTC_RF_star <- sqrt(1/cFORESTC_RF_star))
  (edrFORESTC_SM2_star <- sqrt(1/cFORESTC_SM2_star))
  (edrFORESTC_SM3_star <- sqrt(1/cFORESTC_SM3_star))
  (edrFORESTC_ZM_star <- sqrt(1/cFORESTC_ZM_star))
  (edrFORESTD_HUM_star <- sqrt(1/cFORESTD_HUM_star))
  (edrFORESTD_RF_star <- sqrt(1/cFORESTD_RF_star))
  (edrFORESTD_SM2_star <- sqrt(1/cFORESTD_SM2_star))
  (edrFORESTD_SM3_star <- sqrt(1/cFORESTD_SM3_star))
  (edrFORESTD_ZM_star <- sqrt(1/cFORESTD_ZM_star))
  
  ROAD_HUM[[i]] <- edrROAD_HUM_star
  ROAD_RF[[i]] <- edrROAD_RF_star
  ROAD_SM2[[i]] <- edrROAD_SM2_star
  ROAD_SM3[[i]] <- edrROAD_SM3_star
  ROAD_ZM[[i]] <- edrROAD_ZM_star
  FORESTC_HUM[[i]] <- edrFORESTC_HUM_star
  FORESTC_RF[[i]] <- edrFORESTC_RF_star
  FORESTC_SM2[[i]] <- edrFORESTC_SM2_star
  FORESTC_SM3[[i]] <- edrFORESTC_SM3_star
  FORESTC_ZM[[i]] <- edrFORESTC_ZM_star
  FORESTD_HUM[[i]] <- edrFORESTD_HUM_star
  FORESTD_RF[[i]] <- edrFORESTD_RF_star
  FORESTD_SM2[[i]] <- edrFORESTD_SM2_star
  FORESTD_SM3[[i]] <- edrFORESTD_SM3_star
  FORESTD_ZM[[i]] <- edrFORESTD_ZM_star
}

ROAD_HUM_mat <- do.call(rbind, ROAD_HUM)
ROAD_RF_mat <- do.call(rbind, ROAD_RF)
ROAD_SM2_mat <- do.call(rbind, ROAD_SM2)
ROAD_SM3_mat <- do.call(rbind, ROAD_SM3)
ROAD_ZM_mat <- do.call(rbind, ROAD_ZM)
FORESTC_HUM_mat <- do.call(rbind, FORESTC_HUM)
FORESTC_RF_mat <- do.call(rbind, FORESTC_RF)
FORESTC_SM2_mat <- do.call(rbind, FORESTC_SM2)
FORESTC_SM3_mat <- do.call(rbind, FORESTC_SM3)
FORESTC_ZM_mat <- do.call(rbind, FORESTC_ZM)
FORESTD_HUM_mat <- do.call(rbind, FORESTD_HUM)
FORESTD_RF_mat <- do.call(rbind, FORESTD_RF)
FORESTD_SM2_mat <- do.call(rbind, FORESTD_SM2)
FORESTD_SM3_mat <- do.call(rbind, FORESTD_SM3)
FORESTD_ZM_mat <- do.call(rbind, FORESTD_ZM)

confidence[["T11313Hz.ROAD_HUM"]] <- apply(ROAD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_RF"]] <- apply(ROAD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_SM2"]] <- apply(ROAD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_SM3"]] <- apply(ROAD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.ROAD_ZM"]] <- apply(ROAD_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTC_HUM"]] <- apply(FORESTC_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTC_RF"]] <- apply(FORESTC_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTC_SM2"]] <- apply(FORESTC_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTC_SM3"]] <- apply(FORESTC_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz,FORESTC_ZM"]] <- apply(FORESTC_ZM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTD_HUM"]] <- apply(FORESTD_HUM_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTD_RF"]] <- apply(FORESTD_RF_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTD_SM2"]] <- apply(FORESTD_SM2_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTD_SM3"]] <- apply(FORESTD_SM3_mat, 2, quantile, c(0.05, 0.95))
confidence[["T11313Hz.FORESTD_ZM"]] <- apply(FORESTD_ZM_mat, 2, quantile, c(0.05, 0.95))

####SUMMARIZE CONFIDENCE####
conf1 <- data.frame(confidence)
conf2 <- names(confidence)
names(conf1) <- c(conf2)
conf3 <- unlist(confidence)
conf3 <- t(conf1)
colnames(conf3) <- c("LL90CI", "UL90CI")
conf3 <- data.frame(conf3)
conf3["Interval"] <- abs(conf3$UL90CI - conf3$LL90CI)
write.csv(conf3, file="uncertainty.csv")

####SUMMARIZE EDR####
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Road", Type="RiverForks", EDR=edr[["BEKI.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Road", Type="SM2", EDR=edr[["BEKI.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Road", Type="SM3", EDR=edr[["BEKI.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Road", Type="Zoom", EDR=edr[["BEKI.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Conifer", Type="RiverForks", EDR=edr[["BEKI.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Conifer", Type="SM2", EDR=edr[["BEKI.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Conifer", Type="SM3", EDR=edr[["BEKI.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Conifer", Type="Zoom", EDR=edr[["BEKI.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BEKI.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Deciduous", Type="SM2", EDR=edr[["BEKI.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Deciduous", Type="SM3", EDR=edr[["BEKI.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Deciduous", Type="Zoom", EDR=edr[["BEKI.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Road", Type="RiverForks", EDR=edr[["PISI.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Road", Type="SM2", EDR=edr[["PISI.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Road", Type="SM3", EDR=edr[["PISI.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Road", Type="Zoom", EDR=edr[["PISI.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Conifer", Type="RiverForks", EDR=edr[["PISI.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Conifer", Type="SM2", EDR=edr[["PISI.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Conifer", Type="SM3", EDR=edr[["PISI.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Conifer", Type="Zoom", EDR=edr[["PISI.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Deciduous", Type="RiverForks", EDR=edr[["PISI.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Deciduous", Type="SM2", EDR=edr[["PISI.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Deciduous", Type="SM3", EDR=edr[["PISI.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Deciduous", Type="Zoom", EDR=edr[["PISI.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Road", Type="RiverForks", EDR=edr[["TEWA.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Road", Type="SM2", EDR=edr[["TEWA.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Road", Type="SM3", EDR=edr[["TEWA.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Road", Type="Zoom", EDR=edr[["TEWA.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Conifer", Type="RiverForks", EDR=edr[["TEWA.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Conifer", Type="SM2", EDR=edr[["TEWA.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Conifer", Type="SM3", EDR=edr[["TEWA.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Conifer", Type="Zoom", EDR=edr[["TEWA.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Deciduous", Type="RiverForks", EDR=edr[["TEWA.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Deciduous", Type="SM2", EDR=edr[["TEWA.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Deciduous", Type="SM3", EDR=edr[["TEWA.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="TEWA", Habitat="Deciduous", Type="Zoom", EDR=edr[["TEWA.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Road", Type="RiverForks", EDR=edr[["BLWA.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Road", Type="SM2", EDR=edr[["BLWA.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Road", Type="SM3", EDR=edr[["BLWA.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Road", Type="Zoom", EDR=edr[["BLWA.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Conifer", Type="RiverForks", EDR=edr[["BLWA.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Conifer", Type="SM2", EDR=edr[["BLWA.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Conifer", Type="SM3", EDR=edr[["BLWA.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Conifer", Type="Zoom", EDR=edr[["BLWA.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BLWA.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Deciduous", Type="SM2", EDR=edr[["BLWA.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Deciduous", Type="SM3", EDR=edr[["BLWA.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Deciduous", Type="Zoom", EDR=edr[["BLWA.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Road", Type="RiverForks", EDR=edr[["LISP.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Road", Type="SM2", EDR=edr[["LISP.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Road", Type="SM3", EDR=edr[["LISP.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Road", Type="Zoom", EDR=edr[["LISP.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Conifer", Type="RiverForks", EDR=edr[["LISP.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Conifer", Type="SM2", EDR=edr[["LISP.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Conifer", Type="SM3", EDR=edr[["LISP.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Conifer", Type="Zoom", EDR=edr[["LISP.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Deciduous", Type="RiverForks", EDR=edr[["LISP.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Deciduous", Type="SM2", EDR=edr[["LISP.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Deciduous", Type="SM3", EDR=edr[["LISP.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Deciduous", Type="Zoom", EDR=edr[["LISP.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Road", Type="RiverForks", EDR=edr[["CCSP.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Road", Type="SM2", EDR=edr[["CCSP.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Road", Type="SM3", EDR=edr[["CCSP.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Road", Type="Zoom", EDR=edr[["CCSP.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Conifer", Type="RiverForks", EDR=edr[["CCSP.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Conifer", Type="SM2", EDR=edr[["CCSP.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Conifer", Type="SM3", EDR=edr[["CCSP.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Conifer", Type="Zoom", EDR=edr[["CCSP.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Deciduous", Type="RiverForks", EDR=edr[["CCSP.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Deciduous", Type="SM2", EDR=edr[["CCSP.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Deciduous", Type="SM3", EDR=edr[["CCSP.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Deciduous", Type="Zoom", EDR=edr[["CCSP.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Road", Type="RiverForks", EDR=edr[["BAWW.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Road", Type="SM2", EDR=edr[["BAWW.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Road", Type="SM3", EDR=edr[["BAWW.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Road", Type="Zoom", EDR=edr[["BAWW.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Conifer", Type="RiverForks", EDR=edr[["BAWW.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Conifer", Type="SM2", EDR=edr[["BAWW.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Conifer", Type="SM3", EDR=edr[["BAWW.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Conifer", Type="Zoom", EDR=edr[["BAWW.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BAWW.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Deciduous", Type="SM2", EDR=edr[["BAWW.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Deciduous", Type="SM3", EDR=edr[["BAWW.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Deciduous", Type="Zoom", EDR=edr[["BAWW.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Road", Type="RiverForks", EDR=edr[["WAVI.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Road", Type="SM2", EDR=edr[["WAVI.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Road", Type="SM3", EDR=edr[["WAVI.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Road", Type="Zoom", EDR=edr[["WAVI.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Conifer", Type="RiverForks", EDR=edr[["WAVI.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Conifer", Type="SM2", EDR=edr[["WAVI.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Conifer", Type="SM3", EDR=edr[["WAVI.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Conifer", Type="Zoom", EDR=edr[["WAVI.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Deciduous", Type="RiverForks", EDR=edr[["WAVI.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Deciduous", Type="SM2", EDR=edr[["WAVI.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Deciduous", Type="SM3", EDR=edr[["WAVI.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Deciduous", Type="Zoom", EDR=edr[["WAVI.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Road", Type="RiverForks", EDR=edr[["OVEN.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Road", Type="SM2", EDR=edr[["OVEN.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Road", Type="SM3", EDR=edr[["OVEN.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Road", Type="Zoom", EDR=edr[["OVEN.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Conifer", Type="RiverForks", EDR=edr[["OVEN.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Conifer", Type="SM2", EDR=edr[["OVEN.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Conifer", Type="SM3", EDR=edr[["OVEN.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Conifer", Type="Zoom", EDR=edr[["OVEN.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Deciduous", Type="RiverForks", EDR=edr[["OVEN.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Deciduous", Type="SM2", EDR=edr[["OVEN.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Deciduous", Type="SM3", EDR=edr[["OVEN.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Deciduous", Type="Zoom", EDR=edr[["OVEN.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Road", Type="RiverForks", EDR=edr[["DEJU.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Road", Type="SM2", EDR=edr[["DEJU.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Road", Type="SM3", EDR=edr[["DEJU.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Road", Type="Zoom", EDR=edr[["DEJU.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Conifer", Type="RiverForks", EDR=edr[["DEJU.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Conifer", Type="SM2", EDR=edr[["DEJU.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Conifer", Type="SM3", EDR=edr[["DEJU.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Conifer", Type="Zoom", EDR=edr[["DEJU.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Deciduous", Type="RiverForks", EDR=edr[["DEJU.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Deciduous", Type="SM2", EDR=edr[["DEJU.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Deciduous", Type="SM3", EDR=edr[["DEJU.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Deciduous", Type="Zoom", EDR=edr[["DEJU.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Road", Type="RiverForks", EDR=edr[["BHCO.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Road", Type="SM2", EDR=edr[["BHCO.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Road", Type="SM3", EDR=edr[["BHCO.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Road", Type="Zoom", EDR=edr[["BHCO.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Conifer", Type="RiverForks", EDR=edr[["BHCO.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Conifer", Type="SM2", EDR=edr[["BHCO.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Conifer", Type="SM3", EDR=edr[["BHCO.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Conifer", Type="Zoom", EDR=edr[["BHCO.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BHCO.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Deciduous", Type="SM2", EDR=edr[["BHCO.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Deciduous", Type="SM3", EDR=edr[["BHCO.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Deciduous", Type="Zoom", EDR=edr[["BHCO.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Road", Type="RiverForks", EDR=edr[["4000Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Road", Type="SM2", EDR=edr[["4000Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Road", Type="SM3", EDR=edr[["4000Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Road", Type="Zoom", EDR=edr[["4000Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["4000Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Conifer", Type="SM2", EDR=edr[["4000Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Conifer", Type="SM3", EDR=edr[["4000Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["4000Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["4000Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["4000Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["4000Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["4000Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Road", Type="RiverForks", EDR=edr[["2828Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Road", Type="SM2", EDR=edr[["2828Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Road", Type="SM3", EDR=edr[["2828Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Road", Type="Zoom", EDR=edr[["2828Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["2828Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Conifer", Type="SM2", EDR=edr[["2828Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Conifer", Type="SM3", EDR=edr[["2828Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["2828Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["2828Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["2828Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["2828Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["2828Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Road", Type="RiverForks", EDR=edr[["5656Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Road", Type="SM2", EDR=edr[["5656Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Road", Type="SM3", EDR=edr[["5656Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Road", Type="Zoom", EDR=edr[["5656Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["5656Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Conifer", Type="SM2", EDR=edr[["5656Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Conifer", Type="SM3", EDR=edr[["5656Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["5656Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["5656Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["5656Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["5656Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["5656Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Road", Type="RiverForks", EDR=edr[["8000Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Road", Type="SM2", EDR=edr[["8000Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Road", Type="SM3", EDR=edr[["8000Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Road", Type="Zoom", EDR=edr[["8000Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["8000Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Conifer", Type="SM2", EDR=edr[["8000Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Conifer", Type="SM3", EDR=edr[["8000Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["8000Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["8000Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["8000Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["8000Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["8000Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Road", Type="RiverForks", EDR=edr[["BOOW.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Road", Type="SM2", EDR=edr[["BOOW.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Road", Type="SM3", EDR=edr[["BOOW.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Road", Type="Zoom", EDR=edr[["BOOW.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Conifer", Type="RiverForks", EDR=edr[["BOOW.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Conifer", Type="SM2", EDR=edr[["BOOW.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Conifer", Type="SM3", EDR=edr[["BOOW.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Conifer", Type="Zoom", EDR=edr[["BOOW.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BOOW.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Deciduous", Type="SM2", EDR=edr[["BOOW.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Deciduous", Type="SM3", EDR=edr[["BOOW.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Deciduous", Type="Zoom", EDR=edr[["BOOW.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Road", Type="RiverForks", EDR=edr[["NSWO.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Road", Type="SM2", EDR=edr[["NSWO.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Road", Type="SM3", EDR=edr[["NSWO.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Road", Type="Zoom", EDR=edr[["NSWO.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Conifer", Type="RiverForks", EDR=edr[["NSWO.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Conifer", Type="SM2", EDR=edr[["NSWO.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Conifer", Type="SM3", EDR=edr[["NSWO.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Conifer", Type="Zoom", EDR=edr[["NSWO.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Deciduous", Type="RiverForks", EDR=edr[["NSWO.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Deciduous", Type="SM2", EDR=edr[["NSWO.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Deciduous", Type="SM3", EDR=edr[["NSWO.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Deciduous", Type="Zoom", EDR=edr[["NSWO.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Road", Type="RiverForks", EDR=edr[["BBWA.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Road", Type="SM2", EDR=edr[["BBWA.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Road", Type="SM3", EDR=edr[["BBWA.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Road", Type="Zoom", EDR=edr[["BBWA.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Conifer", Type="RiverForks", EDR=edr[["BBWA.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Conifer", Type="SM2", EDR=edr[["BBWA.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Conifer", Type="SM3", EDR=edr[["BBWA.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Conifer", Type="Zoom", EDR=edr[["BBWA.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BBWA.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Deciduous", Type="SM2", EDR=edr[["BBWA.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Deciduous", Type="SM3", EDR=edr[["BBWA.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Deciduous", Type="Zoom", EDR=edr[["BBWA.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Road", Type="RiverForks", EDR=edr[["LEOW.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Road", Type="SM2", EDR=edr[["LEOW.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Road", Type="SM3", EDR=edr[["LEOW.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Road", Type="Zoom", EDR=edr[["LEOW.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Conifer", Type="RiverForks", EDR=edr[["LEOW.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Conifer", Type="SM2", EDR=edr[["LEOW.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Conifer", Type="SM3", EDR=edr[["LEOW.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Conifer", Type="Zoom", EDR=edr[["LEOW.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Deciduous", Type="RiverForks", EDR=edr[["LEOW.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Deciduous", Type="SM2", EDR=edr[["LEOW.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Deciduous", Type="SM3", EDR=edr[["LEOW.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Deciduous", Type="Zoom", EDR=edr[["LEOW.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Road", Type="RiverForks", EDR=edr[["BADO.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Road", Type="SM2", EDR=edr[["BADO.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Road", Type="SM3", EDR=edr[["BADO.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Road", Type="Zoom", EDR=edr[["BADO.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Conifer", Type="RiverForks", EDR=edr[["BADO.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Conifer", Type="SM2", EDR=edr[["BADO.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Conifer", Type="SM3", EDR=edr[["BADO.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Conifer", Type="Zoom", EDR=edr[["BADO.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Deciduous", Type="RiverForks", EDR=edr[["BADO.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Deciduous", Type="SM2", EDR=edr[["BADO.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Deciduous", Type="SM3", EDR=edr[["BADO.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Deciduous", Type="Zoom", EDR=edr[["BADO.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Road", Type="RiverForks", EDR=edr[["WETO.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Road", Type="SM2", EDR=edr[["WETO.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Road", Type="SM3", EDR=edr[["WETO.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Road", Type="Zoom", EDR=edr[["WETO.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Conifer", Type="RiverForks", EDR=edr[["WETO.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Conifer", Type="SM2", EDR=edr[["WETO.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Conifer", Type="SM3", EDR=edr[["WETO.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Conifer", Type="Zoom", EDR=edr[["WETO.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Deciduous", Type="RiverForks", EDR=edr[["WETO.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Deciduous", Type="SM2", EDR=edr[["WETO.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Deciduous", Type="SM3", EDR=edr[["WETO.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Deciduous", Type="Zoom", EDR=edr[["WETO.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Road", Type="RiverForks", EDR=edr[["CORA.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Road", Type="SM2", EDR=edr[["CORA.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Road", Type="SM3", EDR=edr[["CORA.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Road", Type="Zoom", EDR=edr[["CORA.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Conifer", Type="RiverForks", EDR=edr[["CORA.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Conifer", Type="SM2", EDR=edr[["CORA.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Conifer", Type="SM3", EDR=edr[["CORA.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Conifer", Type="Zoom", EDR=edr[["CORA.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Deciduous", Type="RiverForks", EDR=edr[["CORA.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Deciduous", Type="SM2", EDR=edr[["CORA.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Deciduous", Type="SM3", EDR=edr[["CORA.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Deciduous", Type="Zoom", EDR=edr[["CORA.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Road", Type="RiverForks", EDR=edr[["RBNU.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Road", Type="SM2", EDR=edr[["RBNU.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Road", Type="SM3", EDR=edr[["RBNU.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Road", Type="Zoom", EDR=edr[["RBNU.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Conifer", Type="RiverForks", EDR=edr[["RBNU.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Conifer", Type="SM2", EDR=edr[["RBNU.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Conifer", Type="SM3", EDR=edr[["RBNU.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Conifer", Type="Zoom", EDR=edr[["RBNU.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Deciduous", Type="RiverForks", EDR=edr[["RBNU.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Deciduous", Type="SM2", EDR=edr[["RBNU.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Deciduous", Type="SM3", EDR=edr[["RBNU.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Deciduous", Type="Zoom", EDR=edr[["RBNU.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Road", Type="RiverForks", EDR=edr[["GGOW.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Road", Type="SM2", EDR=edr[["GGOW.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Road", Type="SM3", EDR=edr[["GGOW.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Road", Type="Zoom", EDR=edr[["GGOW.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Conifer", Type="RiverForks", EDR=edr[["GGOW.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Conifer", Type="SM2", EDR=edr[["GGOW.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Conifer", Type="SM3", EDR=edr[["GGOW.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Conifer", Type="Zoom", EDR=edr[["GGOW.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Deciduous", Type="RiverForks", EDR=edr[["GGOW.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Deciduous", Type="SM2", EDR=edr[["GGOW.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Deciduous", Type="SM3", EDR=edr[["GGOW.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Deciduous", Type="Zoom", EDR=edr[["GGOW.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Road", Type="RiverForks", EDR=edr[["WTSP.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Road", Type="SM2", EDR=edr[["WTSP.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Road", Type="SM3", EDR=edr[["WTSP.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Road", Type="Zoom", EDR=edr[["WTSP.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Conifer", Type="RiverForks", EDR=edr[["WTSP.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Conifer", Type="SM2", EDR=edr[["WTSP.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Conifer", Type="SM3", EDR=edr[["WTSP.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Conifer", Type="Zoom", EDR=edr[["WTSP.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Deciduous", Type="RiverForks", EDR=edr[["WTSP.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Deciduous", Type="SM2", EDR=edr[["WTSP.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Deciduous", Type="SM3", EDR=edr[["WTSP.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Deciduous", Type="Zoom", EDR=edr[["WTSP.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Road", Type="RiverForks", EDR=edr[["OSFL.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Road", Type="SM2", EDR=edr[["OSFL.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Road", Type="SM3", EDR=edr[["OSFL.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Road", Type="Zoom", EDR=edr[["OSFL.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Conifer", Type="RiverForks", EDR=edr[["OSFL.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Conifer", Type="SM2", EDR=edr[["OSFL.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Conifer", Type="SM3", EDR=edr[["OSFL.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Conifer", Type="Zoom", EDR=edr[["OSFL.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Deciduous", Type="RiverForks", EDR=edr[["OSFL.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Deciduous", Type="SM2", EDR=edr[["OSFL.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Deciduous", Type="SM3", EDR=edr[["OSFL.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Deciduous", Type="Zoom", EDR=edr[["OSFL.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Road", Type="RiverForks", EDR=edr[["RBGR.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Road", Type="SM2", EDR=edr[["RBGR.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Road", Type="SM3", EDR=edr[["RBGR.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Road", Type="Zoom", EDR=edr[["RBGR.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Conifer", Type="RiverForks", EDR=edr[["RBGR.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Conifer", Type="SM2", EDR=edr[["RBGR.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Conifer", Type="SM3", EDR=edr[["RBGR.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Conifer", Type="Zoom", EDR=edr[["RBGR.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Deciduous", Type="RiverForks", EDR=edr[["RBGR.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Deciduous", Type="SM2", EDR=edr[["RBGR.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Deciduous", Type="SM3", EDR=edr[["RBGR.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Deciduous", Type="Zoom", EDR=edr[["RBGR.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Road", Type="RiverForks", EDR=edr[["CATO.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Road", Type="SM2", EDR=edr[["CATO.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Road", Type="SM3", EDR=edr[["CATO.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Road", Type="Zoom", EDR=edr[["CATO.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Conifer", Type="RiverForks", EDR=edr[["CATO.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Conifer", Type="SM2", EDR=edr[["CATO.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Conifer", Type="SM3", EDR=edr[["CATO.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Conifer", Type="Zoom", EDR=edr[["CATO.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Deciduous", Type="RiverForks", EDR=edr[["CATO.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Deciduous", Type="SM2", EDR=edr[["CATO.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Deciduous", Type="SM3", EDR=edr[["CATO.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Deciduous", Type="Zoom", EDR=edr[["CATO.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Road", Type="RiverForks", EDR=edr[["2000Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Road", Type="SM2", EDR=edr[["2000Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Road", Type="SM3", EDR=edr[["2000Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Road", Type="Zoom", EDR=edr[["2000Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["2000Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Conifer", Type="SM2", EDR=edr[["2000Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Conifer", Type="SM3", EDR=edr[["2000Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["2000Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["2000Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["2000Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["2000Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["2000Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Road", Type="RiverForks", EDR=edr[["1000Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Road", Type="SM2", EDR=edr[["1000Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Road", Type="SM3", EDR=edr[["1000Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Road", Type="Zoom", EDR=edr[["1000Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["1000Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Conifer", Type="SM2", EDR=edr[["1000Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Conifer", Type="SM3", EDR=edr[["1000Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["1000Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["1000Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["1000Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["1000Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["1000Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Road", Type="RiverForks", EDR=edr[["1414Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Road", Type="SM2", EDR=edr[["1414Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Road", Type="SM3", EDR=edr[["1414Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Road", Type="Zoom", EDR=edr[["1414Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["1414Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Conifer", Type="SM2", EDR=edr[["1414Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Conifer", Type="SM3", EDR=edr[["1414Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["1414Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["1414Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["1414Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["1414Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["1414Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Road", Type="RiverForks", EDR=edr[["YERA.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Road", Type="SM2", EDR=edr[["YERA.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Road", Type="SM3", EDR=edr[["YERA.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Road", Type="Zoom", EDR=edr[["YERA.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Conifer", Type="RiverForks", EDR=edr[["YERA.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Conifer", Type="SM2", EDR=edr[["YERA.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Conifer", Type="SM3", EDR=edr[["YERA.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Conifer", Type="Zoom", EDR=edr[["YERA.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Deciduous", Type="RiverForks", EDR=edr[["YERA.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Deciduous", Type="SM2", EDR=edr[["YERA.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Deciduous", Type="SM3", EDR=edr[["YERA.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Deciduous", Type="Zoom", EDR=edr[["YERA.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Road", Type="RiverForks", EDR=edr[["11313Hz.ROAD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Road", Type="SM2", EDR=edr[["11313Hz.ROAD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Road", Type="SM3", EDR=edr[["11313Hz.ROAD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Road", Type="Zoom", EDR=edr[["11313Hz.ROAD_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Conifer", Type="RiverForks", EDR=edr[["11313Hz.FORESTC_RF"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Conifer", Type="SM2", EDR=edr[["11313Hz.FORESTC_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Conifer", Type="SM3", EDR=edr[["11313Hz.FORESTC_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Conifer", Type="Zoom", EDR=edr[["11313Hz.FORESTC_ZM"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Deciduous", Type="RiverForks", EDR=edr[["11313Hz.FORESTD_RF"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Deciduous", Type="SM2", EDR=edr[["11313Hz.FORESTD_SM2"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Deciduous", Type="SM3", EDR=edr[["11313Hz.FORESTD_SM3"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Deciduous", Type="Zoom", EDR=edr[["11313Hz.FORESTD_ZM"]]))

df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Road", Type="Human", EDR=edr[["BEKI.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Road", Type="Human", EDR=edr[["PISI.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="TEwA", Habitat="Road", Type="Human", EDR=edr[["TEWA.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Road", Type="Human", EDR=edr[["BLWA.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Road", Type="Human", EDR=edr[["LISP.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Road", Type="Human", EDR=edr[["CCSP.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Road", Type="Human", EDR=edr[["BAWW.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Road", Type="Human", EDR=edr[["WAVI.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Road", Type="Human", EDR=edr[["OVEN.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Road", Type="Human", EDR=edr[["DEJU.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Road", Type="Human", EDR=edr[["BHCO.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Road", Type="Human", EDR=edr[["4000Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Road", Type="Human", EDR=edr[["2828Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Road", Type="Human", EDR=edr[["5656Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Road", Type="Human", EDR=edr[["8000Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Road", Type="Human", EDR=edr[["BOOW.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Road", Type="Human", EDR=edr[["NSWO.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Road", Type="Human", EDR=edr[["BBWA.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Road", Type="Human", EDR=edr[["LEOW.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Road", Type="Human", EDR=edr[["BADO.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Road", Type="Human", EDR=edr[["WETO.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Road", Type="Human", EDR=edr[["CORA.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Road", Type="Human", EDR=edr[["RBNU.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Road", Type="Human", EDR=edr[["GGOW.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Road", Type="Human", EDR=edr[["WTSP.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Road", Type="Human", EDR=edr[["OSFL.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Road", Type="Human", EDR=edr[["RBGR.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Road", Type="Human", EDR=edr[["CATO.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Road", Type="Human", EDR=edr[["2000Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Road", Type="Human", EDR=edr[["1000Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Road", Type="Human", EDR=edr[["1414Hz.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Road", Type="Human", EDR=edr[["YERA.ROAD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Road", Type="Human", EDR=edr[["11313Hz.ROAD_HUM"]]))

df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Conifer", Type="Human", EDR=edr[["BEKI.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Conifer", Type="Human", EDR=edr[["PISI.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="TEwA", Habitat="Conifer", Type="Human", EDR=edr[["TEWA.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Conifer", Type="Human", EDR=edr[["BLWA.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Conifer", Type="Human", EDR=edr[["LISP.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Conifer", Type="Human", EDR=edr[["CCSP.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Conifer", Type="Human", EDR=edr[["BAWW.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Conifer", Type="Human", EDR=edr[["WAVI.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Conifer", Type="Human", EDR=edr[["OVEN.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Conifer", Type="Human", EDR=edr[["DEJU.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Conifer", Type="Human", EDR=edr[["BHCO.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Conifer", Type="Human", EDR=edr[["4000Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Conifer", Type="Human", EDR=edr[["2828Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Conifer", Type="Human", EDR=edr[["5656Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Conifer", Type="Human", EDR=edr[["8000Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Conifer", Type="Human", EDR=edr[["BOOW.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Conifer", Type="Human", EDR=edr[["NSWO.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Conifer", Type="Human", EDR=edr[["BBWA.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Conifer", Type="Human", EDR=edr[["LEOW.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Conifer", Type="Human", EDR=edr[["BADO.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Conifer", Type="Human", EDR=edr[["WETO.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Conifer", Type="Human", EDR=edr[["CORA.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Conifer", Type="Human", EDR=edr[["RBNU.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Conifer", Type="Human", EDR=edr[["GGOW.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Conifer", Type="Human", EDR=edr[["WTSP.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Conifer", Type="Human", EDR=edr[["OSFL.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Conifer", Type="Human", EDR=edr[["RBGR.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Conifer", Type="Human", EDR=edr[["CATO.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Conifer", Type="Human", EDR=edr[["2000Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Conifer", Type="Human", EDR=edr[["1000Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Conifer", Type="Human", EDR=edr[["1414Hz.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Conifer", Type="Human", EDR=edr[["YERA.FORESTC_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Conifer", Type="Human", EDR=edr[["11313Hz.FORESTC_HUM"]]))

df1 <- rbind(df1, data.frame(Sound="BEKI", Habitat="Deciduous", Type="Human", EDR=edr[["BEKI.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="PISI", Habitat="Deciduous", Type="Human", EDR=edr[["PISI.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="TEwA", Habitat="Deciduous", Type="Human", EDR=edr[["TEWA.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BLWA", Habitat="Deciduous", Type="Human", EDR=edr[["BLWA.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LISP", Habitat="Deciduous", Type="Human", EDR=edr[["LISP.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CCSP", Habitat="Deciduous", Type="Human", EDR=edr[["CCSP.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BAWW", Habitat="Deciduous", Type="Human", EDR=edr[["BAWW.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WAVI", Habitat="Deciduous", Type="Human", EDR=edr[["WAVI.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OVEN", Habitat="Deciduous", Type="Human", EDR=edr[["OVEN.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="DEJU", Habitat="Deciduous", Type="Human", EDR=edr[["DEJU.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BHCO", Habitat="Deciduous", Type="Human", EDR=edr[["BHCO.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="4000Hz", Habitat="Deciduous", Type="Human", EDR=edr[["4000Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2828Hz", Habitat="Deciduous", Type="Human", EDR=edr[["2828Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="5656Hz", Habitat="Deciduous", Type="Human", EDR=edr[["5656Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="8000Hz", Habitat="Deciduous", Type="Human", EDR=edr[["8000Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BOOW", Habitat="Deciduous", Type="Human", EDR=edr[["BOOW.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="NSWO", Habitat="Deciduous", Type="Human", EDR=edr[["NSWO.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BBWA", Habitat="Deciduous", Type="Human", EDR=edr[["BBWA.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="LEOW", Habitat="Deciduous", Type="Human", EDR=edr[["LEOW.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="BADO", Habitat="Deciduous", Type="Human", EDR=edr[["BADO.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WETO", Habitat="Deciduous", Type="Human", EDR=edr[["WETO.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CORA", Habitat="Deciduous", Type="Human", EDR=edr[["CORA.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBNU", Habitat="Deciduous", Type="Human", EDR=edr[["RBNU.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="GGOW", Habitat="Deciduous", Type="Human", EDR=edr[["GGOW.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="WTSP", Habitat="Deciduous", Type="Human", EDR=edr[["WTSP.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="OSFL", Habitat="Deciduous", Type="Human", EDR=edr[["OSFL.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="RBGR", Habitat="Deciduous", Type="Human", EDR=edr[["RBGR.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="CATO", Habitat="Deciduous", Type="Human", EDR=edr[["CATO.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="2000Hz", Habitat="Deciduous", Type="Human", EDR=edr[["2000Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1000Hz", Habitat="Deciduous", Type="Human", EDR=edr[["1000Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="1414Hz", Habitat="Deciduous", Type="Human", EDR=edr[["1414Hz.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="YERA", Habitat="Deciduous", Type="Human", EDR=edr[["YERA.FORESTD_HUM"]]))
df1 <- rbind(df1, data.frame(Sound="11313Hz", Habitat="Deciduous", Type="Human", EDR=edr[["11313Hz.FORESTD_HUM"]]))

write.csv(df1, file = "EDR.csv")

mod.list <-names(model.rank)
birdframe.list <- list()
MASTERBIRDFRAME <- NSWO.1
for(i in 1:length(mod.list)){
  BIRDFRAME <- as.data.frame(model.rank[[mod.list[i]]])[9:13]
  BIRDFRAME$sound <- c(mod.list[i])
  BIRDFRAME$model <- c(row.names(BIRDFRAME))
  MASTERBIRDFRAME <- rbind(MASTERBIRDFRAME, BIRDFRAME)
  #birdframe.list[[mod.list[i]]] <- assign(paste(mod.list[i], ".1", sep=""), BIRDFRAME)
  #birdframe.list[[mod.list[i]]] <- BIRDFRAME
}
setwd("C:/Users/Daniel/Desktop/Analysis/Human vs ARU/GLM/model.sel/")
write.csv(MASTERBIRDFRAME, file="model.sel.csv")
