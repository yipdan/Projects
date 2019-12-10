library(dplyr)
library(psych)
library(lme4)
library(ICC)
library(MuMIn)
library(car)
library(yhat)
library(sjPlot)
library(sjmisc)
library(ggplot2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(gridExtra)
options(max.print = 10000)

setwd("C:/Users/Daniel/Desktop/Analysis/Multiple Scales of Habitat/")

habitat <- read.csv(file = "habitat_detections_cleaned.csv", header = TRUE)
rf <- read.csv(file = "roadforest_detections_cleaned.csv", header = TRUE)
meta <- read.csv(file = "species_parameters.csv", header = TRUE)
str(habitat)
str(rf)
str(meta)

hbt <- rbind(habitat,rf)
hbt$Observer <- as.factor(hbt$Observer)

levels(hbt$Species)["BAOW"] <- "BADO"
hbt$Species[hbt$Species == "BAOW"]  <- "BADO"
levels(hbt$Species)["BLWA"] <- "CMWA"
hbt$Species[hbt$Species == "BLWA"] <- "CMWA"

hbt <- hbt[hbt$Species != "11313Hz", ]
hbt <- hbt[hbt$Species != "8000Hz", ]

hbt <- merge(hbt, meta, by.x = "Species", by.y = "Sound")

hbt$Species <- factor(hbt$Species)

hbt %>%
  group_by(Species) %>%
  summarise(no_rows = length(Species))

hbt %>%
  group_by(Habitat) %>%
  summarise(no_rows = length(Habitat))

hbt$Habitat <- as.character(hbt$Habitat)
hbt$Habitat[hbt$Habitat == "B"] <- "Bog"
hbt$Habitat[hbt$Habitat == "D"] <- "Deciduous"
hbt$Habitat[hbt$Habitat == "G"] <- "Grassland"
hbt$Habitat[hbt$Habitat == "R"] <- "Roadside"
hbt$Habitat[hbt$Habitat == "C"] <- "Conifer"
hbt$Habitat[hbt$Habitat == "F"] <- "Fen"
hbt$Habitat[hbt$Habitat == "E"] <- "Forest Edge"
hbt$Habitat <- as.factor(hbt$Habitat)

describe(hbt$Wind)
describe(hbt$Temp)
describe(hbt$Humidity)

hist(hbt$Wind)
hist(hbt$Humidity)
hist(hbt$Temp)

ICCest(x = Observer, y = Detected, data = hbt)
ICCest(x = Transect, y = Detected, data = hbt)
ICCest(x = Species, y = Detected, data = hbt)
#no random effects required; 
#Observer ICC - 0.008012819
#Transect ICC - 0.08220165

#test variable inflation factor
#multicollinearity not present
viftest <- glm(Detected ~ Distance + Habitat + Temp + Wind + Humidity + MinFreq + Bandwidth, data=hbt, family = binomial(link = logit))
vifout <- vif(viftest) #provivde VIF coefficient, drop coefficients greater than 5
durbinWatsonTest(viftest) #provide D-W statistic + p-value if reporting
cor(hbt$Wind, hbt$Temp)
cor(hbt$Wind, hbt$Humidity)
cor(hbt$Temp, hbt$Humidity)

#scale + centre
hbt$Distance.st <- scale(hbt$Distance, center=TRUE, scale = TRUE)
hbt$MinFreq.st <- scale(hbt$MinFreq, center=TRUE, scale = TRUE)
hbt$Srate.st <- scale(hbt$Srate, center=TRUE, scale = TRUE)
hbt$Bandwidth.st <- scale(hbt$Bandwidth, center=TRUE, scale = TRUE)
hbt$Wind.st <- scale(hbt$Wind, center=TRUE, scale = TRUE)
hbt$Temp.st <- scale(hbt$Temp, center=TRUE, scale = TRUE)
hbt$Humidity.st <- scale(hbt$Humidity, center=TRUE, scale=TRUE)

#pooled logistic models with sound parameters
p1 <- glmer(Detected ~ Distance.st*Habitat + (1|Species), data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p2 <- glmer(Detected ~ Distance.st*Habitat + Srate.st + (1|Species), data = hbt, family = binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p3 <- glmer(Detected ~ Distance.st*Habitat + MinFreq.st + (1|Species), data = hbt, family = binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p4 <- glmer(Detected ~ Distance.st*Habitat + Bandwidth.st + (1|Species), data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p5 <- glmer(Detected ~ Distance.st*Habitat + Srate.st + MinFreq.st + (1|Species), data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p6 <- glmer(Detected ~ Distance.st*Habitat + Srate.st + Bandwidth.st + (1|Species), data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p7 <- glmer(Detected ~ Distance.st*Habitat + MinFreq.st + Bandwidth.st + (1|Species), data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p8 <- glmer(Detected ~ Distance.st*Habitat + Srate.st + MinFreq.st + Bandwidth.st + (1|Species), data = hbt, 
            family=binomial(link=logit), control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))

p9 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st + Bandwidth.st + (1|Species), 
            data = hbt, family=binomial(link=logit),
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p10 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st*Habitat + Bandwidth.st + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p11 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st + Bandwidth.st*Habitat + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p12 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st*Habitat + Bandwidth.st*Habitat + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p13 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st*Habitat + Bandwidth.st*Habitat
             + Temp.st + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p14 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st*Habitat + Bandwidth.st*Habitat
             + Humidity.st + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
p15 <- glmer(Detected ~ Distance.st*Habitat + Wind.st*Habitat + MinFreq.st*Habitat + Bandwidth.st*Habitat 
             + Humidity.st + Temp.st + (1|Species), 
             data = hbt, family=binomial(link=logit),
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))


pmodels <- model.sel(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15)
#best model p15

p.apsOut <- aps(hbt, "Detected", list("Distance","Habitat","Wind","MinFreq","Bandwidth","Humidity","Temp"))
p.comF <- data.frame(commonality(p.apsOut))

####pooled plots
s_effects <- glm(Detected ~ MinFreq*Habitat, data = bt, family = binomial(link=logit))

# set_theme(
#   panel.background = element_rect(fill = 'white'),
#   axis.line = element_line(colour = "black", size = 0.3),
#   axis.line.x = element_line(colour = "black", size = 0.3), 
#   axis.text.x = element_text(size=5),
#   axis.line.y = element_line(colour = "black", size = 0.3), 
#   axis.text.y = element_text(size=5),
#   axis.text=element_text(size=5),
#   axis.title=element_text(size=6),
#   axis.text=element_text(colour = "black", size=5),
#   legend.title=element_blank(),
#   plot.title = element_text(hjust=1, margin = margin(t=5, b=-10), size=10, face="bold"),
#   axis.ticks = element_line(colour = "black", size= 0.3),
#   legend.key.size = unit(0.3, "cm"),
#   legend.text = element_text(size=6)
# )
set_theme(
  base = theme_classic(),
  axis.textsize = .1,
  axis.line.size = 0.3,
  legend.inside = TRUE,
  legend.pos = "bottom",
  title.align = "right",
  axis.textcolor = "black"
  
)

freq.plot <- plot_model(p15, type = "pred", terms = c("MinFreq.st", "Habitat"), 
                        show.legend = FALSE, title = "A", 
                        axis.title = c("Standardized Minimum Frequency", "Detection Probability"),
                        axis.lim = c(0,1),
                        line.size = 0.3)
freq.plot <- freq.plot +
  font_size(axis_title.x = 6, axis_title.y = 6, labels.x = 5, labels.y = 5)
  

band.plot <- plot_model(p15, type = "pred", terms = c("Bandwidth.st", "Habitat"),
                        show.legend = FALSE, title = "B",
                        axis.title = c("Standardized Sound Bandwidth", "Detection Probability"),
                        axis.lim = c(0,1),
                        line.size = 0.3)
band.plot <- band.plot +
  font_size(axis_title.x = 6, axis_title.y = 6, labels.x = 5, labels.y = 5)


tiff(width = 3.5, height = 4, units = "in", res = 600, file = "Figure 1 - marginal effects.tiff")
grid.arrange(freq.plot, band.plot, ncol = 1)
dev.off()

test.plot <- plot_model(p15, type = "pred", terms = c("Bandwidth.st", "Habitat"),
                        show.legend = TRUE, legend.title = "Vegetation", title = "A", 
                        axis.title = c("Standardized Minimum Frequency", "Detection Probability"),
                        axis.lim = c(0,1),
                        line.size = 0.3,
                        legend.size = 10)
test.plot <- test.plot +
  font_size(axis_title.x = 6, axis_title.y = 6, labels.x = 5, labels.y = 5)


tiff(width = 7, height = 4, units = "in", res = 600, file = "fig1-legend.tiff")
test.plot
dev.off()

# wind.plot <- plot_model(p15, type = "pred", terms = c("Wind.st", "Habitat"))
# 
# tiff(width=7,height=4, units="in", res=600, file="Figure 1 - FreqPooled.tiff")
# freq.plot
# dev.off()
# 
# tiff(width=7,height=4, units="in", res=600, file="Figure 2 - BandPooled.tiff")
# band.plot
# dev.off()
# 
# tiff(width=7,height=4, units="in", res=600, file="Figure 3 - WindPooled.tiff")
# wind.plot
# dev.off()


####A1: logistic regression; habitat effects####
#Test ICC for all sounds, no random effects required
ICC.List = list()
Sp_List <- unique(hbt$Species)
for(i in 1:length(Sp_List)){
  data <- hbt[hbt$Species == Sp_List[i], ]
  df.ICC <- data.frame()
  ObsICC <- data.frame(ICCest(x = Observer, y = Detected, data = data))
  TransICC <- data.frame(ICCest(x = Transect, y = Detected, data = data))
  ObsICC[["Species"]] <- as.character(paste(Sp_List[i] , sep=""))
  TransICC[["Species"]] <- as.character(paste(Sp_List[i] , sep=""))
  ObsICC[["Effect"]] <- "Obs"
  TransICC[["Effect"]] <- "Transect"
  df.ICC <- rbind(df.ICC, ObsICC)
  df.ICC <- rbind(df.ICC, TransICC)
  ICC.List[[as.character(paste(Sp_List[i] , sep=""))]] = df.ICC
}

ICCout <- do.call("rbind", ICC.List)
#write.csv(ICCout, file = "ICCoutput.csv")

#AIC candidate model selection
AIC.List = list()
for(i in 1:length(Sp_List)){
  data <- hbt[hbt$Species == Sp_List[i], ]
  plogr1 <- glm(Detected ~ Distance, data = data, family = binomial(link = logit))
  plogr2 <- glm(Detected ~ Distance*Habitat, data = data, family = binomial(link = logit))
  plogr3 <- glm(Detected ~ Distance*Habitat + Temp, data = data, family = binomial(link = logit))
  plogr4 <- glm(Detected ~ Distance*Habitat + Humidity, data = data, family = binomial(link = logit))
  plogr5 <- glm(Detected ~ Distance*Habitat + Temp + Humidity, data = data, family = binomial(link = logit))
  plogr6 <- glm(Detected ~ Distance*Habitat + Wind, data = data, family = binomial(link = logit))
  plogr7 <- glm(Detected ~ Distance*Habitat + Wind + Humidity + Temp, data = data, family = binomial(link = logit))
  plogr8 <- glm(Detected ~ Distance*Habitat + Wind + Temp, data = data, family = binomial(link = logit))
  plogr9 <- glm(Detected ~ Distance*Habitat + Wind + Humidity, data = data, family = binomial(link = logit))
  plogr10 <- glm(Detected ~ Distance*Habitat + Wind*Habitat, data = data, family = binomial(link = logit))
  plogr11 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Humidity + Temp, data = data, family = binomial(link = logit))
  plogr12 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Humidity, data = data, family = binomial(link = logit))
  plogr13 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Temp, data = data, family = binomial(link = logit))
  df.modsel <- data.frame(model.sel(plogr1, plogr2, plogr3, plogr4, plogr5, plogr6, plogr7,
                                    plogr8, plogr9, plogr10, plogr11, plogr12, plogr13))
  df.modsel[["Species"]] <- as.character(paste(Sp_List[i] , sep = ""))
  AIC.List[[as.character(paste(Sp_List[i] , sep=""))]] = df.modsel
}
AICout <- do.call("rbind", AIC.List)
#write.csv(AICout, file = "AICoutputs.csv")

top.model <- data.frame("Species" = c("1000Hz", "1414Hz", "2000Hz", "2828Hz", "4000Hz",
                                      "5656Hz", "BADO", "BAWW", "BBWA", "BEKI", "BHCO",
                                      "BOOW", "CATO", "CCSP", "CMWA", "CORA", "DEJU", 
                                      "GGOW", "LEOW", "LISP", "NSWO", "OSFL", "OVEN", 
                                      "PISI", "RBGR", "RBNU", "TEWA", "WAVI", "WETO", 
                                      "WTSP", "YERA"),
                        "Model" = c(13,12,12,13,10,11,13,11,2,10,10,11,8,10,10,12,10,13,11,
                                    2,13,10,4,5,10,12,10,10,13,11,3))
#for LISP, all candidate models very close in dAIC

#generate model coefficients, commonality analysis
temp <- top.model[top.model$Model == 2, ]
model2 <- unique(temp$Species)
temp <- top.model[top.model$Model == 3, ]
model3 <- unique(temp$Species)
temp <- top.model[top.model$Model == 4, ]  
model4 <- unique(temp$Species)
temp <- top.model[top.model$Model == 5, ]  
model5 <- unique(temp$Species)
temp <- top.model[top.model$Model == 8, ]  
model8 <- unique(temp$Species)
temp <- top.model[top.model$Model == 10, ]  
model10 <- unique(temp$Species)
temp <- top.model[top.model$Model == 11, ]  
model11 <- unique(temp$Species)
temp <- top.model[top.model$Model == 12, ]  
model12 <- unique(temp$Species)
temp <- top.model[top.model$Model == 13, ]  
model13 <- unique(temp$Species)

model.list <- list()
model.commonality <- list()
model.coef <- list()

for(i in 1:length(model2)){
  data <- hbt[hbt$Species == model2[i], ]
  plogr2 <- glm(Detected ~ Distance*Habitat, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model2[i] , sep=""))]] = plogr2
  mc <- data.frame(summary(plogr2)$coefficients)
  mc[["Species"]] <- as.character(paste(model2[i] , sep = ""))
  model.coef[[as.character(paste(model2[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model2[i] , sep = ""))
  model.commonality[[as.character(paste(model2[i] , sep=""))]] = comF
}
for(i in 1:length(model3)){
  data <- hbt[hbt$Species == model3[i], ]
  plogr3 <- glm(Detected ~ Distance*Habitat + Temp, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model3[i] , sep=""))]] = plogr3
  mc <- data.frame(summary(plogr3)$coefficients)
  mc[["Species"]] <- as.character(paste(model3[i] , sep = ""))
  model.coef[[as.character(paste(model3[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Temp"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model3[i] , sep = ""))
  model.commonality[[as.character(paste(model3[i] , sep=""))]] = comF
}
for(i in 1:length(model4)){
  data <- hbt[hbt$Species == model4[i], ]
  plogr4 <- glm(Detected ~ Distance*Habitat + Humidity, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model4[i] , sep=""))]] = plogr4
  mc <- data.frame(summary(plogr4)$coefficients)
  mc[["Species"]] <- as.character(paste(model4[i] , sep = ""))
  model.coef[[as.character(paste(model4[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Humidity"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model4[i] , sep = ""))
  model.commonality[[as.character(paste(model4[i] , sep=""))]] = comF
}
for(i in 1:length(model5)){
  data <- hbt[hbt$Species == model5[i], ]
  plogr5 <- glm(Detected ~ Distance*Habitat + Temp + Humidity, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model5[i] , sep=""))]] = plogr5
  mc <- data.frame(summary(plogr5)$coefficients)
  mc[["Species"]] <- as.character(paste(model5[i] , sep = ""))
  model.coef[[as.character(paste(model5[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Temp", "Humidity"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model5[i] , sep = ""))
  model.commonality[[as.character(paste(model5[i] , sep=""))]] = comF
}
for(i in 1:length(model8)){
  data <- hbt[hbt$Species == model8[i], ]
  plogr8 <- glm(Detected ~ Distance*Habitat + Wind + Temp, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model8[i] , sep=""))]] = plogr8
  mc <- data.frame(summary(plogr8)$coefficients)
  mc[["Species"]] <- as.character(paste(model8[i] , sep = ""))
  model.coef[[as.character(paste(model8[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Temp", "Wind"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model8[i] , sep = ""))
  model.commonality[[as.character(paste(model8[i] , sep=""))]] = comF
}
for(i in 1:length(model10)){
  data <- hbt[hbt$Species == model10[i], ]
  plogr10 <- glm(Detected ~ Distance*Habitat + Wind*Habitat, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model10[i] , sep=""))]] = plogr10
  mc <- data.frame(summary(plogr10)$coefficients)
  mc[["Species"]] <- as.character(paste(model10[i] , sep = ""))
  model.coef[[as.character(paste(model10[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Wind"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model10[i] , sep = ""))
  model.commonality[[as.character(paste(model10[i] , sep=""))]] = comF
}
for(i in 1:length(model11)){
  data <- hbt[hbt$Species == model11[i], ]
  plogr11 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Humidity + Temp, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model11[i] , sep=""))]] = plogr11
  mc <- data.frame(summary(plogr11)$coefficients)
  mc[["Species"]] <- as.character(paste(model11[i] , sep = ""))
  model.coef[[as.character(paste(model11[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Humidity", "Temp", "Wind"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model11[i] , sep = ""))
  model.commonality[[as.character(paste(model11[i] , sep=""))]] = comF
}
for(i in 1:length(model12)){
  data <- hbt[hbt$Species == model12[i], ]
  plogr12 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Humidity, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model12[i] , sep=""))]] = plogr12
  mc <- data.frame(summary(plogr12)$coefficients)
  mc[["Species"]] <- as.character(paste(model12[i] , sep = ""))
  model.coef[[as.character(paste(model12[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Humidity", "Wind"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model12[i] , sep = ""))
  model.commonality[[as.character(paste(model12[i] , sep=""))]] = comF
}
for(i in 1:length(model13)){
  data <- hbt[hbt$Species == model13[i], ]
  plogr13 <- glm(Detected ~ Distance*Habitat + Wind*Habitat + Temp, data = data, family = binomial(link = logit))
  model.list[[as.character(paste(model13[i] , sep=""))]] = plogr13
  mc <- data.frame(summary(plogr13)$coefficients)
  mc[["Species"]] <- as.character(paste(model13[i] , sep = ""))
  model.coef[[as.character(paste(model13[i] , sep=""))]] = mc
  
  apsOut <- aps(data, "Detected", list("Distance","Habitat", "Temp", "Wind"))
  comF <- data.frame(commonality(apsOut))
  comF[["Species"]] <- as.character(paste(model13[i] , sep = ""))
  model.commonality[[as.character(paste(model13[i] , sep=""))]] = comF
}

comFout <- do.call("rbind", model.commonality)
#write.csv(comFout, file = "CommonalityOutput.csv")

modelCoef <- do.call("rbind", model.coef)
#write.csv(modelCoef, file = "ModelCoefficients.csv")


W1000 <- plot_model(model.list[["1000Hz"]], type = "pred", terms = c("Wind", "Habitat"))
W1414 <- plot_model(model.list[["1414Hz"]], type = "pred", terms = c("Wind", "Habitat"))
W2000 <- plot_model(model.list[["2000Hz"]], type = "pred", terms = c("Wind", "Habitat"))
W2828 <- plot_model(model.list[["2828Hz"]], type = "pred", terms = c("Wind", "Habitat"))
W4000 <- plot_model(model.list[["4000Hz"]], type = "pred", terms = c("Wind", "Habitat"))
W5656 <- plot_model(model.list[["5656Hz"]], type = "pred", terms = c("Wind", "Habitat"))
WBADO <- plot_model(model.list[["BADO"]], type = "pred", terms = c("Wind", "Habitat"))
WBAWW <- plot_model(model.list[["BAWW"]], type = "pred", terms = c("Wind", "Habitat"))
WBBWA <- plot_model(model.list[["BBWA"]], type = "pred", terms = c("Wind", "Habitat"))
WBEKI <- plot_model(model.list[["BEKI"]], type = "pred", terms = c("Wind", "Habitat"))
WBHCO <- plot_model(model.list[["BHCO"]], type = "pred", terms = c("Wind", "Habitat"))
WBOOW <- plot_model(model.list[["BOOW"]], type = "pred", terms = c("Wind", "Habitat"))
WCATO <- plot_model(model.list[["CATO"]], type = "pred", terms = c("Wind", "Habitat"))
WCCSP <- plot_model(model.list[["CCSP"]], type = "pred", terms = c("Wind", "Habitat"))
WCMWA <- plot_model(model.list[["CMWA"]], type = "pred", terms = c("Wind", "Habitat"))
WCORA <- plot_model(model.list[["CORA"]], type = "pred", terms = c("Wind", "Habitat"))
WDEJU <- plot_model(model.list[["DEJU"]], type = "pred", terms = c("Wind", "Habitat"))
WGGOW <- plot_model(model.list[["GGOW"]], type = "pred", terms = c("Wind", "Habitat"))
WLEOW <- plot_model(model.list[["LEOW"]], type = "pred", terms = c("Wind", "Habitat"))
WLISP <- plot_model(model.list[["LISP"]], type = "pred", terms = c("Wind", "Habitat"))
WNSWO <- plot_model(model.list[["NSWO"]], type = "pred", terms = c("Wind", "Habitat"))
WOSFL <- plot_model(model.list[["OSFL"]], type = "pred", terms = c("Wind", "Habitat"))
WOVEN <- plot_model(model.list[["OVEN"]], type = "pred", terms = c("Wind", "Habitat"))
WPISI <- plot_model(model.list[["PISI"]], type = "pred", terms = c("Wind", "Habitat"))
WRBGR <- plot_model(model.list[["RBGR"]], type = "pred", terms = c("Wind", "Habitat"))
WRBNU <- plot_model(model.list[["RBNU"]], type = "pred", terms = c("Wind", "Habitat"))
WTEWA <- plot_model(model.list[["TEWA"]], type = "pred", terms = c("Wind", "Habitat"))
WWAVI <- plot_model(model.list[["WAVI"]], type = "pred", terms = c("Wind", "Habitat"))
WWETO <- plot_model(model.list[["WETO"]], type = "pred", terms = c("Wind", "Habitat"))
WWTSP <- plot_model(model.list[["WTSP"]], type = "pred", terms = c("Wind", "Habitat"))
WYERA <- plot_model(model.list[["YERA"]], type = "pred", terms = c("Wind", "Habitat"))


plot_model(pooled, type = "pred", terms = c("Wind", "Habitat"))

######################################################
####Monte Carlo Simulation and pairwise comparison####
######################################################
#requires previous analysis to be run first (model selection needed)
hbt.pw <- hbt
hbt.pw$x <- -hbt.pw$Distance^2
values.pw <- data.frame()

for(i in 1:length(model2)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model2[i], ]
  m <- glm(Detected ~ x + x:Habitat -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x)
    m_star <- glm(y ~ x + x:Habitat -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model2[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model3)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model3[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Temp -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Temp = data$Temp)
    m_star <- glm(y ~ x + x:Habitat + x:Temp -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model3[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model4)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model4[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Humidity -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Humidity = data$Humidity)
    m_star <- glm(y ~ x + x:Habitat + x:Humidity -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model4[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model5)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model5[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Temp + x:Humidity -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Temp = data$Temp, Humidity = data$Humidity)
    m_star <- glm(y ~ x + x:Habitat + x:Temp + x:Humidity -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9] + cf[10]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9] + cf[10]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9] + cf[10]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9] + cf[10]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9] + cf[10]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9] + cf[10]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9] + cf[10]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model5[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model8)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model8[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Temp + x:Wind -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Temp = data$Temp, Wind = data$Wind)
    m_star <- glm(y ~ x + x:Habitat + x:Temp + x:Wind -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9] + cf[10]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9] + cf[10]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9] + cf[10]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9] + cf[10]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9] + cf[10]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9] + cf[10]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9] + cf[10]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model8[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model10)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model10[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Wind -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Wind = data$Wind)
    m_star <- glm(y ~ x + x:Habitat + x:Wind -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model10[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model11)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model11[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Temp + x:Wind + x:Humidity -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Temp = data$Temp, Wind = data$Wind, Humidity = data$Humidity)
    m_star <- glm(y ~ x + x:Habitat + x:Temp + x:Wind + x:Humidity -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9] + cf[10] + cf[11]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9] + cf[10] + cf[11]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9] + cf[10] + cf[11]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9] + cf[10] + cf[11]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9] + cf[10] + cf[11]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9] + cf[10] + cf[11]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9] + cf[10] + cf[11]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model11[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model12)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model12[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Humidity + x:Wind -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Humidity = data$Humidity, Wind = data$Wind)
    m_star <- glm(y ~ x + x:Habitat + x:Humidity + x:Wind -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9] + cf[10]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9] + cf[10]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9] + cf[10]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9] + cf[10]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9] + cf[10]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9] + cf[10]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9] + cf[10]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model12[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

for(i in 1:length(model13)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == model13[i], ]
  m <- glm(Detected ~ x + x:Habitat + x:Temp + x:Wind -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f),1, f)
    df_star <- data.frame(y = y_star, Habitat = data$Habitat, x = data$x, Temp = data$Temp, Wind = data$Wind)
    m_star <- glm(y ~ x + x:Habitat + x:Temp + x:Wind -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1] + cf[9] + cf[10]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2] + cf[9] + cf[10]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3] + cf[9] + cf[10]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4] + cf[9] + cf[10]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5] + cf[9] + cf[10]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6] + cf[9] + cf[10]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7] + cf[9] + cf[10]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  values <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
                  sp = as.character(paste(model13[i], sep = "")))
  values.pw <- rbind(values.pw, values)
}

values.pw.cleaned <- values.pw
values.pw.cleaned$x <- as.numeric(as.character(values.pw.cleaned$x))
values.pw.cleaned <- values.pw.cleaned[complete.cases(values.pw.cleaned), ]

mc.check <- count(values.pw.cleaned, sp, type)
mc.keep <- mc.check[mc.check$n > 950, ]
mc.keep$temp <- paste(mc.keep$sp, mc.keep$type)
values.pw.cleaned$temp <- paste(values.pw.cleaned$sp, values.pw.cleaned$type)

values.pw.cleaned <- values.pw.cleaned %>%
  filter(temp %in% mc.keep$temp)
values.pw.cleaned$sp <- droplevels(values.pw.cleaned$sp)

veg.list <- unique(values.pw.cleaned$type)
sp.list <- unique(values.pw.cleaned$sp)
mc.sim <- data.frame()

for(i in 1:length(sp.list)){
  sp.data <- values.pw.cleaned[values.pw.cleaned$sp == sp.list[i], ]
  for(j in 1:length(veg.list)){
    veg.sp.data <- sp.data[sp.data$type == veg.list[j], ]
    a <- mean(veg.sp.data$x) #mean EDR
    b <- quantile(veg.sp.data$x, 0.05) #5% CI
    c <- quantile(veg.sp.data$x, 0.95) #95% CI
    d <- as.character(paste(veg.list[j], sep = "")) #veg type
    e <- as.character(paste(sp.list[i], sep = "")) #sp type
    mc.sim.temp <- data.frame(mean = a, lower = b, upper = c, type = d, sp = e)
    mc.sim <- rbind(mc.sim, mc.sim.temp)
  }
}

mc.sim.clean <- mc.sim[complete.cases(mc.sim), ]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "CATO"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "BOOW"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "1000Hz"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "2828Hz"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "BADO"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "GGOW"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "NSWO"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "WETO"),]
mc.sim.clean <- mc.sim.clean[!(mc.sim.clean$sp == "RBGR" & mc.sim.clean$type == "grass"),]

plot.list <- list()
sp.list <- unique(mc.sim.clean$sp)

for(i in 1:length(sp.list)){
  plot.data <- mc.sim.clean[mc.sim.clean$sp == sp.list[i], ]
  mc.sim.plot <- ggplot(plot.data, aes(type, mean, factor = "type")) +
    geom_point(aes(x=type, y=mean), size = 0.4) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width=0.2, position=position_dodge(0.9), size = 0.3) +
    ylab("EDR Estimates") + xlab("Vegetation Type") +
    theme(legend.position="none") +
    theme(panel.background = element_rect(fill = 'white')) +
    theme(axis.line = element_line(colour = "black", size = 0.3)) +
    theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6), axis.ticks.x = element_blank()) +
    theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
    theme(legend.title=element_blank()) +
    theme(plot.title = element_text(hjust=0.03, size=8)) +
    theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
    theme(legend.key.size = unit(0.3, "cm")) +
    theme(legend.text = element_text(size=6)) +
    ggtitle(as.character(paste(sp.list[i], sep = "")))
  
  plot.list[[as.character(paste(sp.list[i], sep = ""))]] = mc.sim.plot
  #file_name = paste(i, ".tiff", sep = "")
  #tiff(width=7,height=4, units="in", res=600, file=file_name)
  #print(mc.sim.plot)
  #dev.off()
}

pairwise <- mc.sim.clean %>%
  group_by(sp) %>%
  summarise(no_rows = length(sp))

tiff(width = 7, height = 4, units = "in", res = 600, file = "Figure 3 - habitat differences.tiff")
grid.arrange(plot.list[["LISP"]], plot.list[["OVEN"]], plot.list[["BEKI"]], plot.list[["RBNU"]])
dev.off()

####summarize commonality####
pooled.c <- data.frame(p.comF)
pooled.c <- cbind(predictors = rownames(pooled.c), pooled.c)
rownames(pooled.c) <- 1:nrow(pooled.c)

model.commonality[["pooled"]] = pooled.c

summarize.commonality <- data.frame()
for(i in 1:length(model.commonality)){
  data <- data.frame(model.commonality[i])
  data <- cbind(predictors = rownames(data), data)
  rownames(data) <- 1:nrow(data)
  data <- data
  dist <- as.numeric(data[data$predictors == "Distance", ][3])
  hab <- as.numeric(data[data$predictors == "Habitat", ][3])
  temp <- as.numeric(data[data$predictors == "Temp", ][3])
  wind <- as.numeric(data[data$predictors == "Wind", ][3])
  hum <- as.numeric(data[data$predictors == "Humidity", ][3])
  sp <- as.character(paste(names(model.commonality[i]), sep = ""))
  temp <- data.frame(Distance = dist, Habitat = hab, Temp = temp, Wind = wind, 
                     Humidity = hum, Species = sp)
  summarize.commonality <- rbind(summarize.commonality, temp)
}

summarize.commonality <- merge(summarize.commonality, meta, 
                               by.x = "Species", by.y = "Sound")

dist.c.plot <- ggplot(summarize.commonality, aes(MinFreq, Distance)) +
  geom_point(aes(x=MinFreq, y=Distance), size = 0.4) +
  stat_smooth(method = glm, colour = "black", size = 0.3) +
  ylab("% total unique variation explained") + xlab("Minimum Frequency (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("A")
dist.c.plot

veg.c.plot <- ggplot(summarize.commonality, aes(MinFreq, Habitat)) +
  geom_point(aes(x=MinFreq, y=Habitat), size = 0.4) +
  stat_smooth(method = glm, colour = "black", size = 0.3) +
  ylab("% total unique variation explained") + xlab("Minimum Frequency (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("B")
veg.c.plot

tiff(width = 3.5, height = 4, units = "in", res = 600, file = "Figure 2 - commonality - frequency effects.tiff")
grid.arrange(dist.c.plot, veg.c.plot, ncol = 1)
dev.off()

summary(glm(Distance ~ MinFreq, data = summarize.commonality))
summary(glm(Habitat ~ MinFreq, data = summarize.commonality))
summary(glm(Temp ~ MinFreq, data = summarize.commonality))
summary(glm(Wind ~ MinFreq, data = summarize.commonality))
summary(glm(Humidity ~ MinFreq, data = summarize.commonality))

summary(glm(Distance ~ Bandwidth, data = summarize.commonality))
summary(glm(Habitat ~ Bandwidth, data = summarize.commonality))
summary(glm(Temp ~ Bandwidth, data = summarize.commonality))
summary(glm(Wind ~ Bandwidth, data = summarize.commonality))
summary(glm(Humidity ~ Bandwidth, data = summarize.commonality))

wind.c.plot <- ggplot(summarize.commonality, aes(Bandwidth, Wind)) +
  geom_point(aes(x=Bandwidth, y=Wind)) +
  stat_smooth(method = glm) +
  ylab("% total unique variation explained") + xlab("Bandwidth (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.title.x=element_text(size=8), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("Wind")
wind.c.plot

dist.c.plot2 <- ggplot(summarize.commonality, aes(Bandwidth, Distance)) +
  geom_point(aes(x=Bandwidth, y=Distance)) +
  stat_smooth(method = glm) +
  ylab("% total unique variation explained") + xlab("Bandwidth (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.title.x=element_text(size=8), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("Distance")
dist.c.plot2


####testing area surveyed biases####
#species pooled
edr.pooled.spp <- list()
data <- hbt.pw
data <- data[grep("Hz", data$Species, invert = TRUE), ]
m <- glm(Detected ~ x -1, data = data, family = binomial(link = cloglog))
f <- fitted(m)
for(i in 1:1000){
  y_star <- rbinom(length(f), 1, f)
  df_star <- data.frame(y = y_star, x = data$x)
  m_star <- glm(y ~ x -1, data = df_star, family = binomial(link = cloglog))
  cf <- coef(m_star)
  edr.pooled.spp[[i]] <- sqrt(1/(cf[1]))
}
edr.pooled.spp_mat <- do.call(rbind, edr.pooled.spp)

#species split
edr.list <- list()
edr.split.spp <- data.frame()
sp.list <- unique(hbt.pw$Species)
for(i in 1:length(sp.list)){
  data <- hbt.pw[hbt.pw$Species == sp.list[i], ]
  m <- glm(Detected ~ x -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f), 1, f)
    df_star <- data.frame(y = y_star, x = data$x)
    m_star <- glm(y ~ x -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    edr.list[[j]] <- sqrt(1/(cf[1]))
  }
  edr.split.spp_mat <- cbind(do.call(rbind, edr.list), type = as.character(paste(sp.list[i], sep = "")))
  edr.split.spp <- rbind(edr.split.spp, edr.split.spp_mat)
}

#species split, habitat open/closed
hbt.oc <- hbt.pw
hbt.oc$oc <- ifelse(hbt.oc$Habitat == "B" | hbt.oc$Habitat == "C" | hbt.oc$Habitat == "D", "closed", "open")
edr.oc.spp <- data.frame()
for(i in 1:length(sp.list)){
  edr.list.o <- list()
  edr.list.c <- list()
  data <- hbt.oc[hbt.oc$Species == sp.list[i],]
  m <- glm(Detected ~ x + x:oc -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f), 1, f)
    df_star <- data.frame(y = y_star, x = data$x, oc = data$oc)
    m_star <- glm(y ~ x + x:oc -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    edr.list.o[[j]] <- sqrt(1/(cf[1]))
    edr.list.c[[j]] <- sqrt(1/(cf[1] + cf[2]))
  }
  edr.list.o_mat <- cbind(do.call(rbind, edr.list.o), type = "open")
  edr.list.c_mat <- cbind(do.call(rbind, edr.list.c), type = "closed")
  temp <- cbind(rbind(edr.list.o_mat, edr.list.c_mat), sp = as.character(paste(sp.list[i], sep = "")))
  edr.oc.spp <- rbind(edr.oc.spp, temp)
}

#species split, habitat split
edr.hbt.spp <- data.frame()
for(i in 1:length(sp.list)){
  ROAD <- list()
  BOG <- list()
  CONIFER <- list()
  DECIDUOUS <- list()
  FEN <- list()
  GRASS <- list()
  EDGE <- list()
  data <- hbt.pw[hbt.pw$Species == sp.list[i],]
  m <- glm(Detected ~ x + x:Habitat -1, data = data, family = binomial(link = cloglog))
  f <- fitted(m)
  for(j in 1:1000){
    y_star <- rbinom(length(f), 1, f)
    df_star <- data.frame(y = y_star, x = data$x, hbt = data$Habitat)
    m_star <- glm(y ~ x + x:hbt -1, data = df_star, family = binomial(link = cloglog))
    cf <- coef(m_star)
    ROAD[[j]] <- sqrt(1/(cf[1]))
    BOG[[j]] <- sqrt(1/(cf[1] + cf[2]))
    CONIFER[[j]] <- sqrt(1/(cf[1] + cf[3]))
    DECIDUOUS[[j]] <- sqrt(1/(cf[1] + cf[4]))
    FEN[[j]] <- sqrt(1/(cf[1] + cf[5]))
    GRASS[[j]] <- sqrt(1/(cf[1] + cf[6]))
    EDGE[[j]] <- sqrt(1/(cf[1] + cf[7]))
  }
  ROAD_mat <- cbind(do.call(rbind, ROAD), type = "road")
  BOG_mat <- cbind(do.call(rbind, BOG), type = "bog")
  CONIFER_mat <- cbind(do.call(rbind, CONIFER), type = "conifer")
  DECIDUOUS_mat <- cbind(do.call(rbind, DECIDUOUS), type = "deciduous")
  FEN_mat <- cbind(do.call(rbind, FEN), type = "fen")
  GRASS_mat <- cbind(do.call(rbind, GRASS), type = "grass")
  EDGE_mat <- cbind(do.call(rbind, EDGE), type = "edge")
  temp <- cbind(rbind(ROAD_mat, BOG_mat, CONIFER_mat, DECIDUOUS_mat, FEN_mat, GRASS_mat, EDGE_mat), 
        sp = as.character(paste(sp.list[i], sep = "")))
  edr.hbt.spp <- rbind(edr.hbt.spp, temp)
}

#edr.pooled.spp_mat = all species pooled
#edr.split.spp = separate species, all habitat pooled
#edr.oc.spp = separate species, habitat split open/closed
#edr.hbt.spp = separate species, habitat split by class

allspp.allveg <- data.frame(edr.pooled.spp_mat)
sepspp.allveg <- edr.split.spp
sepspp.ocveg <- edr.oc.spp
sepspp.sepveg <- edr.hbt.spp

sepspp.allveg <- sepspp.allveg[grep("Hz", sepspp.allveg$type, invert = TRUE), ]
sepspp.ocveg <- sepspp.ocveg[grep("Hz", sepspp.ocveg$sp, invert = TRUE), ]
sepspp.sepveg <- sepspp.sepveg[grep("Hz", sepspp.sepveg$sp, invert = TRUE), ]

sepspp.allveg$x <- as.numeric(as.character(sepspp.allveg$x))
sepspp.allveg <- sepspp.allveg[complete.cases(sepspp.allveg), ]
sepspp.allveg.check <- count(sepspp.allveg, type)

sepspp.ocveg$x <- as.numeric(as.character(sepspp.ocveg$x))
sepspp.ocveg <- sepspp.ocveg[complete.cases(sepspp.ocveg), ]
sepspp.ocveg.check <- count(sepspp.ocveg, sp, type)

sepspp.sepveg$x <- as.numeric(as.character(sepspp.sepveg$x))
sepspp.sepveg <- sepspp.sepveg[complete.cases(sepspp.sepveg), ]
sepspp.sepveg.check <- count(sepspp.sepveg, sp, type)
sepspp.sepveg.keep <- sepspp.sepveg.check[sepspp.sepveg.check$n > 950, ]
sepspp.sepveg.keep$temp <- paste(sepspp.sepveg.keep$sp, sepspp.sepveg.keep$type)
sepspp.sepveg$temp <- paste(sepspp.sepveg$sp, sepspp.sepveg$type)
sepspp.sepveg <- sepspp.sepveg %>%
  filter( temp %in% sepspp.sepveg.keep$temp)
sepspp.sepveg$sp <- droplevels(sepspp.sepveg$sp)

allspp.allveg$a <- pi*(allspp.allveg$x)^2
sepspp.allveg$a <- pi*(sepspp.allveg$x)^2
sepspp.ocveg$a <- pi*(sepspp.ocveg$x)^2
sepspp.sepveg$a <- pi*(sepspp.sepveg$x)^2

sp.list <- unique(sepspp.allveg$type)
sepspp.allveg.mc <- data.frame()

for(i in 1:length(sp.list)){
  sp.data <- sepspp.allveg[sepspp.allveg$type == sp.list[i],]
  a <- mean(sp.data$a)
  b <- sd(sp.data$a)
  c <- as.character(paste(sp.list[i], sep = ""))
  temp <- data.frame(mean = a, sd = b, sp = c)
  sepspp.allveg.mc <- rbind(sepspp.allveg.mc, temp)
}

sp.list <- unique(sepspp.ocveg$sp)
veg.list <- unique(sepspp.ocveg$type)
sepspp.ocveg.mc <- data.frame()

for(i in 1:length(sp.list)){
  sp.data <- sepspp.ocveg[sepspp.ocveg$sp == sp.list[i],]
  for(j in 1:length(veg.list)){
    veg.data <- sp.data[sp.data$type == veg.list[j],]
    a <- mean(veg.data$a)
    b <- sd(veg.data$a)
    c <- as.character(paste(sp.list[i], sep = ""))
    d <- as.character(paste(veg.list[j], sep = ""))
    temp <- data.frame(mean = a, sd = b, sp = c, veg = d)
    sepspp.ocveg.mc <- rbind(sepspp.ocveg.mc, temp)
  }
}

sp.list <- unique(sepspp.sepveg$sp)
veg.list <- unique(sepspp.sepveg$type)
sepspp.sepveg.mc <- data.frame()

for(i in 1:length(sp.list)){
  sp.data <- sepspp.sepveg[sepspp.sepveg$sp == sp.list[i],]
  for(j in 1:length(veg.list)){
    veg.data <- sp.data[sp.data$type == veg.list[j],]
    a <- mean(veg.data$a)
    b <- sd(veg.data$a)
    c <- as.character(paste(sp.list[i], sep = ""))
    d <- as.character(paste(veg.list[j], sep = ""))
    temp <- data.frame(mean = a, sd = b, sp = c, veg = d)
    sepspp.sepveg.mc <- rbind(sepspp.sepveg.mc, temp)
  }
}

allspp.allveg.mc <- data.frame(mean = mean(allspp.allveg$a), sd = sd(allspp.allveg$a))

#all pooled mean = 911826m^2, sd = 11464.58
sepspp.allveg.mc <- read.csv(file = "sepspp_allveg.csv")
sepspp.ocveg.mc <- read.csv(file = "sepspp_ocveg.csv")
sepspp.sepveg.mc <- read.csv(file = "sepspp_sepveg.csv")

# spp.bias <- sepspp.allveg.mc
# spp.bias$bias <- spp.bias$mean/911826
# 
# oc.bias1 <- sepspp.ocveg.mc
# sp.list <- unique(oc.bias1$sp)
# oc.bias <- data.frame()
# for(i in 1:length(sp.list)){
#   data <- oc.bias1[oc.bias1$sp == sp.list[i],]
#   data1 <- spp.bias[spp.bias$sp == sp.list[i],]
#   x <- data1$mean
#   data$bias <- data$mean/x
#   oc.bias <- rbind(oc.bias, data)
# }
# 
# veg.bias1 <- sepspp.sepveg.mc
# veg.bias1$x1[(veg.bias1$veg == "road") | (veg.bias1$veg == "edge") | (veg.bias1$veg == "fen") | (veg.bias1$veg == "grass")] <- "open"
# veg.bias1$x1[(veg.bias1$veg == "deciduous") | (veg.bias1$veg == "conifer") | (veg.bias1$veg == "bog")] <- "closed"
# veg.bias1$x2 <- paste(veg.bias1$sp, veg.bias1$x1, sep = ".")
# oc.bias$x2 <- paste(oc.bias$sp, oc.bias$veg, sep = ".")
# sp.list <- unique(veg.bias1$x2)
# veg.bias <- data.frame()
# for(i in 1:length(sp.list)){
#   data <- veg.bias1[veg.bias1$x2 == sp.list[i],]
#   data1 <- oc.bias[oc.bias$x2 == sp.list[i],]
#   x <- data1$mean
#   data$bias <- data$mean/x
#   veg.bias <- rbind(veg.bias, data)
# }
# veg.bias1 <- veg.bias
# sp.list <- unique(veg.bias1$sp)
# veg.bias <- data.frame()
# for(i in 1:length(sp.list)){
#   data <- veg.bias1[veg.bias1$sp == sp.list[i],]
#   data1 <- spp.bias[spp.bias$sp == sp.list[i],]
#   x <- data1$mean
#   data$bias1 <- data$mean/x
#   veg.bias <- rbind(veg.bias, data)
# }
# veg.bias$bias2 <- veg.bias$mean/911826
# 
# spp.bias = as.data.table(spp.bias)
# spp.bias[bias < 1, magnitude := (-(1/bias))+1]
# spp.bias[bias >= 1, magnitude := bias-1]
# spp.bias$x <- "spp"
# 
# oc.bias = as.data.table(oc.bias)
# oc.bias[bias < 1, magnitude := (-(1/bias))+1]
# oc.bias[bias >= 1, magnitude := bias-1]
# oc.bias$x <- "spp"
# 
# veg.bias = as.data.table(veg.bias)
# veg.bias[bias < 1, magnitude := (-(1/bias))+1]
# veg.bias[bias >= 1, magnitude := bias-1]
# veg.bias[bias1 < 1, magnitude1 := (-(1/bias1))+1]
# veg.bias[bias1 >= 1, magnitude1 := bias1-1]
# veg.bias$x <- "spp"
# 
# oc.bias$abs.magnitude <- abs(oc.bias$magnitude)
# veg.bias$abs.magnitude <- abs(veg.bias$magnitude)
# veg.bias$abs.magnitude1 <- abs(veg.bias$magnitude1)
# 
# spp.plot <- ggplot(spp.bias, aes(x, magnitude)) +
#   geom_jitter(width = 0.01, size = 0.3) +
#   geom_text_repel(size=1.3, box.padding=unit(0.1, "lines"), segment.size=0.2, aes(label = spp.bias$sp)) +
#   #geom_point(size=0.3) +
#   #ggtitle("OVEN Manual Measure") +
#   ylab("Estimate Bias") + xlab("By Species") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme(legend.position="none")+
#   theme(panel.background = element_rect(fill = 'white')) +
#   theme(axis.line = element_line(colour = "black", size = 0.3)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   #theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
#   theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
#   theme(legend.title=element_blank()) +
#   theme(plot.title = element_text(hjust=0.03, size=8)) +
#   theme(axis.ticks = element_line(colour = "black", size= 0.3))+
#   theme(legend.key.size = unit(0.3, "cm"))+
#   theme(legend.text = element_text(size=6)) +
#   coord_cartesian(ylim = c(-8,5))
# spp.plot
# oc.plot <- ggplot(oc.bias, aes(x, magnitude, col = veg)) +
#   geom_jitter(width = 0.01, size = 0.3) +
#   geom_text_repel(size=1.3, box.padding=unit(0.1, "lines"), segment.size=0.2, aes(label = oc.bias$sp)) +
#   #geom_point(size=0.3) +
#   #ggtitle("OVEN Manual Measure") +
#   ylab("") + xlab("By canopy closure") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme(legend.position="right")+
#   theme(panel.background = element_rect(fill = 'white')) +
#   theme(axis.line = element_line(colour = "black", size = 0.3)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
#   #theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
#   #theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
#   theme(legend.title=element_blank()) +
#   theme(plot.title = element_text(hjust=0.03, size=8)) +
#   theme(axis.ticks = element_line(colour = "black", size= 0.3))+
#   theme(legend.key.size = unit(0.3, "cm"))+
#   theme(legend.text = element_text(size=6)) +
#   coord_cartesian(ylim = c(-8,5))
# oc.plot
# veg.plot <- ggplot(veg.bias, aes(veg, magnitude1, col = veg)) +
#   geom_jitter(width = 0.01, size = 0.3) +
#   geom_text_repel(size=1.3, box.padding=unit(0.1, "lines"), segment.size=0.2, aes(label = veg.bias$sp)) +
#   #geom_point(size=0.3) +
#   #ggtitle("OVEN Manual Measure") +
#   ylab("") + xlab("By habitat class") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme(legend.position="right")+
#   theme(panel.background = element_rect(fill = 'white')) +
#   theme(axis.line = element_line(colour = "black", size = 0.3)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
#   #theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
#   #theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
#   theme(legend.title=element_blank()) +
#   theme(plot.title = element_text(hjust=0.03, size=8)) +
#   theme(axis.ticks = element_line(colour = "black", size= 0.3))+
#   theme(legend.key.size = unit(0.3, "cm"))+
#   theme(legend.text = element_text(size=6)) +
#   coord_cartesian(ylim = c(-8,5))
# veg.plot
# 
# grid.arrange(spp.plot, oc.plot, veg.plot, ncol=3)
# 
# tiff(width=7,height=4, units="in", res=600, file="bias.tiff")
# grid.arrange(spp.plot, oc.plot, veg.plot, ncol=3)
# dev.off() 
# 
# #mean edr for each category
# oc.o <- oc.bias[oc.bias$veg == "open", ]
# oc.c <- oc.bias[oc.bias$veg == "closed", ]
# vop <- cbind(mean(oc.o$magnitude), sd(oc.o$magnitude), "open")
# vcl <- cbind(mean(oc.c$magnitude), sd(oc.c$magnitude), "closed")
# 
# veg.r <- veg.bias[veg.bias$veg == "road",]
# veg.e <- veg.bias[veg.bias$veg == "edge",]
# veg.g <- veg.bias[veg.bias$veg == "grass",]
# veg.f <- veg.bias[veg.bias$veg == "fen",]
# veg.c <- veg.bias[veg.bias$veg == "conifer",]
# veg.d <- veg.bias[veg.bias$veg == "deciduous",]
# veg.b <- veg.bias[veg.bias$veg == "bog",]
# 
# vr <- cbind(mean(veg.r$magnitude1), sd(veg.r$magnitude1), "road")
# ve <- cbind(mean(veg.e$magnitude1), sd(veg.e$magnitude1), "edge")
# vg <- cbind(mean(veg.g$magnitude1, na.rm = TRUE), sd(veg.g$magnitude1, na.rm = TRUE), "grass")
# vf <- cbind(mean(veg.f$magnitude1), sd(veg.f$magnitude1), "fen")
# vc <- cbind(mean(veg.c$magnitude1), sd(veg.c$magnitude1), "conifer")
# vd <- cbind(mean(veg.d$magnitude1), sd(veg.d$magnitude1), "deciduous")
# vb <- cbind(mean(veg.b$magnitude1), sd(veg.b$magnitude1), "bog")
# 
# correct <- data.frame()
# correct <- rbind(vop,vcl,vr,ve,vg,vf,vc,vd,vb)
# 
# 
# meta <- read.csv(file = "species_parameters.csv", header = TRUE)
# veg.bias <- merge(veg.bias, meta, by.x = "sp", by.y = "Sound")
# veg.bias$Bandwidth.st <- scale(veg.bias$Bandwidth, center=TRUE, scale = TRUE)
# veg.bias$MinFreq.st <- scale(veg.bias$MinFreq, center=TRUE, scale = TRUE)
# 
# m1 <- glm(bias1 ~ veg + MinFreq + Bandwidth, data = veg.bias)
# m2 <- glm(log(bias1) ~ sp + veg, data = veg.bias)
# m3 <- glm(bias1 ~ veg + MinFreq + Bandwidth, data = veg.bias, family=gaussian(link=log))
# ms <- model.sel(m1,m3)
# 
# null <- glm(bias1 ~ veg, data = veg.bias, family=gaussian(link=log))
# m1 <- glm(bias1 ~ veg + MinFreq, data = veg.bias, family=gaussian(link=log))
# m2 <- glm(bias1 ~ veg + Bandwidth, data = veg.bias, family=gaussian(link=log))
# m3 <- glm(bias1 ~ veg + MinFreq + Bandwidth, data = veg.bias, family = gaussian(link=log))
# ms <- model.sel(null,m1,m2,m3)
# 
# #LOOCV
# sp.list <- unique(veg.bias$sp)
# v.predict <- data.frame()
# v.r2 <- data.frame()
# for(i in 1:length(sp.list)){
#   train <- veg.bias[veg.bias$sp != sp.list[i], ]
#   test <- veg.bias[veg.bias$sp == sp.list[i], ]
#   m <- glm(bias1 ~ veg + MinFreq + Bandwidth, data = train, family = gaussian(link = log))
#   p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
#   v.predict <- rbind(v.predict, p.test)
#   validate_model <- lm(bias1 ~ fitted, data = p.test)
#   adj.r2 <- cbind(summary(validate_model)$adj.r.squared, as.character(paste(sp.list[i])))
#   v.r2 <- rbind(v.r2, adj.r2)
# }
# 
# v.predict[fitted < 1, m.predict := (-(1/fitted))+1]
# v.predict[fitted >= 1, m.predict := fitted-1]
# 
# v.predict$error <- abs(v.predict$magnitude1 - v.predict$m.predict)
# 
# p.plot <- ggplot(v.predict, aes(veg, m.predict, col = veg)) +
#   geom_jitter(width = 0.01, size = 0.3) +
#   geom_text_repel(size=1.3, box.padding=unit(0.1, "lines"), segment.size=0.2, aes(label = veg.bias$sp)) +
#   #geom_point(size=0.3) +
#   #ggtitle("OVEN Manual Measure") +
#   ylab("") + xlab("By habitat class") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme(legend.position="right")+
#   theme(panel.background = element_rect(fill = 'white')) +
#   theme(axis.line = element_line(colour = "black", size = 0.3)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
#   #theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
#   #theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
#   theme(legend.title=element_blank()) +
#   theme(plot.title = element_text(hjust=0.03, size=8)) +
#   theme(axis.ticks = element_line(colour = "black", size= 0.3))+
#   theme(legend.key.size = unit(0.3, "cm"))+
#   theme(legend.text = element_text(size=6)) +
#   coord_cartesian(ylim = c(-8,5))
# p.plot
# 
# grid.arrange(veg.plot, p.plot, ncol = 2)
# 
# veg.bias$edr <- sqrt(veg.bias$mean/pi)
# m1 <- glm(edr ~ veg + MinFreq + Bandwidth, data = veg.bias)
# m3 <- glm(edr ~ veg + MinFreq + Bandwidth, data = veg.bias, family=gaussian(link=log))
# ms <- model.sel(m1,m3)
# 
# null <- glm(edr ~ veg, data = veg.bias, family=gaussian(link=log))
# m1 <- glm(edr ~ veg + MinFreq, data = veg.bias, family=gaussian(link=log))
# m2 <- glm(edr ~ veg + Bandwidth, data = veg.bias, family=gaussian(link=log))
# m3 <- glm(edr ~ veg + MinFreq + Bandwidth, data = veg.bias, family = gaussian(link=log))
# m4 <- glm(edr ~ veg*MinFreq + Bandwidth, data = veg.bias, family = gaussian(link=log))
# m5 <- glm(edr ~ veg*Bandwidth + MinFreq, data = veg.bias, family = gaussian(link=log))
# m6 <- glm(edr ~ veg*MinFreq + veg*Bandwidth, data = veg.bias, family = gaussian(link=log))
# ms <- model.sel(null,m1,m2,m3,m4,m5,m6)
# 
# #LOOCV
# sp.list <- unique(veg.bias$sp)
# e.predict <- data.frame()
# e.r2 <- data.frame()
# for(i in 1:length(sp.list)){
#   train <- veg.bias[veg.bias$sp != sp.list[i], ]
#   test <- veg.bias[veg.bias$sp == sp.list[i], ]
#   m <- glm(edr ~ veg + MinFreq + Bandwidth, data = train, family = gaussian(link = log))
#   p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
#   e.predict <- rbind(e.predict, p.test)
#   validate_model <- lm(edr ~ fitted, data = p.test)
#   adj.r2 <- summary(validate_model)$adj.r.squared
#   e.r2 <- rbind(e.r2, adj.r2)
# }
# 
# e.predict$error <- abs(e.predict$edr - e.predict$fitted)
# 
# oc.bias <- merge(oc.bias, meta, by.x = "sp", by.y = "Sound")
# oc.bias$Bandwidth.st <- scale(oc.bias$Bandwidth, center=TRUE, scale = TRUE)
# oc.bias$MinFreq.st <- scale(oc.bias$MinFreq, center=TRUE, scale = TRUE)
# 
# m1 <- glm(bias ~ veg + MinFreq + Bandwidth, data = oc.bias)
# m2 <- glm(log(bias) ~ sp + veg, data = oc.bias)
# m3 <- glm(bias ~ veg + MinFreq + Bandwidth, data = oc.bias, family=gaussian(link=log))
# ms <- model.sel(m1,m3)
# 
# null <- glm(bias ~ veg, data = oc.bias, family=gaussian(link=log))
# m1 <- glm(bias ~ veg + MinFreq, data = oc.bias, family=gaussian(link=log))
# m2 <- glm(bias ~ veg + Bandwidth, data = oc.bias, family=gaussian(link=log))
# m3 <- glm(bias ~ veg + MinFreq + Bandwidth, data = oc.bias, family = gaussian(link=log))
# ms <- model.sel(null,m1,m2,m3)
# 
# #LOOCV
# sp.list <- unique(oc.bias$sp)
# oc.predict <- data.frame()
# oc.r2 <- data.frame()
# for(i in 1:length(sp.list)){
#   train <- oc.bias[oc.bias$sp != sp.list[i], ]
#   test <- oc.bias[oc.bias$sp == sp.list[i], ]
#   m <- glm(bias ~ veg + MinFreq + Bandwidth, data = train, family = gaussian(link = log))
#   p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
#   oc.predict <- rbind(oc.predict, p.test)
#   validate_model <- lm(bias ~ fitted, data = p.test)
#   adj.r2 <- summary(validate_model)$adj.r.squared
#   oc.r2 <- rbind(oc.r2, adj.r2)
# }
# 
# oc.predict[fitted < 1, m.predict := (-(1/fitted))+1]
# oc.predict[fitted >= 1, m.predict := fitted-1]
# 
# oc.predict$error <- abs(oc.predict$magnitude - oc.predict$m.predict)
# 
# poc.plot <- ggplot(oc.predict, aes(x, m.predict, col = veg)) +
#   geom_jitter(width = 0.01, size = 0.3) +
#   geom_text_repel(size=1.3, box.padding=unit(0.1, "lines"), segment.size=0.2, aes(label = oc.predict$sp)) +
#   #geom_point(size=0.3) +
#   #ggtitle("OVEN Manual Measure") +
#   ylab("") + xlab("By habitat class") +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   theme(legend.position="right")+
#   theme(panel.background = element_rect(fill = 'white')) +
#   theme(axis.line = element_line(colour = "black", size = 0.3)) +
#   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
#   theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()) +
#   #theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6)) +
#   #theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
#   theme(legend.title=element_blank()) +
#   theme(plot.title = element_text(hjust=0.03, size=8)) +
#   theme(axis.ticks = element_line(colour = "black", size= 0.3))+
#   theme(legend.key.size = unit(0.3, "cm"))+
#   theme(legend.text = element_text(size=6)) +
#   coord_cartesian(ylim = c(-8,5))
# poc.plot
# 
# grid.arrange(oc.plot, poc.plot, ncol = 2)
# 
# ####LOOCV EDR####
# sp.list <- unique(sepspp.sepveg$sp)
# veg.list <- unique(sepspp.sepveg$type)
# sepspp.sepveg.edr <- data.frame()
# 
# for(i in 1:length(sp.list)){
#   sp.data <- sepspp.sepveg[sepspp.sepveg$sp == sp.list[i],]
#   for(j in 1:length(veg.list)){
#     veg.data <- sp.data[sp.data$type == veg.list[j],]
#     a <- mean(veg.data$x)
#     b <- sd(veg.data$x)
#     c <- as.character(paste(sp.list[i], sep = ""))
#     d <- as.character(paste(veg.list[j], sep = ""))
#     temp <- data.frame(mean = a, sd = b, sp = c, veg = d)
#     sepspp.sepveg.edr <- rbind(sepspp.sepveg.edr, temp)
#   }
# }
# 
# meta <- read.csv(file = "species_parameters.csv", header = TRUE)
# sepspp.sepveg.edr <- merge(sepspp.sepveg.edr, meta, by.x = "sp", by.y = "Sound")
# sepspp.sepveg.edr$Bandwidth.st <- scale(sepspp.sepveg.edr$Bandwidth, center=TRUE, scale = TRUE)
# sepspp.sepveg.edr$MinFreq.st <- scale(sepspp.sepveg.edr$MinFreq, center=TRUE, scale = TRUE)
# 
# m1 <- glm(mean ~ veg + MinFreq + Bandwidth, data = sepspp.sepveg.edr)
# m2 <- glm(log(mean) ~ sp + veg, data = sepspp.sepveg.edr)
# m3 <- glm(mean ~ veg + MinFreq + Bandwidth, data = sepspp.sepveg.edr, family=gaussian(link=log))
# ms <- model.sel(m1,m3)
# 
# null <- glm(mean ~ veg, data = sepspp.sepveg.edr, family=gaussian(link=log))
# m1 <- glm(mean ~ veg + MinFreq, data = sepspp.sepveg.edr, family=gaussian(link=log))
# m2 <- glm(mean ~ veg + Bandwidth, data = sepspp.sepveg.edr, family=gaussian(link=log))
# m3 <- glm(mean ~ veg + MinFreq + Bandwidth, data = sepspp.sepveg.edr, family = gaussian(link=log))
# ms <- model.sel(null,m1,m2,m3)
# 
# sp.list <- unique(sepspp.sepveg.edr$sp)
# v.predict <- data.frame()
# v.r2 <- data.frame()
# for(i in 1:length(sp.list)){
#   train <- sepspp.sepveg.edr[sepspp.sepveg.edr$sp != sp.list[i], ]
#   test <- sepspp.sepveg.edr[sepspp.sepveg.edr$sp == sp.list[i], ]
#   m <- glm(mean ~ veg + MinFreq + Bandwidth, data = train, family = gaussian(link = log))
#   p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
#   v.predict <- rbind(v.predict, p.test)
#   validate_model <- lm(mean ~ fitted, data = p.test)
#   adj.r2 <- cbind(summary(validate_model)$adj.r.squared, as.character(paste(sp.list[i])))
#   v.r2 <- rbind(v.r2, adj.r2)
# }
# v.r2$V1 <- as.numeric(as.character(v.r2$V1))
# v.predict$error <- abs(v.predict$mean - v.predict$fitted)
# v.predict$p.error <- (v.predict$error/v.predict$fitted)

####try2

sp.list <- unique(sepspp.sepveg.mc$sp)
v.corr <- data.frame()
for(i in 1:length(sp.list)){
  data <- sepspp.sepveg.mc[sepspp.sepveg.mc$sp == sp.list[i], ]
  ref <- data[data$veg == "edge", ]
  data$corr <- ref$mean/data$mean #reference level is always numerator, levels to be corrected are in denominator
  v.corr <- rbind(v.corr, data)
}

sp.list <- unique(sepspp.ocveg.mc$sp)
c.corr <- data.frame()
for(i in 1:length(sp.list)){
  data <- sepspp.ocveg.mc[sepspp.ocveg.mc$sp == sp.list[i], ]
  ref <- data[data$veg == "open", ]
  data$corr <- ref$mean/data$mean #reference level is always numerator, levels to be corrected are in denominator
  c.corr <- rbind(c.corr, data)
}

v.corr <- v.corr[v.corr$veg != "edge", ]
c.corr <- c.corr[c.corr$veg != "open", ]

#model sel
meta <- read.csv(file = "species_parameters.csv", header = TRUE)
v.corr <- merge(v.corr, meta, by.x = "sp", by.y = "Sound")
v.corr$Bandwidth.st <- scale(v.corr$Bandwidth, center=TRUE, scale = TRUE)
v.corr$MinFreq.st <- scale(v.corr$MinFreq, center=TRUE, scale = TRUE)

m1 <- glm(corr ~ veg + MinFreq + Bandwidth, data = v.corr)
m2 <- glm(corr ~ veg + MinFreq + Bandwidth, data = v.corr, family=gaussian(link=log))
ms <- model.sel(m1,m2)

null <- glm(corr ~ veg, data = v.corr, family=gaussian(link=log))
m1 <- glm(corr ~ veg + MinFreq, data = v.corr, family=gaussian(link=log))
m2 <- glm(corr ~ veg + Bandwidth, data = v.corr, family=gaussian(link=log))
m3 <- glm(corr ~ veg + MinFreq + Bandwidth, data = v.corr, family = gaussian(link=log))
ms <- model.sel(null,m1,m2,m3)

sp.list <- unique(v.corr$sp)
v.predict <- data.frame()
v.r2 <- data.frame()
for(i in 1:length(sp.list)){
  train <- v.corr[v.corr$sp != sp.list[i], ]
  test <- v.corr[v.corr$sp == sp.list[i], ]
  m <- glm(corr ~ veg + MinFreq + Bandwidth, data = train, family = gaussian(link = log))
  p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
  v.predict <- rbind(v.predict, p.test)
  validate_model <- lm(mean ~ fitted, data = p.test)
  adj.r2 <- cbind(summary(validate_model)$adj.r.squared, as.character(paste(sp.list[i])))
  v.r2 <- rbind(v.r2, adj.r2)
}

v.r2$V1 <- as.numeric(as.character(v.r2$V1))
v.predict$error <- abs(v.predict$corr - v.predict$fitted)
v.predict$p.error <- (v.predict$error/v.predict$fitted)

c.corr <- merge(c.corr, meta, by.x = "sp", by.y = "Sound")
c.corr$Bandwidth.st <- scale(c.corr$Bandwidth, center=TRUE, scale = TRUE)
c.corr$MinFreq.st <- scale(c.corr$MinFreq, center=TRUE, scale = TRUE)

m1 <- glm(corr ~ MinFreq + Bandwidth, data = c.corr)
m2 <- glm(corr ~ MinFreq + Bandwidth, data = c.corr, family=gaussian(link=log))
ms <- model.sel(m1,m2)

m1 <- glm(corr ~ MinFreq, data = c.corr)
m2 <- glm(corr ~ Bandwidth, data = c.corr)
m3 <- glm(corr ~ MinFreq + Bandwidth, data = c.corr)
ms <- model.sel(m1,m2,m3)

sp.list <- unique(c.corr$sp)
c.predict <- data.frame()
c.r2 <- data.frame()
for(i in 1:length(sp.list)){
  train <- c.corr[c.corr$sp != sp.list[i], ]
  test <- c.corr[c.corr$sp == sp.list[i], ]
  m <- glm(corr ~ MinFreq, data = train)
  p.test <- cbind(test, fitted = predict(m, newdata = test, type = "response", allow.new.levels = TRUE, se.fit = FALSE))
  c.predict <- rbind(c.predict, p.test)
  validate_model <- lm(mean ~ fitted, data = p.test)
  adj.r2 <- cbind(summary(validate_model)$adj.r.squared, as.character(paste(sp.list[i])))
  c.r2 <- rbind(c.r2, adj.r2)
}

c.r2$V1 <- as.numeric(as.character(c.r2$V1))
c.predict$error <- abs(c.predict$corr - c.predict$fitted)
c.predict$p.error <- (c.predict$error/c.predict$fitted)

# > describe(v.predict$error)
# vars   n mean   sd median trimmed  mad min  max range skew kurtosis   se
# X1    1 148 0.69 0.74   0.45    0.57 0.45   0 4.13  4.13 2.08     5.54 0.06
# > describe(v.predict$p.error)
# vars   n mean   sd median trimmed  mad min max range skew kurtosis   se
# X1    1 148 0.34 0.47   0.24    0.26 0.21   0 3.5   3.5 4.72    26.58 0.04

# > describe(c.predict$error)
# vars  n mean   sd median trimmed  mad  min  max range skew kurtosis  se
# X1    1 25 0.77 0.52   0.75    0.76 0.71 0.08 1.66  1.58 0.18    -1.58 0.1
# > describe(c.predict$p.error)
# vars  n mean   sd median trimmed  mad  min  max range skew kurtosis   se
# X1    1 25 0.17 0.11   0.17    0.16 0.15 0.02 0.38  0.36 0.19    -1.41 0.02

spp.error <- data.frame()
veg.error <- data.frame()
sp.list <- unique(v.predict$sp)
veg.list <- unique(v.predict$veg)
for(i in 1:length(sp.list)){
  data <- v.predict[v.predict$sp == sp.list[i], ]
  a <- mean(data$error, na.rm = TRUE)
  b <- mean(data$p.error, na.rm = TRUE)
  c <- paste(as.character(sp.list[i], sep = ""))
  d <- sd(data$error, na.rm = TRUE)
  f <- sd(data$p.error, na.rm = TRUE)
  temp <- data.frame(sp = c, err = a, rel = b, err.sd = d, rel.sd = f)
  spp.error <- rbind(spp.error, temp)
}
for(i in 1:length(veg.list)){
  data <- v.predict[v.predict$veg == veg.list[i], ]
  a <- mean(data$error, na.rm = TRUE)
  b <- mean(data$p.error, na.rm = TRUE)
  c <- paste(as.character(veg.list[i], sep = ""))
  d <- sd(data$error, na.rm = TRUE)
  f <- sd(data$p.error, na.rm = TRUE)
  temp <- data.frame(sp = c, err = a, rel = b, err.sd = d, rel.sd = f)
  veg.error <- rbind(veg.error, temp)
}

spp.error <- merge(spp.error, meta, by.x = "sp", by.y = "Sound")

veg.error$sp <- factor(veg.error$sp, c("bog","conifer","deciduous","fen","grass","road"))

veg.err.plot <- ggplot(veg.error, aes(sp, err, factor = "sp")) +
  geom_point(aes(x=sp, y=err), size = 0.4) +
  geom_errorbar(aes(ymin=err-err.sd, ymax=err+err.sd), width=0.2, position=position_dodge(0.9), size = 0.3) +
  ylab("Mean prediction error") + xlab("Vegetation Type") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust=0.03, size=8)) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("A")
veg.err.plot

veg.rel.plot <- ggplot(veg.error, aes(sp, rel, factor = "sp")) +
  geom_point(aes(x=sp, y=rel), size = 0.4) +
  geom_errorbar(aes(ymin=rel-rel.sd, ymax=rel+rel.sd), width=0.2, position=position_dodge(0.9), size = 0.3) +
  ylab("Mean relative prediction error") + xlab("Vegetation Type") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust=0.03, size=8)) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("")
veg.rel.plot

spp.err.plot <- ggplot(spp.error, aes(MinFreq, err, factor = "sp")) +
  geom_point(aes(x=MinFreq, y=err), size = 0.4) +
  geom_errorbar(aes(ymin=err-err.sd, ymax=err+err.sd), width=0.2, position=position_dodge(0.9), size = 0.3) +
  geom_text_repel(aes(label = sp), segment.alpha = 0, size = 1) +
  ylab("Mean prediction error") + xlab("Frequency (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust=0.03, size=8)) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("B")
spp.err.plot

spp.rel.plot <- ggplot(spp.error, aes(MinFreq, rel, factor = "sp")) +
  geom_point(aes(x=MinFreq, y=rel), size = 0.4) +
  geom_errorbar(aes(ymin=rel-rel.sd, ymax=rel+rel.sd), width=0.2, position=position_dodge(0.9), size = 0.3) +
  geom_text_repel(aes(label = sp), segment.alpha = 0, size = 1) +
  ylab("Mean relative prediction error") + xlab("Frequency (Hz)") +
  theme(legend.position="none") +
  theme(panel.background = element_rect(fill = 'white')) +
  theme(axis.line = element_line(colour = "black", size = 0.3)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.3), axis.text.x = element_text(size=5), axis.title.x=element_text(size=6), axis.ticks.x = element_blank()) +
  theme(axis.line.y = element_line(colour = "black", size = 0.3), axis.text.y = element_text(size=5), axis.title.y=element_text(size=6)) +
  theme(legend.title=element_blank()) +
  theme(plot.title = element_text(hjust=0.03, size=8)) +
  theme(axis.ticks = element_line(colour = "black", size= 0.3)) +
  theme(legend.key.size = unit(0.3, "cm")) +
  theme(legend.text = element_text(size=6)) +
  ggtitle("")
spp.rel.plot


tiff(width = 7, height = 4, units = "in", res = 600, file = "Figure 4 - prediction error.tiff")
grid.arrange(veg.err.plot, veg.rel.plot, spp.err.plot, spp.rel.plot, ncol = 2)
dev.off()

sp.list <- unique(sepspp.sepveg.mc$sp)
v.percent <- data.frame()
for(i in 1:length(sp.list)){
  data <- sepspp.sepveg.mc[sepspp.sepveg.mc$sp == sp.list[i], ]
  ref <- data[data$veg == "edge", ]
  data$corr <- ((data$mean-ref$mean)/data$mean)*100 #reference level is always numerator, levels to be corrected are in denominator
  v.percent <- rbind(v.percent, data)
}

sp.list <- unique(sepspp.ocveg.mc$sp)
c.percent <- data.frame()
for(i in 1:length(sp.list)){
  data <- sepspp.ocveg.mc[sepspp.ocveg.mc$sp == sp.list[i], ]
  ref <- data[data$veg == "open", ]
  data$corr <- ((ref$mean-data$mean)/data$mean)*100 #reference level is always numerator, levels to be corrected are in denominator
  c.percent <- rbind(c.percent, data)
}

v.percent <- v.percent[v.percent$veg != "edge", ]
c.percent <- c.percent[c.percent$veg != "open", ]

summary.error <- data.frame()
summary.v.percent <- v.percent
summary.v.percent <- summary.v.percent[complete.cases(summary.v.percent),]
veg.list <- unique(summary.v.percent$veg)
for(i in 1:length(veg.list)){
  data <- summary.v.percent[summary.v.percent$veg == veg.list[i], ]
  a <- mean(data$corr)
  b <- sd(data$corr)
  c <- paste(as.character(veg.list[i], sep = ""))
  temp <- data.frame(veg = c, mean = a, sd = b)
  summary.error <- rbind(summary.error, temp)
}
