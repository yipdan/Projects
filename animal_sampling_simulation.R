library(mrds)
library(reshape2)
library(dplyr)
library(unmarked)
library(Distance)

setwd("C:/Users/Daniel/Desktop/Analysis/Density/")
####Set Simulation Parameters####
#f,g,h,i,j,k,l are iterator variables
#enter simulation parameters
load(file = "CONI_ddf.rda") #CONI detection function

sample_xmax = 2000 #how big is sample area on x-axis in metres
sample_ymax = 2000 #how big is sample area on y-axis in metres
dens = 0.4 #no ind./ha.
buffer = 100 #mininmum spacing between individuals (simulating territorial exclusion) in metres
s_iter = 50 #number of sampling locations to simulate
reloc = 10 #number of relocations/revisits to sample
s_length = 10 #length of survey in minutes to calculate availability
dist_fun = m #detection function (distance model containing ddf object - use require(mrds))
avail_fun = 0.7933194 #availability function - calculated from availability.R
i_dist = "uniform" #animal distribution - TO BE IMPLEMENTED

#movement parameters (see CONI.MCPandUSE.R for calculations)
r = 799.4406 #territory radius parameter
p.move.100 = 0.25 #probability of distribution within territory
p.move.75 = 0.25
p.move.50 = 0.25
p.move.25 = 0.25


#define sampling schemes
#single point
sample_scheme_singlept <- data.frame(locx = 1000, locy = 1000)
plot_singlept <- plot(sample_scheme_singlept$locx, sample_scheme_singlept$locy)

#cluster points
sample_scheme_cluster <- data.frame(station = c("NW", "SW", "NE", "SE", "CT"), 
                                    locx = c(700, 700, 1300, 1300, 1000), 
                                    locy = c(1300, 700, 1300, 700, 1000))
plot_cluster <- plot(sample_scheme_cluster$locx, sample_scheme_cluster$locy)

#grid points
sample_scheme_grid <- data.frame(station = c("1","2","3","4","5","6","7","8","9","10","11",
                                             "12","13","14","15","16"), locx = c(550, 850, 1150, 1450,
                                                                                 550, 850, 1150, 1450,
                                                                                 550, 850, 1150, 1450,
                                                                                 550, 850, 1150, 1450),
                                 locy = c(1450, 1450, 1450, 1450, 1150, 1150, 1150, 1150,
                                          850, 850, 850, 850, 550, 550, 550, 550))
plot_grid <- plot(sample_scheme_grid$locx, sample_scheme_grid$locy)

####Simulate CONI Locations####
#randomly sample density from poisson distribution
#randomly assign locations to individuals (random uniform distribution)
#buffer individuals by min separation value
convert_ha = (sample_xmax*sample_ymax)/10000
indiv <- convert_ha*dens
singlept_data <- data.frame(survey.no = numeric(), ID=character(), reloc = character(), distance = numeric(),
                            locx = numeric(), locy = numeric())
cluster_data <- data.frame(survey.no = numeric(), station = character(), ID=character(), reloc = character(),
                           distance = numeric(), locx = numeric(), locy = numeric())
grid_data <- data.frame(survey.no = numeric(), station = character(), ID=character(), reloc = character(),
                        distance = numeric(), locx = numeric(), locy = numeric())

for(h in 1:s_iter){
  sample_density <- rpois(1, indiv)
  
  sim_loc <- data.frame(survey.no = numeric(), ID = numeric(), reloc = character(), locx = numeric(), locy = numeric())
  
  sim_locx <- runif(1, 0, sample_xmax)
  sim_locy <- runif(1, 0, sample_ymax)
  sim_loc_temp <- data.frame(survey.no = h, ID = sample_density, reloc = "start", locx = sim_locx, locy = sim_locy)
  sim_loc <- rbind(sim_loc, sim_loc_temp)
  
  i = 1
  while(i < sample_density){
    sim_locx <- runif(1, 0, sample_xmax)
    sim_locy <- runif(1, 0, sample_ymax)
    sim_loc_temp <- data.frame(survey.no = h, ID = i, reloc = "start", locx = sim_locx, locy = sim_locy)
    
    dist_check <- sqrt(((abs(sim_loc$locx - sim_loc_temp$locx))^2) + 
                         (abs(sim_loc$locy - sim_loc_temp$locy))^2)
    
    if(all(dist_check > buffer)){
      sim_loc <- rbind(sim_loc, sim_loc_temp)
      i = i + 1
    }
  }
  
  #plot(sim_loc$locx, sim_loc$locy)
  
  ####Movement####
  #simulate movment, N number of CONI relocations
  for(g in 1:nrow(sim_loc)) {
    row <- sim_loc[g,]
    for(f in 1:reloc){
      ro <- sample(c(1,2,3,4), 1, replace = TRUE, prob = c(p.move.100, p.move.75, 
                                                           p.move.50, p.move.25))
      
      if(ro == 1){
        min = 578.3313
        max = 799.4406
      }
      
      if(ro == 2){
        min = 396.9598
        max = 578.3313
      }
      
      if(ro == 3){
        min = 240.4634
        max = 396.9598
      }
      
      if(ro == 4){
        min = 0
        max = 240.4634
      }
      
      rad.random <- runif(1, min = min, max = max)
      theta <- runif(1, min = 0, max = 2*pi)
      
      xdiff = sin(theta)*rad.random
      ydiff = cos(theta)*rad.random
      
      mlocx = row$locx + xdiff
      mlocy = row$locy + ydiff
      
      sim_loc_temp <- data.frame(survey.no = row$survey.no, ID = row$ID, reloc = as.factor(f), locx = mlocx, locy = mlocy)
      sim_loc <- rbind(sim_loc, sim_loc_temp)
    }
  }
  
  ####Availability####
  #function to predict availability of CONI
  for(j in 1:s_length){
    sim_loc[[as.character(paste("a", j, sep=""))]] <- apply(sim_loc, MARGIN = 1, function(x) rbinom(1, 1, avail_fun))
  }
  
  #remove any individuals who did not have detections in any of the sample minutes
  sim_loc <- sim_loc[unique(row(sim_loc[-c(1:5)])[sim_loc[-c(1:5)] == "1"]),]

                              
  #calculate distances to animals
  #single point
  distance_singlept <- sqrt(((abs(sim_loc$locx - sample_scheme_singlept$locx))^2) + 
                              (abs(sim_loc$locy - sample_scheme_singlept$locy))^2)
  singlept_data_temp <- data.frame(sim_loc, ID = as.character(paste(h, "-", sim_loc$ID, sep="")), distance = distance_singlept)
  singlept_data <- rbind(singlept_data, singlept_data_temp)
  
  #cluster points
  for(k in 1:nrow(sample_scheme_cluster)) {
    row <- sample_scheme_cluster[k,]
    distance_cluster <- sqrt(((abs(sim_loc$locx - row$locx))^2) + 
                               (abs(sim_loc$locy - row$locy))^2)
    cluster_data_temp <- data.frame(sim_loc, station = row$station, ID = as.character(paste(h, "-", sim_loc$ID, sep="")), distance = distance_cluster)
    cluster_data <- rbind(cluster_data, cluster_data_temp)
  }
  
  #grid points
  for(l in 1:nrow(sample_scheme_grid)) {
    row <- sample_scheme_grid[l,]
    distance_grid <- sqrt(((abs(sim_loc$locx - row$locx))^2) + 
                            (abs(sim_loc$locy - row$locy))^2)
    grid_data_temp <- data.frame(sim_loc, station = row$station, ID = as.character(paste(h, "-", sim_loc$ID, sep="")), distance = distance_grid)
    grid_data <- rbind(grid_data, grid_data_temp)
  }
  print(as.character(paste((h/s_iter)*100, "%", "      Simulating animal locations, movement, and availability with  ***", s_length, 
                           "***  availability blocks and  ***", reloc, "***  revisits...")))
}

singlept_data <- singlept_data[singlept_data$distance <= 600, ]
cluster_data <- cluster_data[cluster_data$distance <= 600, ]
grid_data <- grid_data[grid_data$distance <= 600, ]

####simulate error and bias of RSL for animals####
#path = df$known distance, meanE = mean error, sdE = sd of error
rtnorm <- function(path, meanE, sdE){
  aa <- rnorm(n = 1, mean = path + meanE, sd = sdE)
  
  while(any(aa<0)) { aa <- rnorm(n = 1, mean = path + meanE, sd = sdE) }
  aa
}

####simulate RSL for CONI hand measure
# Bin Min Max         Mean       SD
# 1   10 299 499 -101.4126769 66.04969
# 2    9 227 299  -16.0708441 60.68472
# 3    8 195 227    7.2520734 61.14422
# 4    7 165 195    9.0299374 53.53200
# 5    6 135 165  -10.6811944 51.37228
# 6    5 105 135   -0.1347489 46.97572
# 7    4  76 105   10.1682685 40.31009
# 8    3  47  76   20.0687247 40.93959
# 9    2  20  47   36.3447640 39.69040
# 10   1   0  20   39.8799178 33.03501

singlept_data.b1 <- singlept_data[singlept_data$distance <= 20, ]
singlept_data.b2 <- singlept_data[singlept_data$distance > 20 & singlept_data$distance <= 47, ]
singlept_data.b3 <- singlept_data[singlept_data$distance > 47 & singlept_data$distance <= 76, ]
singlept_data.b4 <- singlept_data[singlept_data$distance > 76 & singlept_data$distance <= 105, ]
singlept_data.b5 <- singlept_data[singlept_data$distance > 105 & singlept_data$distance <= 135, ]
singlept_data.b6 <- singlept_data[singlept_data$distance > 135 & singlept_data$distance <= 165, ]
singlept_data.b7 <- singlept_data[singlept_data$distance > 165 & singlept_data$distance <= 195, ]
singlept_data.b8 <- singlept_data[singlept_data$distance > 195 & singlept_data$distance <= 227, ]
singlept_data.b9 <- singlept_data[singlept_data$distance > 227 & singlept_data$distance <= 299, ]
singlept_data.b10 <- singlept_data[singlept_data$distance > 299, ]

singlept_data.pred.dist <- data.frame()

generateStdev <- function(rowSegment, meanE, sdE){
  for(k in 1:nrow(rowSegment)){
    row <- singlept_data.b1[k,]
    row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, meanE, sdE)))
    singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
    print(as.character(paste((k/nrow(singlept_data.b1))*100, "%", "      Simulating RSL for B1 (singlept.manual)...")))
  }
}
optimizr <- function(rowsSegment, meanE, sdE){
manual <- numeric (nrow(rowsSegment))
for(k in 1:nrow(rowsSegment)){
  row <- rowsSegment[k,]
  manual[k]  <- rtnorm(row$distance, meanE, sdE)

print(as.character(paste((k/nrow(singlept_data.b1))*100, "%", "      Simulating RSL for B1 (singlept.manual)...")))
}
rowsSegment$manual <- manual
singlept_data.pred.dist <- rbind(singlept_data.pred.dist, rowsSegment)
}

optimizr(singlept_data.b10, 39.8799178, 33.03501)
# for(k in 1:nrow(singlept_data.b1)){
#   row <- singlept_data.b1[k,]
#   row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 39.8799178, 33.03501)))
#   singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
#   print(as.character(paste((k/nrow(singlept_data.b1))*100, "%", "      Simulating RSL for B1 (singlept.manual)...")))
# }

generateStdev(singlept_data.b2, 36.3447640, 39.69040)
# for(k in 1:nrow(singlept_data.b2)){
#   row <- singlept_data.b2[k,]
#   row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 36.3447640, 39.69040)))
#   singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
#   print(as.character(paste((k/nrow(singlept_data.b2))*100, "%", "      Simulating RSL for B2 (singlept.manual)...")))
# }

for(k in 1:nrow(singlept_data.b3)){
  row <- singlept_data.b3[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 20.0687247, 40.93959)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b3))*100, "%", "      Simulating RSL for B3 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b4)){
  row <- singlept_data.b4[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 10.1682685, 40.31009)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b4))*100, "%", "      Simulating RSL for B4 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b5)){
  row <- singlept_data.b5[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -0.1347489, 46.97572)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b5))*100, "%", "      Simulating RSL for B5 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b6)){
  row <- singlept_data.b6[k,]
  row$singlept_data <- apply(row, 1, function(x) (rtnorm(row$distance, -10.6811944, 51.37228)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b6))*100, "%", "      Simulating RSL for B6 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b7)){
  row <- singlept_data.b7[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 9.0299374, 53.53200)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b7))*100, "%", "      Simulating RSL for B7 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b8)){
  row <- singlept_data.b8[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 7.2520734, 61.14422)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b8))*100, "%", "      Simulating RSL for B8 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b9)){
  row <- singlept_data.b9[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -16.0708441, 60.68472)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b9))*100, "%", "      Simulating RSL for B9 (singlept.manual)...")))
}

for(k in 1:nrow(singlept_data.b10)){
  row <- singlept_data.b10[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -101.4126769, 66.04969)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b10))*100, "%", "      Simulating RSL for B10 (singlept.manual)...")))
}

singlept_data <- merge(singlept_data, singlept_data.pred.dist)
####RSL simulation of CONI CNN
# Bin Min Max       Mean       SD
# 1   10 300 501 -75.894581 62.88787
# 2    9 225 300 -24.922001 51.20347
# 3    8 195 225  -5.988805 51.42895
# 4    7 160 195  10.145396 49.16111
# 5    6 130 160  19.134906 48.04585
# 6    5 102 130  11.084199 37.46374
# 7    4  75 102   7.379064 37.06011
# 8    3  47  75  14.698139 38.08312
# 9    2  22  47  14.035220 25.57165
# 10   1   0  22  33.541653 37.81236

singlept_data.b1 <- singlept_data[singlept_data$distance <= 22, ]
singlept_data.b2 <- singlept_data[singlept_data$distance > 22 & singlept_data$distance <= 47, ]
singlept_data.b3 <- singlept_data[singlept_data$distance > 47 & singlept_data$distance <= 75, ]
singlept_data.b4 <- singlept_data[singlept_data$distance > 75 & singlept_data$distance <= 102, ]
singlept_data.b5 <- singlept_data[singlept_data$distance > 102 & singlept_data$distance <= 130, ]
singlept_data.b6 <- singlept_data[singlept_data$distance > 130 & singlept_data$distance <= 160, ]
singlept_data.b7 <- singlept_data[singlept_data$distance > 160 & singlept_data$distance <= 195, ]
singlept_data.b8 <- singlept_data[singlept_data$distance > 195 & singlept_data$distance <= 225, ]
singlept_data.b9 <- singlept_data[singlept_data$distance > 225 & singlept_data$distance <= 300, ]
singlept_data.b10 <- singlept_data[singlept_data$distance > 300, ]

singlept_data.pred.dist <- data.frame()

for(k in 1:nrow(singlept_data.b1)){
  row <- singlept_data.b1[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 33.541653, 37.81236)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b1))*100, "%", "      Simulating RSL for B1 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b2)){
  row <- singlept_data.b2[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.035220, 25.57165)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b2))*100, "%", "      Simulating RSL for B2 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b3)){
  row <- singlept_data.b3[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.698139, 38.08312)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b3))*100, "%", "      Simulating RSL for B3 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b4)){
  row <- singlept_data.b4[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 7.379064, 37.06011)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b4))*100, "%", "      Simulating RSL for B4 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b5)){
  row <- singlept_data.b5[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 11.084199, 37.46374)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b5))*100, "%", "      Simulating RSL for B5 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b6)){
  row <- singlept_data.b6[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 19.134906, 48.04585)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b6))*100, "%", "      Simulating RSL for B6 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b7)){
  row <- singlept_data.b7[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 10.145396, 49.16111)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b7))*100, "%", "      Simulating RSL for B7 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b8)){
  row <- singlept_data.b8[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -5.988805, 51.42895)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b8))*100, "%", "      Simulating RSL for B8 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b9)){
  row <- singlept_data.b9[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -24.922001, 51.20347)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b9))*100, "%", "      Simulating RSL for B9 (singlept.recognizer)...")))
}

for(k in 1:nrow(singlept_data.b10)){
  row <- singlept_data.b10[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -75.894581, 62.88787)))
  singlept_data.pred.dist <- rbind(singlept_data.pred.dist, row)
  print(as.character(paste((k/nrow(singlept_data.b10))*100, "%", "      Simulating RSL for B10 (singlept.recognizer)...")))
}

singlept_data <- merge(singlept_data, singlept_data.pred.dist)

####simulate RSL for CONI hand measure
# Bin Min Max         Mean       SD
# 1   10 299 499 -101.4126769 66.04969
# 2    9 227 299  -16.0708441 60.68472
# 3    8 195 227    7.2520734 61.14422
# 4    7 165 195    9.0299374 53.53200
# 5    6 135 165  -10.6811944 51.37228
# 6    5 105 135   -0.1347489 46.97572
# 7    4  76 105   10.1682685 40.31009
# 8    3  47  76   20.0687247 40.93959
# 9    2  20  47   36.3447640 39.69040
# 10   1   0  20   39.8799178 33.03501

cluster_data.b1 <- cluster_data[cluster_data$distance <= 20, ]
cluster_data.b2 <- cluster_data[cluster_data$distance > 20 & cluster_data$distance <= 47, ]
cluster_data.b3 <- cluster_data[cluster_data$distance > 47 & cluster_data$distance <= 76, ]
cluster_data.b4 <- cluster_data[cluster_data$distance > 76 & cluster_data$distance <= 105, ]
cluster_data.b5 <- cluster_data[cluster_data$distance > 105 & cluster_data$distance <= 135, ]
cluster_data.b6 <- cluster_data[cluster_data$distance > 135 & cluster_data$distance <= 165, ]
cluster_data.b7 <- cluster_data[cluster_data$distance > 165 & cluster_data$distance <= 195, ]
cluster_data.b8 <- cluster_data[cluster_data$distance > 195 & cluster_data$distance <= 227, ]
cluster_data.b9 <- cluster_data[cluster_data$distance > 227 & cluster_data$distance <= 299, ]
cluster_data.b10 <- cluster_data[cluster_data$distance > 299, ]

cluster_data.pred.dist <- data.frame()

for(k in 1:nrow(cluster_data.b1)){
  row <- cluster_data.b1[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 39.8799178, 33.03501)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b1))*100, "%", "      Simulating RSL for B1 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b2)){
  row <- cluster_data.b2[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 36.3447640, 39.69040)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b2))*100, "%", "      Simulating RSL for B2 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b3)){
  row <- cluster_data.b3[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 20.0687247, 40.93959)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b3))*100, "%", "      Simulating RSL for B3 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b4)){
  row <- cluster_data.b4[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 10.1682685, 40.31009)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b4))*100, "%", "      Simulating RSL for B4 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b5)){
  row <- cluster_data.b5[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -0.1347489, 46.97572)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b5))*100, "%", "      Simulating RSL for B5 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b6)){
  row <- cluster_data.b6[k,]
  row$cluster_data <- apply(row, 1, function(x) (rtnorm(row$distance, -10.6811944, 51.37228)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b6))*100, "%", "      Simulating RSL for B6 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b7)){
  row <- cluster_data.b7[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 9.0299374, 53.53200)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b7))*100, "%", "      Simulating RSL for B7 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b8)){
  row <- cluster_data.b8[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 7.2520734, 61.14422)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b8))*100, "%", "      Simulating RSL for B8 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b9)){
  row <- cluster_data.b9[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -16.0708441, 60.68472)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b9))*100, "%", "      Simulating RSL for B9 (cluster.manual)...")))
}

for(k in 1:nrow(cluster_data.b10)){
  row <- cluster_data.b10[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -101.4126769, 66.04969)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b10))*100, "%", "      Simulating RSL for B10 (cluster.manual)...")))
}

cluster_data <- merge(cluster_data, cluster_data.pred.dist)
####RSL simulation of CONI CNN
# Bin Min Max       Mean       SD
# 1   10 300 501 -75.894581 62.88787
# 2    9 225 300 -24.922001 51.20347
# 3    8 195 225  -5.988805 51.42895
# 4    7 160 195  10.145396 49.16111
# 5    6 130 160  19.134906 48.04585
# 6    5 102 130  11.084199 37.46374
# 7    4  75 102   7.379064 37.06011
# 8    3  47  75  14.698139 38.08312
# 9    2  22  47  14.035220 25.57165
# 10   1   0  22  33.541653 37.81236

cluster_data.b1 <- cluster_data[cluster_data$distance <= 22, ]
cluster_data.b2 <- cluster_data[cluster_data$distance > 22 & cluster_data$distance <= 47, ]
cluster_data.b3 <- cluster_data[cluster_data$distance > 47 & cluster_data$distance <= 75, ]
cluster_data.b4 <- cluster_data[cluster_data$distance > 75 & cluster_data$distance <= 102, ]
cluster_data.b5 <- cluster_data[cluster_data$distance > 102 & cluster_data$distance <= 130, ]
cluster_data.b6 <- cluster_data[cluster_data$distance > 130 & cluster_data$distance <= 160, ]
cluster_data.b7 <- cluster_data[cluster_data$distance > 160 & cluster_data$distance <= 195, ]
cluster_data.b8 <- cluster_data[cluster_data$distance > 195 & cluster_data$distance <= 225, ]
cluster_data.b9 <- cluster_data[cluster_data$distance > 225 & cluster_data$distance <= 300, ]
cluster_data.b10 <- cluster_data[cluster_data$distance > 300, ]

cluster_data.pred.dist <- data.frame()

for(k in 1:nrow(cluster_data.b1)){
  row <- cluster_data.b1[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 33.541653, 37.81236)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b1))*100, "%", "      Simulating RSL for B1 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b2)){
  row <- cluster_data.b2[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.035220, 25.57165)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b2))*100, "%", "      Simulating RSL for B2 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b3)){
  row <- cluster_data.b3[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.698139, 38.08312)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b3))*100, "%", "      Simulating RSL for B3 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b4)){
  row <- cluster_data.b4[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 7.379064, 37.06011)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b4))*100, "%", "      Simulating RSL for B4 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b5)){
  row <- cluster_data.b5[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 11.084199, 37.46374)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b5))*100, "%", "      Simulating RSL for B5 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b6)){
  row <- cluster_data.b6[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 19.134906, 48.04585)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b6))*100, "%", "      Simulating RSL for B6 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b7)){
  row <- cluster_data.b7[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 10.145396, 49.16111)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b7))*100, "%", "      Simulating RSL for B7 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b8)){
  row <- cluster_data.b8[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -5.988805, 51.42895)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b8))*100, "%", "      Simulating RSL for B8 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b9)){
  row <- cluster_data.b9[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -24.922001, 51.20347)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b9))*100, "%", "      Simulating RSL for B9 (cluster.recognizer)...")))
}

for(k in 1:nrow(cluster_data.b10)){
  row <- cluster_data.b10[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -75.894581, 62.88787)))
  cluster_data.pred.dist <- rbind(cluster_data.pred.dist, row)
  print(as.character(paste((k/nrow(cluster_data.b10))*100, "%", "      Simulating RSL for B10 (cluster.recognizer)...")))
}

cluster_data <- merge(cluster_data, cluster_data.pred.dist)

####simulate RSL for CONI hand measure
# Bin Min Max         Mean       SD
# 1   10 299 499 -101.4126769 66.04969
# 2    9 227 299  -16.0708441 60.68472
# 3    8 195 227    7.2520734 61.14422
# 4    7 165 195    9.0299374 53.53200
# 5    6 135 165  -10.6811944 51.37228
# 6    5 105 135   -0.1347489 46.97572
# 7    4  76 105   10.1682685 40.31009
# 8    3  47  76   20.0687247 40.93959
# 9    2  20  47   36.3447640 39.69040
# 10   1   0  20   39.8799178 33.03501

grid_data.b1 <- grid_data[grid_data$distance <= 20, ]
grid_data.b2 <- grid_data[grid_data$distance > 20 & grid_data$distance <= 47, ]
grid_data.b3 <- grid_data[grid_data$distance > 47 & grid_data$distance <= 76, ]
grid_data.b4 <- grid_data[grid_data$distance > 76 & grid_data$distance <= 105, ]
grid_data.b5 <- grid_data[grid_data$distance > 105 & grid_data$distance <= 135, ]
grid_data.b6 <- grid_data[grid_data$distance > 135 & grid_data$distance <= 165, ]
grid_data.b7 <- grid_data[grid_data$distance > 165 & grid_data$distance <= 195, ]
grid_data.b8 <- grid_data[grid_data$distance > 195 & grid_data$distance <= 227, ]
grid_data.b9 <- grid_data[grid_data$distance > 227 & grid_data$distance <= 299, ]
grid_data.b10 <- grid_data[grid_data$distance > 299, ]

grid_data.pred.dist <- data.frame()

for(k in 1:nrow(grid_data.b1)){
  row <- grid_data.b1[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 39.8799178, 33.03501)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b1))*100, "%", "      Simulating RSL for B1 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b2)){
  row <- grid_data.b2[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 36.3447640, 39.69040)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b2))*100, "%", "      Simulating RSL for B2 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b3)){
  row <- grid_data.b3[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 20.0687247, 40.93959)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b3))*100, "%", "      Simulating RSL for B3 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b4)){
  row <- grid_data.b4[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 10.1682685, 40.31009)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b4))*100, "%", "      Simulating RSL for B4 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b5)){
  row <- grid_data.b5[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -0.1347489, 46.97572)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b5))*100, "%", "      Simulating RSL for B5 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b6)){
  row <- grid_data.b6[k,]
  row$grid_data <- apply(row, 1, function(x) (rtnorm(row$distance, -10.6811944, 51.37228)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b6))*100, "%", "      Simulating RSL for B6 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b7)){
  row <- grid_data.b7[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 9.0299374, 53.53200)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b7))*100, "%", "      Simulating RSL for B7 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b8)){
  row <- grid_data.b8[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, 7.2520734, 61.14422)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b8))*100, "%", "      Simulating RSL for B8 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b9)){
  row <- grid_data.b9[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -16.0708441, 60.68472)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b9))*100, "%", "      Simulating RSL for B9 (grid.manual)...")))
}

for(k in 1:nrow(grid_data.b10)){
  row <- grid_data.b10[k,]
  row$manual <- apply(row, 1, function(x) (rtnorm(row$distance, -101.4126769, 66.04969)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b10))*100, "%", "      Simulating RSL for B10 (grid.manual)...")))
}

grid_data <- merge(grid_data, grid_data.pred.dist)
####RSL simulation of CONI CNN
# Bin Min Max       Mean       SD
# 1   10 300 501 -75.894581 62.88787
# 2    9 225 300 -24.922001 51.20347
# 3    8 195 225  -5.988805 51.42895
# 4    7 160 195  10.145396 49.16111
# 5    6 130 160  19.134906 48.04585
# 6    5 102 130  11.084199 37.46374
# 7    4  75 102   7.379064 37.06011
# 8    3  47  75  14.698139 38.08312
# 9    2  22  47  14.035220 25.57165
# 10   1   0  22  33.541653 37.81236

grid_data.b1 <- grid_data[grid_data$distance <= 22, ]
grid_data.b2 <- grid_data[grid_data$distance > 22 & grid_data$distance <= 47, ]
grid_data.b3 <- grid_data[grid_data$distance > 47 & grid_data$distance <= 75, ]
grid_data.b4 <- grid_data[grid_data$distance > 75 & grid_data$distance <= 102, ]
grid_data.b5 <- grid_data[grid_data$distance > 102 & grid_data$distance <= 130, ]
grid_data.b6 <- grid_data[grid_data$distance > 130 & grid_data$distance <= 160, ]
grid_data.b7 <- grid_data[grid_data$distance > 160 & grid_data$distance <= 195, ]
grid_data.b8 <- grid_data[grid_data$distance > 195 & grid_data$distance <= 225, ]
grid_data.b9 <- grid_data[grid_data$distance > 225 & grid_data$distance <= 300, ]
grid_data.b10 <- grid_data[grid_data$distance > 300, ]

grid_data.pred.dist <- data.frame()

for(k in 1:nrow(grid_data.b1)){
  row <- grid_data.b1[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 33.541653, 37.81236)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b1))*100, "%", "      Simulating RSL for B1 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b2)){
  row <- grid_data.b2[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.035220, 25.57165)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b2))*100, "%", "      Simulating RSL for B2 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b3)){
  row <- grid_data.b3[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 14.698139, 38.08312)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b3))*100, "%", "      Simulating RSL for B3 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b4)){
  row <- grid_data.b4[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 7.379064, 37.06011)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b4))*100, "%", "      Simulating RSL for B4 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b5)){
  row <- grid_data.b5[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 11.084199, 37.46374)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b5))*100, "%", "      Simulating RSL for B5 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b6)){
  row <- grid_data.b6[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 19.134906, 48.04585)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b6))*100, "%", "      Simulating RSL for B6 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b7)){
  row <- grid_data.b7[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, 10.145396, 49.16111)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b7))*100, "%", "      Simulating RSL for B7 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b8)){
  row <- grid_data.b8[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -5.988805, 51.42895)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b8))*100, "%", "      Simulating RSL for B8 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b9)){
  row <- grid_data.b9[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -24.922001, 51.20347)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b9))*100, "%", "      Simulating RSL for B9 (grid.recognizer)...")))
}

for(k in 1:nrow(grid_data.b10)){
  row <- grid_data.b10[k,]
  row$recognizer <- apply(row, 1, function(x) (rtnorm(row$distance, -75.894581, 62.88787)))
  grid_data.pred.dist <- rbind(grid_data.pred.dist, row)
  print(as.character(paste((k/nrow(grid_data.b10))*100, "%", "      Simulating RSL for B10 (grid.recognizer)...")))
}

grid_data <- merge(grid_data, grid_data.pred.dist)

load(file = "singlept_data.rda")
load(file = "cluster_data.rda")
load(file = "grid_data.rda")


####Detectability#### 
#function to predict probability of detection ~ distance using mrds conventional distance sampling model, assign detected/not detected
predictDistFunction <- function(d, xy)
  pmax(0, splinefun(xy)(d))


getDistFunction <- function (x, which = 2, breaks = NULL, nc = NULL, jitter.v = rep(0,
                                                                                   3), showpoints = TRUE, subset = NULL, pl.col = "black", bw.col = grey(0),
                            black.white = FALSE, pl.den = rep(20, 1), pl.ang = rep(-45,
                                                                                   1), main = NULL, pages = 0, pdf = FALSE, ylim = NULL,
                            xlab = "Distance", ...)
{
  model <- x$ddf
  lower <- 0
  vname <- "distance"
  dat <- model$data
  if (pdf & !model$meta.data$point) {
    warning("Ignoring pdf=TRUE for line transect data")
    pdf <- FALSE
  }
  show <- rep(FALSE, 2)
  show[which] <- TRUE
  if (black.white) {
    byval1 <- bw.col[1]
  }
  else {
    byval1 <- pl.col[1]
  }
  denval1 <- pl.den[1]
  angval1 <- pl.ang[1]
  width <- model$meta.data$width
  left <- model$meta.data$left
  ddfobj <- model$ds$aux$ddfobj
  point <- model$ds$aux$point
  if (is.null(model$ds$aux$int.range)) {
    int.range <- c(0, width)
  }
  else {
    int.range <- model$ds$aux$int.range
  }
  if (is.matrix(int.range)) {
    max.range <- as.vector(int.range[1, ])
    int.range <- int.range[2:dim(int.range)[1], ]
    range.varies <- TRUE
  }
  else {
    max.range <- int.range
    normalize <- FALSE
    range.varies <- FALSE
  }
  if (range.varies & showpoints) {
    warning("Point values can be misleading for g(x) when the range varies")
  }
  if (!is.null(substitute(subset))) {
    selected <- eval(substitute(subset), ddfobj$xmat)
  }
  else {
    selected <- rep(TRUE, nrow(ddfobj$xmat))
  }
  if (all(!selected)) {
    stop("Specified subset is empty.")
  }
  if (is.matrix(int.range)) {
    int.range <- int.range[selected, ]
  }
  xmat <- ddfobj$xmat[selected, ]
  if (!is.null(ddfobj$scale)) {
    z <- ddfobj$scale$dm[selected, , drop = FALSE]
  }
  else {
    z <- matrix(1, nrow = 1, ncol = 1)
  }
  if (length(model$fitted) == 1) {
    pdot <- rep(model$fitted, sum(as.numeric(selected)))
  }
  else {
    pdot <- model$fitted[selected]
    Nhat <- sum(1/pdot)
  }
  zdim <- dim(z)[2]
  n <- length(xmat$distance)
  if (!is.null(breaks)) {
    nc <- length(breaks) - 1
  }
  if (is.null(nc)) {
    nc <- round(sqrt(n), 0)
  }
  hascov <- FALSE
  if (!ddfobj$intercept.only) {
    hascov <- TRUE
  }
  if (!hascov) {
    xgrid <- seq(0, width, length.out = 101)
    zgrid <- matrix(rep(z[1, ], length(xgrid)), byrow = TRUE,
                    ncol = sum(zdim))
  }
  if (is.null(breaks)) {
    if (is.null(model$meta.data$binned)) {
      binned <- FALSE
    }
    else {
      binned <- model$meta.data$binned
    }
    if (binned) {
      breaks <- model$ds$aux$breaks
      nc <- length(breaks) - 1
    }
    else {
      breaks <- c(max(0, (max.range[1])), max.range[1] +
                    ((max.range[2] - max.range[1])/nc) * (1:nc))
      if (breaks[1] > left) {
        breaks <- c(left, breaks)
        nc <- nc + 1
      }
    }
  }
  breaks <- mrds:::test.breaks(breaks, model$meta.data$left, width)
  nc <- length(breaks) - 1
  lower <- min(breaks)
  upper <- max(breaks)
  dat <- dat[selected, ]
  keep <- dat[, vname] >= lower & dat[, vname] <= upper
  hist.obj <- hist(dat[, vname][keep], breaks = breaks, plot = FALSE)
  ymax <- max(hist.obj$counts)
  if (normalize & !point) {
    bindata <- function(x, r, breaks) {
      return(hist(r[r >= x[1] & r <= x[2]], breaks = breaks,
                  plot = FALSE)$counts)
    }
    sumit <- function(x, n, wt) {
      return(sum(x/(wt * n)))
    }
    expected.counts <- apply(int.range, 1, bindata, r = (0:1000) *
                               width/1001, breaks = breaks)
    expected.counts <- apply(expected.counts, 1, sumit, n = 1001,
                             wt = pdot)
  }
  else {
    if (!point) {
      expected.counts <- (breaks[2:(nc + 1)] - breaks[1:nc]) *
        (Nhat/breaks[nc + 1])
    }
    else {
      if (!pdf) {
        expected.counts <- -apply(matrix(c(breaks[2:(nc +
                                                       1)]^2, breaks[1:nc]^2), ncol = 2, nrow = nc),
                                  1, diff) * (Nhat/breaks[nc + 1]^2)
      }
      else {
        expected.counts <- sum(hist.obj$counts)
      }
    }
  }
  if (!(pdf & point)) {
    hist.obj$density <- hist.obj$counts/expected.counts
    hist.obj$density[expected.counts == 0] <- 0
  }
  hist.obj$equidist <- FALSE
  if (pages != 1 & sum(show) != 1) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  else if (sum(show) != 1) {
    opar <- par(mfrow = c(1, sum(show)))
    on.exit(par(opar))
  }
  if (show[1]) {
    if (is.null(ylim))
      ylim <- c(0, ymax)
    mrds:::histline(hist.obj$counts, breaks = breaks, lineonly = FALSE,
                    ylim = ylim, xlab = xlab, ylab = "Frequency", angle = angval1,
                    density = denval1, col = byval1, ...)
  }
  if (show[2]) {
    hist_area <- sum(hist.obj$density * diff(breaks))
    if (point & pdf) {
      point_vals <- distpdf(xmat$distance, ddfobj, width = width,
                            point = TRUE, standardize = TRUE)/integratepdf(ddfobj,
                                                                           select = selected, width = width, int.range = int.range,
                                                                           standardize = TRUE, point = TRUE)
    }
    else {
      point_vals <- detfct(xmat$distance, ddfobj, select = selected,
                           width = width)
    }
    if (is.null(ylim))
      ylim <- c(0, max(hist.obj$density, max(point_vals)))
    if (pdf) {
      ylab <- "Probability density"
      det.plot <- FALSE
    }
    else {
      ylab <- "Detection probability"
      det.plot <- TRUE
    }
    mrds:::histline(hist.obj$density, breaks = breaks, lineonly = FALSE,
                    xlab = xlab, ylab = ylab, ylim = ylim, angle = angval1,
                    density = denval1, col = byval1, det.plot = det.plot,
                    ...)
    if (hascov) {
      finebr <- seq(0, width, length.out = 101)
      xgrid <- NULL
      linevalues <- NULL
      newdat <- xmat
      for (i in 1:(length(finebr) - 1)) {
        x <- (finebr[i] + finebr[i + 1])/2
        xgrid <- c(xgrid, x)
        newdat$distance <- rep(x, nrow(newdat))
        detfct.values <- detfct(newdat$distance, ddfobj,
                                select = selected, width = width)
        if (!normalize & range.varies) {
          detfct.values[x < int.range[, 1] | x > int.range[,
                                                           2]] <- 0
        }
        if (point & pdf) {
          r_gr <- distpdf(newdat$distance, ddfobj, width = width,
                          point = TRUE, standardize = TRUE)
          int_r_gr <- integratepdf(ddfobj, select = selected,
                                   width = width, int.range = int.range, standardize = TRUE,
                                   point = TRUE)
          pdf_vals <- r_gr/int_r_gr
          vals <- pdf_vals * hist_area
          linevalues <- c(linevalues, sum(vals/pdot)/sum(1/pdot))
        }
        else {
          linevalues <- c(linevalues, sum(detfct.values/pdot)/sum(1/pdot))
        }
      }
    }
    else {
      if (!is.null(ddfobj$scale)) {
        ddfobj$scale$dm <- ddfobj$scale$dm[rep(1, length(xgrid)),
                                           , drop = FALSE]
      }
      if (!is.null(ddfobj$shape)) {
        ddfobj$shape$dm <- ddfobj$shape$dm[rep(1, length(xgrid)),
                                           , drop = FALSE]
      }
      if (ddfobj$type == "gamma") {
        xgrid[1] <- xgrid[1] + 1e-06
      }
      if (point & pdf) {
        r_gr <- distpdf(xgrid, ddfobj, width = width,
                        point = TRUE, standardize = TRUE)
        int_r_gr <- integratepdf(ddfobj, select = TRUE,
                                 width = width, int.range = int.range, standardize = TRUE,
                                 point = TRUE)[1]
        pdf_vals <- r_gr/int_r_gr
        linevalues <- pdf_vals * hist_area
      }
      else {
        linevalues <- detfct(xgrid, ddfobj, width = width)
      }
    }
    lines(xgrid, linevalues, col = byval1, ...)
    if (showpoints) {
      jitter.p <- rnorm(length(point_vals), 1, jitter.v[1])
      points(xmat$distance, point_vals * jitter.p, col = byval1,
             ...)
    }
    if (!is.null(main)) {
      title(main, cex.main = 0.8)
    }
  }
  invisible(cbind(xgrid=xgrid, linevalues=linevalues))
}

xy <- getDistFunction(m)

Dist_singlept <- singlept_data$distance
singlept_data$probability <- p <- predictDistFunction(Dist_singlept, xy)
singlept_data$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))


Dist_cluster <- cluster_data$distance
cluster_data$probability <- p <- predictDistFunction(Dist_cluster, xy)
cluster_data$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))


Dist_grid <- grid_data$distance
grid_data$probability <- p <- predictDistFunction(Dist_grid, xy)
grid_data$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

#generate datasets of detections
detections_singlept <- singlept_data[singlept_data$detected == 1, ]
detections_cluster <- cluster_data[cluster_data$detected == 1, ]
detections_grid <- grid_data[grid_data$detected == 1, ]






#********************************#
#********************************#
####DISTANCE SAMPLING####
#********************************#
#********************************#
recognizer_singlept <- singlept_data
recognizer_dist <- singlept_data$recognizer
recognizer_singlept$probability <- p <- predictDistFunction(recognizer_dist, xy)
recognizer_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

manual_singlept <- singlept_data
manual_dist <- singlept_data$manual
manual_singlept$probability <- p <- predictDistFunction(manual_dist, xy)
manual_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

detections_singlept_r <- recognizer_singlept[recognizer_singlept$detected == 1, ]
detections_singlept_m <- manual_singlept[manual_singlept$detected == 1, ]

ds_sim_r <- detections_singlept_r
ds_sim_m <- detections_singlept_m
ds_sim_k <- detections_singlept

ds_sim_r$distance <- ds_sim_r$recognizer
ds_sim_m$distance <- ds_sim_m$manual

ds_sim_r <- ds_sim_r[ds_sim_r$reloc == "start", ]
ds_sim_m <- ds_sim_m[ds_sim_m$reloc == "start", ]
ds_sim_k <- ds_sim_k[ds_sim_k$reloc == "start", ]

ds_sim_k["Sample.Label"] <- NA
ds_sim_k$Sample.Label <- ds_sim_k$survey.no
ds_sim_k["Effort"] <- 1
ds_sim_k["Region.Label"] <- "Sim"
ds_sim_k["Area"] <- 0 ####in km2, 13.0072ha, 130072 m2

ds_sim_r["Sample.Label"] <- NA
ds_sim_r$Sample.Label <- ds_sim_r$survey.no
ds_sim_r["Effort"] <- 1
ds_sim_r["Region.Label"] <- "Sim"
ds_sim_r["Area"] <- 0 ####in km2, 13.0072ha, 130072 m2

ds_sim_m["Sample.Label"] <- NA
ds_sim_m$Sample.Label <- ds_sim_m$survey.no
ds_sim_m["Effort"] <- 1
ds_sim_m["Region.Label"] <- "Sim"
ds_sim_m["Area"] <- 0 ####in km2, 13.0072ha, 130072 m2

ds_sim_r1 <- ds_sim_r[unique(row(ds_sim_r[-c(1:5, 7:25)])[ds_sim_r[-c(1:5, 7:25)] == "1"]),]
ds_sim_m1 <- ds_sim_m[unique(row(ds_sim_m[-c(1:5, 7:25)])[ds_sim_m[-c(1:5, 7:25)] == "1"]),]
ds_sim_k1 <- ds_sim_k[unique(row(ds_sim_k[-c(1:5, 7:25)])[ds_sim_k[-c(1:5, 7:25)] == "1"]),]

ds_sim_r2 <- ds_sim_r[unique(row(ds_sim_r[-c(1:5, 8:25)])[ds_sim_r[-c(1:5, 8:25)] == "1"]),]
ds_sim_m2 <- ds_sim_m[unique(row(ds_sim_m[-c(1:5, 8:25)])[ds_sim_m[-c(1:5, 8:25)] == "1"]),]
ds_sim_k2 <- ds_sim_k[unique(row(ds_sim_k[-c(1:5, 8:25)])[ds_sim_k[-c(1:5, 8:25)] == "1"]),]

ds_sim_r3 <- ds_sim_r[unique(row(ds_sim_r[-c(1:5, 9:25)])[ds_sim_r[-c(1:5, 9:25)] == "1"]),]
ds_sim_m3 <- ds_sim_m[unique(row(ds_sim_m[-c(1:5, 9:25)])[ds_sim_m[-c(1:5, 9:25)] == "1"]),]
ds_sim_k3 <- ds_sim_k[unique(row(ds_sim_k[-c(1:5, 9:25)])[ds_sim_k[-c(1:5, 9:25)] == "1"]),]

ds_sim_r5 <- ds_sim_r[unique(row(ds_sim_r[-c(1:5, 11:25)])[ds_sim_r[-c(1:5, 11:25)] == "1"]),]
ds_sim_m5 <- ds_sim_m[unique(row(ds_sim_m[-c(1:5, 11:25)])[ds_sim_m[-c(1:5, 11:25)] == "1"]),]
ds_sim_k5 <- ds_sim_k[unique(row(ds_sim_k[-c(1:5, 11:25)])[ds_sim_k[-c(1:5, 11:25)] == "1"]),]

ds_sim_r10 <- ds_sim_r
ds_sim_m10 <- ds_sim_m
ds_sim_k10 <- ds_sim_k

#KNOWN
known.hn1 <- ds(ds_sim_k, key="hn", adjustment="cos", transect="point", convert.units=0.01)
known.hn2 <- ds(ds_sim_k, key="hn", adjustment="herm", transect="point", convert.units=0.01)
known.hn3 <- ds(ds_sim_k, key="hn", adjustment="poly", transect="point", convert.units=0.01)
known.unifcos <- ds(ds_sim_k, key="unif", adjustment="cos", mono="strict", transect="point", convert.units=0.001)
known.hazard1 <- ds(ds_sim_k, key="hr", adjustment="cos", transect="point", convert.units=0.01)
known.hazard2 <- ds(ds_sim_k, key="hr", adjustment="poly", transect="point", convert.units=0.01)
#best model known.hn1 (half-normal key function, cosine(2,3) adjustments)
known.10 <- known.hn1
known.5 <- ds(ds_sim_k5, key="hn", adjustment="cos", transect="point", convert.units=0.01)
known.3 <- ds(ds_sim_k3, key="hn", adjustment="cos", transect="point", convert.units=0.01)
known.2 <- ds(ds_sim_k2, key="hn", adjustment="cos", transect="point", convert.units=0.01)
known.1 <- ds(ds_sim_k1, key="hn", adjustment="cos", transect="point", convert.units=0.01)

#MANUAL
m.hn1 <- ds(ds_sim_m, key="hn", adjustment="cos", transect="point", convert.units=0.01)
m.hn2 <- ds(ds_sim_m, key="hn", adjustment="herm", transect="point", convert.units=0.01)
m.hn3 <- ds(ds_sim_m, key="hn", adjustment="poly", transect="point", convert.units=0.01)
m.unifcos <- ds(ds_sim_m, key="unif", adjustment="cos", mono="strict", transect="point", convert.units=0.001)
m.hazard1 <- ds(ds_sim_m, key="hr", adjustment="cos", transect="point", convert.units=0.01)
m.hazard2 <- ds(ds_sim_m, key="hr", adjustment="poly", transect="point", convert.units=0.01)
#best model m.hn1 (half-normal key function, cosine(2,3) adjustments)

manual.10 <- m.hn1
manual.5 <- ds(ds_sim_m5, key="hn", adjustment="cos", transect="point", convert.units=0.01)
manual.3 <- ds(ds_sim_m3, key="hn", adjustment="cos", transect="point", convert.units=0.01)
manual.2 <- ds(ds_sim_m2, key="hn", adjustment="cos", transect="point", convert.units=0.01)
manual.1 <- ds(ds_sim_m1, key="hn", adjustment="cos", transect="point", convert.units=0.01)

#RECOGNIZER
r.hn1 <- ds(ds_sim_r, key="hn", adjustment="cos", transect="point", convert.units=0.01)
r.hn2 <- ds(ds_sim_r, key="hn", adjustment="herm", transect="point", convert.units=0.01)
r.hn3 <- ds(ds_sim_r, key="hn", adjustment="poly", transect="point", convert.units=0.01)
r.unifcos <- ds(ds_sim_r, key="unif", adjustment="cos", mono="strict", transect="point", convert.units=0.001)
r.hazard1 <- ds(ds_sim_r, key="hr", adjustment="cos", transect="point", convert.units=0.01)
r.hazard2 <- ds(ds_sim_r, key="hr", adjustment="poly", transect="point", convert.units=0.01)
#best model m.hn1 (half-normal key function, cosine(2,3) adjustments)

recognizer.10 <- r.hn1
recognizer.5 <- ds(ds_sim_r5, key="hn", adjustment="cos", transect="point", convert.units=0.01)
recognizer.3 <- ds(ds_sim_r3, key="hn", adjustment="cos", transect="point", convert.units=0.01)
recognizer.2 <- ds(ds_sim_r2, key="hn", adjustment="cos", transect="point", convert.units=0.01)
recognizer.1 <- ds(ds_sim_r1, key="hn", adjustment="cos", transect="point", convert.units=0.01)

candidate.ds <- list()
candidate.ds[["known.1"]] <- known.1
candidate.ds[["known.2"]] <- known.2
candidate.ds[["known.3"]] <- known.3
candidate.ds[["known.5"]] <- known.5
candidate.ds[["known.10"]] <- known.10
candidate.ds[["manual.1"]] <- manual.1
candidate.ds[["manual.2"]] <- manual.2
candidate.ds[["manual.3"]] <- manual.3
candidate.ds[["manual.5"]] <- manual.5
candidate.ds[["manual.10"]] <- manual.10
candidate.ds[["recognizer.1"]] <- recognizer.1
candidate.ds[["recognizer.2"]] <- recognizer.2
candidate.ds[["recognizer.3"]] <- recognizer.3
candidate.ds[["recognizer.5"]] <- recognizer.5
candidate.ds[["recognizer.10"]] <- recognizer.10

results.ds <- data.frame()
f=1
for(i in candidate.ds){
  est <- i$dht$individuals$D$Estimate[1]
  lcl <- i$dht$individuals$D$lcl[1]
  ucl <- i$dht$individuals$D$ucl[1]
  se <- i$dht$individuals$D$se[1]
  df <- i$dht$individuals$D$df[1]
  listname <- names(candidate.ds[f])
  slen <- strsplit(listname, "[.]")[[1]][2]
  ms <- strsplit(listname, "[.]")[[1]][1]
  
  results.ds.temp <- data.frame(method = "distance", measure = ms, length = slen, estimate = est,
                                lower = lcl, upper = ucl, se = se, df = df)
  results.ds <- rbind(results.ds, results.ds.temp)
  f=f+1
}
rownames(results.ds) <- c()



#********************************#
#********************************#
####N-MIXTURE MODELS####
#********************************#
#********************************#
recognizer_singlept <- singlept_data
recognizer_dist <- singlept_data$recognizer
recognizer_singlept$probability <- p <- predictDistFunction(recognizer_dist, xy)
recognizer_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

manual_singlept <- singlept_data
manual_dist <- singlept_data$manual
manual_singlept$probability <- p <- predictDistFunction(manual_dist, xy)
manual_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

detections_singlept_r <- recognizer_singlept[recognizer_singlept$detected == 1, ]
detections_singlept_m <- manual_singlept[manual_singlept$detected == 1, ]

nmix_sim_r <- detections_singlept_r
nmix_sim_m <- detections_singlept_m
nmix_sim_k <- detections_singlept

nmix_sim_r$distance <- nmix_sim_r$recognizer
nmix_sim_m$distance <- nmix_sim_m$manual

nmix_sim_r1 <- nmix_sim_r[unique(row(nmix_sim_r[-c(1:5, 7:21)])[nmix_sim_r[-c(1:5, 7:21)] == "1"]),]
nmix_sim_m1 <- nmix_sim_m[unique(row(nmix_sim_m[-c(1:5, 7:21)])[nmix_sim_m[-c(1:5, 7:21)] == "1"]),]
nmix_sim_k1 <- nmix_sim_k[unique(row(nmix_sim_k[-c(1:5, 7:21)])[nmix_sim_k[-c(1:5, 7:21)] == "1"]),]

nmix_sim_r2 <- nmix_sim_r[unique(row(nmix_sim_r[-c(1:5, 8:21)])[nmix_sim_r[-c(1:5, 8:21)] == "1"]),]
nmix_sim_m2 <- nmix_sim_m[unique(row(nmix_sim_m[-c(1:5, 8:21)])[nmix_sim_m[-c(1:5, 8:21)] == "1"]),]
nmix_sim_k2 <- nmix_sim_k[unique(row(nmix_sim_k[-c(1:5, 8:21)])[nmix_sim_k[-c(1:5, 8:21)] == "1"]),]

nmix_sim_r3 <- nmix_sim_r[unique(row(nmix_sim_r[-c(1:5, 9:21)])[nmix_sim_r[-c(1:5, 9:21)] == "1"]),]
nmix_sim_m3 <- nmix_sim_m[unique(row(nmix_sim_m[-c(1:5, 9:21)])[nmix_sim_m[-c(1:5, 9:21)] == "1"]),]
nmix_sim_k3 <- nmix_sim_k[unique(row(nmix_sim_k[-c(1:5, 9:21)])[nmix_sim_k[-c(1:5, 9:21)] == "1"]),]

nmix_sim_r5 <- nmix_sim_r[unique(row(nmix_sim_r[-c(1:5, 11:21)])[nmix_sim_r[-c(1:5, 11:21)] == "1"]),]
nmix_sim_m5 <- nmix_sim_m[unique(row(nmix_sim_m[-c(1:5, 11:21)])[nmix_sim_m[-c(1:5, 11:21)] == "1"]),]
nmix_sim_k5 <- nmix_sim_k[unique(row(nmix_sim_k[-c(1:5, 11:21)])[nmix_sim_k[-c(1:5, 11:21)] == "1"]),]

nmix_sim_r10 <- nmix_sim_r
nmix_sim_m10 <- nmix_sim_m
nmix_sim_k10 <- nmix_sim_k

nmix_sim_r1_2 <- nmix_sim_r1[(nmix_sim_r1$reloc == "start")|(nmix_sim_r1$reloc == "1"),]
nmix_sim_m1_2 <- nmix_sim_m1[(nmix_sim_m1$reloc == "start")|(nmix_sim_m1$reloc == "1"),]
nmix_sim_k1_2 <- nmix_sim_k1[(nmix_sim_k1$reloc == "start")|(nmix_sim_k1$reloc == "1"),]

nmix_sim_r1_3 <- nmix_sim_r1[(nmix_sim_r1$reloc == "start")|(nmix_sim_r1$reloc == "1")|(nmix_sim_r1$reloc == "2"),]
nmix_sim_m1_3 <- nmix_sim_m1[(nmix_sim_m1$reloc == "start")|(nmix_sim_m1$reloc == "1")|(nmix_sim_m1$reloc == "2"),]
nmix_sim_k1_3 <- nmix_sim_k1[(nmix_sim_k1$reloc == "start")|(nmix_sim_k1$reloc == "1")|(nmix_sim_k1$reloc == "2"),]

nmix_sim_r1_4 <- nmix_sim_r1[(nmix_sim_r1$reloc == "start")|(nmix_sim_r1$reloc == "1")|(nmix_sim_r1$reloc == "2")|(nmix_sim_r1$reloc == "3"),]
nmix_sim_m1_4 <- nmix_sim_m1[(nmix_sim_m1$reloc == "start")|(nmix_sim_m1$reloc == "1")|(nmix_sim_m1$reloc == "2")|(nmix_sim_m1$reloc == "3"),]
nmix_sim_k1_4 <- nmix_sim_k1[(nmix_sim_k1$reloc == "start")|(nmix_sim_k1$reloc == "1")|(nmix_sim_k1$reloc == "2")|(nmix_sim_k1$reloc == "3"),]

nmix_sim_r1_5 <- nmix_sim_r1[(nmix_sim_r1$reloc == "start")|(nmix_sim_r1$reloc == "1")|(nmix_sim_r1$reloc == "2")|(nmix_sim_r1$reloc == "3")|(nmix_sim_r1$reloc == "4"),]
nmix_sim_m1_5 <- nmix_sim_m1[(nmix_sim_m1$reloc == "start")|(nmix_sim_m1$reloc == "1")|(nmix_sim_m1$reloc == "2")|(nmix_sim_m1$reloc == "3")|(nmix_sim_m1$reloc == "4"),]
nmix_sim_k1_5 <- nmix_sim_k1[(nmix_sim_k1$reloc == "start")|(nmix_sim_k1$reloc == "1")|(nmix_sim_k1$reloc == "2")|(nmix_sim_k1$reloc == "3")|(nmix_sim_k1$reloc == "4"),]

nmix_sim_r1_10 <- nmix_sim_r1[(nmix_sim_r1$reloc != "10"),]
nmix_sim_m1_10 <- nmix_sim_m1[(nmix_sim_m1$reloc != "10"),]
nmix_sim_k1_10 <- nmix_sim_k1[(nmix_sim_k1$reloc != "10"),]

nmix_sim_r2_2 <- nmix_sim_r2[(nmix_sim_r2$reloc == "start")|(nmix_sim_r2$reloc == "1"),]
nmix_sim_m2_2 <- nmix_sim_m2[(nmix_sim_m2$reloc == "start")|(nmix_sim_m2$reloc == "1"),]
nmix_sim_k2_2 <- nmix_sim_k2[(nmix_sim_k2$reloc == "start")|(nmix_sim_k2$reloc == "1"),]

nmix_sim_r2_3 <- nmix_sim_r2[(nmix_sim_r2$reloc == "start")|(nmix_sim_r2$reloc == "1")|(nmix_sim_r2$reloc == "2"),]
nmix_sim_m2_3 <- nmix_sim_m2[(nmix_sim_m2$reloc == "start")|(nmix_sim_m2$reloc == "1")|(nmix_sim_m2$reloc == "2"),]
nmix_sim_k2_3 <- nmix_sim_k2[(nmix_sim_k2$reloc == "start")|(nmix_sim_k2$reloc == "1")|(nmix_sim_k2$reloc == "2"),]

nmix_sim_r2_4 <- nmix_sim_r2[(nmix_sim_r2$reloc == "start")|(nmix_sim_r2$reloc == "1")|(nmix_sim_r2$reloc == "2")|(nmix_sim_r2$reloc == "3"),]
nmix_sim_m2_4 <- nmix_sim_m2[(nmix_sim_m2$reloc == "start")|(nmix_sim_m2$reloc == "1")|(nmix_sim_m2$reloc == "2")|(nmix_sim_m2$reloc == "3"),]
nmix_sim_k2_4 <- nmix_sim_k2[(nmix_sim_k2$reloc == "start")|(nmix_sim_k2$reloc == "1")|(nmix_sim_k2$reloc == "2")|(nmix_sim_k2$reloc == "3"),]

nmix_sim_r2_5 <- nmix_sim_r2[(nmix_sim_r2$reloc == "start")|(nmix_sim_r2$reloc == "1")|(nmix_sim_r2$reloc == "2")|(nmix_sim_r2$reloc == "3")|(nmix_sim_r2$reloc == "4"),]
nmix_sim_m2_5 <- nmix_sim_m2[(nmix_sim_m2$reloc == "start")|(nmix_sim_m2$reloc == "1")|(nmix_sim_m2$reloc == "2")|(nmix_sim_m2$reloc == "3")|(nmix_sim_m2$reloc == "4"),]
nmix_sim_k2_5 <- nmix_sim_k2[(nmix_sim_k2$reloc == "start")|(nmix_sim_k2$reloc == "1")|(nmix_sim_k2$reloc == "2")|(nmix_sim_k2$reloc == "3")|(nmix_sim_k2$reloc == "4"),]

nmix_sim_r2_10 <- nmix_sim_r2[(nmix_sim_r2$reloc != "10"),]
nmix_sim_m2_10 <- nmix_sim_m2[(nmix_sim_m2$reloc != "10"),]
nmix_sim_k2_10 <- nmix_sim_k2[(nmix_sim_k2$reloc != "10"),]

nmix_sim_r3_2 <- nmix_sim_r3[(nmix_sim_r3$reloc == "start")|(nmix_sim_r3$reloc == "1"),]
nmix_sim_m3_2 <- nmix_sim_m3[(nmix_sim_m3$reloc == "start")|(nmix_sim_m3$reloc == "1"),]
nmix_sim_k3_2 <- nmix_sim_k3[(nmix_sim_k3$reloc == "start")|(nmix_sim_k3$reloc == "1"),]

nmix_sim_r3_3 <- nmix_sim_r3[(nmix_sim_r3$reloc == "start")|(nmix_sim_r3$reloc == "1")|(nmix_sim_r3$reloc == "2"),]
nmix_sim_m3_3 <- nmix_sim_m3[(nmix_sim_m3$reloc == "start")|(nmix_sim_m3$reloc == "1")|(nmix_sim_m3$reloc == "2"),]
nmix_sim_k3_3 <- nmix_sim_k3[(nmix_sim_k3$reloc == "start")|(nmix_sim_k3$reloc == "1")|(nmix_sim_k3$reloc == "2"),]

nmix_sim_r3_4 <- nmix_sim_r3[(nmix_sim_r3$reloc == "start")|(nmix_sim_r3$reloc == "1")|(nmix_sim_r3$reloc == "2")|(nmix_sim_r3$reloc == "3"),]
nmix_sim_m3_4 <- nmix_sim_m3[(nmix_sim_m3$reloc == "start")|(nmix_sim_m3$reloc == "1")|(nmix_sim_m3$reloc == "2")|(nmix_sim_m3$reloc == "3"),]
nmix_sim_k3_4 <- nmix_sim_k3[(nmix_sim_k3$reloc == "start")|(nmix_sim_k3$reloc == "1")|(nmix_sim_k3$reloc == "2")|(nmix_sim_k3$reloc == "3"),]

nmix_sim_r3_5 <- nmix_sim_r3[(nmix_sim_r3$reloc == "start")|(nmix_sim_r3$reloc == "1")|(nmix_sim_r3$reloc == "2")|(nmix_sim_r3$reloc == "3")|(nmix_sim_r3$reloc == "4"),]
nmix_sim_m3_5 <- nmix_sim_m3[(nmix_sim_m3$reloc == "start")|(nmix_sim_m3$reloc == "1")|(nmix_sim_m3$reloc == "2")|(nmix_sim_m3$reloc == "3")|(nmix_sim_m3$reloc == "4"),]
nmix_sim_k3_5 <- nmix_sim_k3[(nmix_sim_k3$reloc == "start")|(nmix_sim_k3$reloc == "1")|(nmix_sim_k3$reloc == "2")|(nmix_sim_k3$reloc == "3")|(nmix_sim_k3$reloc == "4"),]

nmix_sim_r3_10 <- nmix_sim_r3[(nmix_sim_r3$reloc != "10"),]
nmix_sim_m3_10 <- nmix_sim_m3[(nmix_sim_m3$reloc != "10"),]
nmix_sim_k3_10 <- nmix_sim_k3[(nmix_sim_k3$reloc != "10"),]

nmix_sim_r5_2 <- nmix_sim_r5[(nmix_sim_r5$reloc == "start")|(nmix_sim_r5$reloc == "1"),]
nmix_sim_m5_2 <- nmix_sim_m5[(nmix_sim_m5$reloc == "start")|(nmix_sim_m5$reloc == "1"),]
nmix_sim_k5_2 <- nmix_sim_k5[(nmix_sim_k5$reloc == "start")|(nmix_sim_k5$reloc == "1"),]

nmix_sim_r5_3 <- nmix_sim_r5[(nmix_sim_r5$reloc == "start")|(nmix_sim_r5$reloc == "1")|(nmix_sim_r5$reloc == "2"),]
nmix_sim_m5_3 <- nmix_sim_m5[(nmix_sim_m5$reloc == "start")|(nmix_sim_m5$reloc == "1")|(nmix_sim_m5$reloc == "2"),]
nmix_sim_k5_3 <- nmix_sim_k5[(nmix_sim_k5$reloc == "start")|(nmix_sim_k5$reloc == "1")|(nmix_sim_k5$reloc == "2"),]

nmix_sim_r5_4 <- nmix_sim_r5[(nmix_sim_r5$reloc == "start")|(nmix_sim_r5$reloc == "1")|(nmix_sim_r5$reloc == "2")|(nmix_sim_r5$reloc == "3"),]
nmix_sim_m5_4 <- nmix_sim_m5[(nmix_sim_m5$reloc == "start")|(nmix_sim_m5$reloc == "1")|(nmix_sim_m5$reloc == "2")|(nmix_sim_m5$reloc == "3"),]
nmix_sim_k5_4 <- nmix_sim_k5[(nmix_sim_k5$reloc == "start")|(nmix_sim_k5$reloc == "1")|(nmix_sim_k5$reloc == "2")|(nmix_sim_k5$reloc == "3"),]

nmix_sim_r5_5 <- nmix_sim_r5[(nmix_sim_r5$reloc == "start")|(nmix_sim_r5$reloc == "1")|(nmix_sim_r5$reloc == "2")|(nmix_sim_r5$reloc == "3")|(nmix_sim_r5$reloc == "4"),]
nmix_sim_m5_5 <- nmix_sim_m5[(nmix_sim_m5$reloc == "start")|(nmix_sim_m5$reloc == "1")|(nmix_sim_m5$reloc == "2")|(nmix_sim_m5$reloc == "3")|(nmix_sim_m5$reloc == "4"),]
nmix_sim_k5_5 <- nmix_sim_k5[(nmix_sim_k5$reloc == "start")|(nmix_sim_k5$reloc == "1")|(nmix_sim_k5$reloc == "2")|(nmix_sim_k5$reloc == "3")|(nmix_sim_k5$reloc == "4"),]

nmix_sim_r5_10 <- nmix_sim_r5[(nmix_sim_r5$reloc != "10"),]
nmix_sim_m5_10 <- nmix_sim_m5[(nmix_sim_m5$reloc != "10"),]
nmix_sim_k5_10 <- nmix_sim_k5[(nmix_sim_k5$reloc != "10"),]

nmix_sim_r10_2 <- nmix_sim_r10[(nmix_sim_r10$reloc == "start")|(nmix_sim_r10$reloc == "1"),]
nmix_sim_m10_2 <- nmix_sim_m10[(nmix_sim_m10$reloc == "start")|(nmix_sim_m10$reloc == "1"),]
nmix_sim_k10_2 <- nmix_sim_k10[(nmix_sim_k10$reloc == "start")|(nmix_sim_k10$reloc == "1"),]

nmix_sim_r10_3 <- nmix_sim_r10[(nmix_sim_r10$reloc == "start")|(nmix_sim_r10$reloc == "1")|(nmix_sim_r10$reloc == "2"),]
nmix_sim_m10_3 <- nmix_sim_m10[(nmix_sim_m10$reloc == "start")|(nmix_sim_m10$reloc == "1")|(nmix_sim_m10$reloc == "2"),]
nmix_sim_k10_3 <- nmix_sim_k10[(nmix_sim_k10$reloc == "start")|(nmix_sim_k10$reloc == "1")|(nmix_sim_k10$reloc == "2"),]

nmix_sim_r10_4 <- nmix_sim_r10[(nmix_sim_r10$reloc == "start")|(nmix_sim_r10$reloc == "1")|(nmix_sim_r10$reloc == "2")|(nmix_sim_r10$reloc == "3"),]
nmix_sim_m10_4 <- nmix_sim_m10[(nmix_sim_m10$reloc == "start")|(nmix_sim_m10$reloc == "1")|(nmix_sim_m10$reloc == "2")|(nmix_sim_m10$reloc == "3"),]
nmix_sim_k10_4 <- nmix_sim_k10[(nmix_sim_k10$reloc == "start")|(nmix_sim_k10$reloc == "1")|(nmix_sim_k10$reloc == "2")|(nmix_sim_k10$reloc == "3"),]

nmix_sim_r10_5 <- nmix_sim_r10[(nmix_sim_r10$reloc == "start")|(nmix_sim_r10$reloc == "1")|(nmix_sim_r10$reloc == "2")|(nmix_sim_r10$reloc == "3")|(nmix_sim_r10$reloc == "4"),]
nmix_sim_m10_5 <- nmix_sim_m10[(nmix_sim_m10$reloc == "start")|(nmix_sim_m10$reloc == "1")|(nmix_sim_m10$reloc == "2")|(nmix_sim_m10$reloc == "3")|(nmix_sim_m10$reloc == "4"),]
nmix_sim_k10_5 <- nmix_sim_k10[(nmix_sim_k10$reloc == "start")|(nmix_sim_k10$reloc == "1")|(nmix_sim_k10$reloc == "2")|(nmix_sim_k10$reloc == "3")|(nmix_sim_k10$reloc == "4"),]

nmix_sim_r10_10 <- nmix_sim_r10[(nmix_sim_r10$reloc != "10"),]
nmix_sim_m10_10 <- nmix_sim_m10[(nmix_sim_m10$reloc != "10"),]
nmix_sim_k10_10 <- nmix_sim_k10[(nmix_sim_k10$reloc != "10"),]

candidate.nmix <- list()
candidate.nmix[["known.1.2"]] <- nmix_sim_k1_2
candidate.nmix[["known.1.3"]] <- nmix_sim_k1_3
candidate.nmix[["known.1.4"]] <- nmix_sim_k1_4
candidate.nmix[["known.1.5"]] <- nmix_sim_k1_5
candidate.nmix[["known.1.10"]] <- nmix_sim_k1_10
candidate.nmix[["known.2.2"]] <- nmix_sim_k2_2
candidate.nmix[["known.2.3"]] <- nmix_sim_k2_3
candidate.nmix[["known.2.4"]] <- nmix_sim_k2_4
candidate.nmix[["known.2.5"]] <- nmix_sim_k2_5
candidate.nmix[["known.2.10"]] <- nmix_sim_k2_10
candidate.nmix[["known.3.2"]] <- nmix_sim_k3_2
candidate.nmix[["known.3.3"]] <- nmix_sim_k3_3
candidate.nmix[["known.3.4"]] <- nmix_sim_k3_4
candidate.nmix[["known.3.5"]] <- nmix_sim_k3_5
candidate.nmix[["known.3.10"]] <- nmix_sim_k3_10
candidate.nmix[["known.5.2"]] <- nmix_sim_k5_2
candidate.nmix[["known.5.3"]] <- nmix_sim_k5_3
candidate.nmix[["known.5.4"]] <- nmix_sim_k5_4
candidate.nmix[["known.5.5"]] <- nmix_sim_k5_5
candidate.nmix[["known.5.10"]] <- nmix_sim_k5_10
candidate.nmix[["known.10.2"]] <- nmix_sim_k10_2
candidate.nmix[["known.10.3"]] <- nmix_sim_k10_3
candidate.nmix[["known.10.4"]] <- nmix_sim_k10_4
candidate.nmix[["known.10.5"]] <- nmix_sim_k10_5
candidate.nmix[["known.10.10"]] <- nmix_sim_k10_10

candidate.nmix[["manual.1.2"]] <- nmix_sim_m1_2
candidate.nmix[["manual.1.3"]] <- nmix_sim_m1_3
candidate.nmix[["manual.1.4"]] <- nmix_sim_m1_4
candidate.nmix[["manual.1.5"]] <- nmix_sim_m1_5
candidate.nmix[["manual.1.10"]] <- nmix_sim_m1_10
candidate.nmix[["manual.2.2"]] <- nmix_sim_m2_2
candidate.nmix[["manual.2.3"]] <- nmix_sim_m2_3
candidate.nmix[["manual.2.4"]] <- nmix_sim_m2_4
candidate.nmix[["manual.2.5"]] <- nmix_sim_m2_5
candidate.nmix[["manual.2.10"]] <- nmix_sim_m2_10
candidate.nmix[["manual.3.2"]] <- nmix_sim_m3_2
candidate.nmix[["manual.3.3"]] <- nmix_sim_m3_3
candidate.nmix[["manual.3.4"]] <- nmix_sim_m3_4
candidate.nmix[["manual.3.5"]] <- nmix_sim_m3_5
candidate.nmix[["manual.3.10"]] <- nmix_sim_m3_10
candidate.nmix[["manual.5.2"]] <- nmix_sim_m5_2
candidate.nmix[["manual.5.3"]] <- nmix_sim_m5_3
candidate.nmix[["manual.5.4"]] <- nmix_sim_m5_4
candidate.nmix[["manual.5.5"]] <- nmix_sim_m5_5
candidate.nmix[["manual.5.10"]] <- nmix_sim_m5_10
candidate.nmix[["manual.10.2"]] <- nmix_sim_m10_2
candidate.nmix[["manual.10.3"]] <- nmix_sim_m10_3
candidate.nmix[["manual.10.4"]] <- nmix_sim_m10_4
candidate.nmix[["manual.10.5"]] <- nmix_sim_m10_5
candidate.nmix[["manual.10.10"]] <- nmix_sim_m10_10

candidate.nmix[["recognizer.1.2"]] <- nmix_sim_r1_2
candidate.nmix[["recognizer.1.3"]] <- nmix_sim_r1_3
candidate.nmix[["recognizer.1.4"]] <- nmix_sim_r1_4
candidate.nmix[["recognizer.1.5"]] <- nmix_sim_r1_5
candidate.nmix[["recognizer.1.10"]] <- nmix_sim_r1_10
candidate.nmix[["recognizer.2.2"]] <- nmix_sim_r2_2
candidate.nmix[["recognizer.2.3"]] <- nmix_sim_r2_3
candidate.nmix[["recognizer.2.4"]] <- nmix_sim_r2_4
candidate.nmix[["recognizer.2.5"]] <- nmix_sim_r2_5
candidate.nmix[["recognizer.2.10"]] <- nmix_sim_r2_10
candidate.nmix[["recognizer.3.2"]] <- nmix_sim_r3_2
candidate.nmix[["recognizer.3.3"]] <- nmix_sim_r3_3
candidate.nmix[["recognizer.3.4"]] <- nmix_sim_r3_4
candidate.nmix[["recognizer.3.5"]] <- nmix_sim_r3_5
candidate.nmix[["recognizer.3.10"]] <- nmix_sim_r3_10
candidate.nmix[["recognizer.5.2"]] <- nmix_sim_r5_2
candidate.nmix[["recognizer.5.3"]] <- nmix_sim_r5_3
candidate.nmix[["recognizer.5.4"]] <- nmix_sim_r5_4
candidate.nmix[["recognizer.5.5"]] <- nmix_sim_r5_5
candidate.nmix[["recognizer.5.10"]] <- nmix_sim_r5_10
candidate.nmix[["recognizer.10.2"]] <- nmix_sim_r10_2
candidate.nmix[["recognizer.10.3"]] <- nmix_sim_r10_3
candidate.nmix[["recognizer.10.4"]] <- nmix_sim_r10_4
candidate.nmix[["recognizer.10.5"]] <- nmix_sim_r10_5
candidate.nmix[["recognizer.10.10"]] <- nmix_sim_r10_10

esa = (pi * ((max(m$ddf$data$distance) * sqrt(m$ddf$fitted[1]))^2))/10000

results.nmix <- data.frame()
f=1
for(i in candidate.nmix){
  nmix.df <- i %>% count(survey.no, reloc)
  nmix.df <- dcast(nmix.df, survey.no ~ reloc)
  nmix.df[is.na(nmix.df)] <- 0
  colnames(nmix.df)[(names(nmix.df) == "start")] <- "0"
  colnames(nmix.df) <- paste("V", colnames(nmix.df), sep = "")
  colnames(nmix.df)[(names(nmix.df) == "Vsurvey.no")] <- "survey.no"
  nmix.df <- within(nmix.df, rm("survey.no"))
  
  nmix.UMF <- unmarkedFramePCount(nmix.df, siteCovs = NULL, obsCov = NULL)
  nmix.fm <- pcount(formula = ~1 ~ 1, data = nmix.UMF, K = 50)
  
  nmix.dens <- nmix.fm@estimates@estimates[["state"]]@estimates / esa
  se <- SE(nmix.fm)[1]
  listname <- names(candidate.nmix[f])
  slen <- strsplit(listname, "[.]")[[1]][2]
  svisit <- strsplit(listname, "[.]")[[1]][3]
  ms <- strsplit(listname, "[.]")[[1]][1]
  
  results.nmix.temp <- data.frame(method = "n-mixture", measure = ms, visits = svisit, length = slen, 
                                  estimate = nmix.dens, se = se)
  results.nmix <- rbind(results.nmix, results.nmix.temp)
  print(f)
  f=f+1
  
}
rownames(results.nmix) <- c()


#********************************#
#********************************#
####SINGLE VISIT APPROACH####
#********************************#
#********************************#
recognizer_singlept <- singlept_data
recognizer_dist <- singlept_data$recognizer
recognizer_singlept$probability <- p <- predictDistFunction(recognizer_dist, xy)
recognizer_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

manual_singlept <- singlept_data
manual_dist <- singlept_data$manual
manual_singlept$probability <- p <- predictDistFunction(manual_dist, xy)
manual_singlept$detected <- rbinom(length(p), size=1, prob=pmax(0.0001, pmin(0.9999, p)))

detections_singlept_r <- recognizer_singlept[recognizer_singlept$detected == 1, ]
detections_singlept_m <- manual_singlept[manual_singlept$detected == 1, ]

sv_sim_r <- detections_singlept_r
sv_sim_m <- detections_singlept_m
sv_sim_k <- detections_singlept

sv_sim_r$distance <- sv_sim_r$recognizer
sv_sim_m$distance <- sv_sim_m$manual

sv_sim_r <- sv_sim_r[(sv_sim_r$reloc == "start"),]
sv_sim_m <- sv_sim_m[(sv_sim_m$reloc == "start"),]
sv_sim_k <- sv_sim_k[(sv_sim_k$reloc == "start"),]

sv_sim_r1 <- sv_sim_r[unique(row(sv_sim_r[-c(1:5, 7:21)])[sv_sim_r[-c(1:5, 7:21)] == "1"]),]
sv_sim_m1 <- sv_sim_m[unique(row(sv_sim_m[-c(1:5, 7:21)])[sv_sim_m[-c(1:5, 7:21)] == "1"]),]
sv_sim_k1 <- sv_sim_k[unique(row(sv_sim_k[-c(1:5, 7:21)])[sv_sim_k[-c(1:5, 7:21)] == "1"]),]

sv_sim_r2 <- sv_sim_r[unique(row(sv_sim_r[-c(1:5, 8:21)])[sv_sim_r[-c(1:5, 8:21)] == "1"]),]
sv_sim_m2 <- sv_sim_m[unique(row(sv_sim_m[-c(1:5, 8:21)])[sv_sim_m[-c(1:5, 8:21)] == "1"]),]
sv_sim_k2 <- sv_sim_k[unique(row(sv_sim_k[-c(1:5, 8:21)])[sv_sim_k[-c(1:5, 8:21)] == "1"]),]

sv_sim_r3 <- sv_sim_r[unique(row(sv_sim_r[-c(1:5, 9:21)])[sv_sim_r[-c(1:5, 9:21)] == "1"]),]
sv_sim_m3 <- sv_sim_m[unique(row(sv_sim_m[-c(1:5, 9:21)])[sv_sim_m[-c(1:5, 9:21)] == "1"]),]
sv_sim_k3 <- sv_sim_k[unique(row(sv_sim_k[-c(1:5, 9:21)])[sv_sim_k[-c(1:5, 9:21)] == "1"]),]

sv_sim_r5 <- sv_sim_r[unique(row(sv_sim_r[-c(1:5, 11:21)])[sv_sim_r[-c(1:5, 11:21)] == "1"]),]
sv_sim_m5 <- sv_sim_m[unique(row(sv_sim_m[-c(1:5, 11:21)])[sv_sim_m[-c(1:5, 11:21)] == "1"]),]
sv_sim_k5 <- sv_sim_k[unique(row(sv_sim_k[-c(1:5, 11:21)])[sv_sim_k[-c(1:5, 11:21)] == "1"]),]

sv_sim_r10 <- sv_sim_r
sv_sim_m10 <- sv_sim_m
sv_sim_k10 <- sv_sim_k

candidate.sv <- list()
candidate.sv[["known.1"]] <- sv_sim_k1
candidate.sv[["known.2"]] <- sv_sim_k2
candidate.sv[["known.3"]] <- sv_sim_k3
candidate.sv[["known.5"]] <- sv_sim_k5
candidate.sv[["known.10"]] <- sv_sim_k10
candidate.sv[["manual.1"]] <- sv_sim_m1
candidate.sv[["manual.2"]] <- sv_sim_m2
candidate.sv[["manual.3"]] <- sv_sim_m3
candidate.sv[["manual.5"]] <- sv_sim_m5
candidate.sv[["manual.10"]] <- sv_sim_m10
candidate.sv[["recognizer.1"]] <- sv_sim_r1
candidate.sv[["recognizer.2"]] <- sv_sim_r2
candidate.sv[["recognizer.3"]] <- sv_sim_r3
candidate.sv[["recognizer.5"]] <- sv_sim_r5
candidate.sv[["recognizer.10"]] <- sv_sim_r10

esa = (pi * ((max(m$ddf$data$distance) * sqrt(m$ddf$fitted[1]))^2))/10000

stderr <- function(x) sqrt(var(x)/length(x))
results.sv <- data.frame()
f=1
for(i in candidate.sv){
  sv.df <- i %>% count(survey.no, reloc)
  sv.df <- dcast(sv.df, survey.no ~ reloc)
  sv.df[is.na(sv.df)] <- 0
  
  est <- ((sum(sv.df$start))/s_iter) / esa
  se <- stderr(sv.df$start)
  listname <- names(candidate.sv[f])
  slen <- strsplit(listname, "[.]")[[1]][2]
  ms <- strsplit(listname, "[.]")[[1]][1]
  
  results.sv.temp <- data.frame(method = "distance", measure = ms, length = slen, estimate = est,
                                se = se)
  results.sv <- rbind(results.sv, results.sv.temp)
  f=f+1
}
rownames(results.sv) <- c()



#********************************#
#********************************#
####SPATIALLY-EXPLICIT CAPTURE-RECAPTURE####
#********************************#
#********************************#


















