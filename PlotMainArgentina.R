library(odin)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(coda)
library(pracma)
require("pals")
require("lubridate")
require("fields")
paluse <- brewer.spectral(9)
stateid <- "Argentina"
popset <- 44940000
setwd(" ")
arg <- read.csv("Argentina_manual.csv") # from https://ais.paho.org/phip/viz/ed_flu.asp

#cleaning
arg$pctpos <- as.numeric(arg$Moving.Average.of.2_P_VSR)
arg$cases <- (arg$pctpos)*1000
arg$year <- arg$Year
arg$week <- arg$EW
arg$date <- as.Date(paste0(arg$year, "-",arg$week,"-", 1), "%Y-%U-%u")
arg$month <- month(arg$date) 
arg <- arg %>% group_by(month, year) %>% summarize(sum = sum(cases))
arg <- arg[order(arg$year, arg$month),]
argst <- arg
arg <- arg[arg$year >= 2016 & arg$year < 2020,]
argav <- arg %>% group_by(month) %>% summarize(mean_cases = mean(sum))

setwd(" ")
source("functions.R")
##########################################################
# load model object
x <- odin::odin("odin_model.R")

##########################################################
init_conds_from_file <- 0 # choose whether to read in some existing ICs
save_init_conds <- 1 # choose whether to save final model state as ICs for next time
max_t <- 2000

# mixing matrix - this is an adapted matrix from Mossong et al
mixing <- as.matrix(read.csv("data/mixing_75.csv", header = TRUE), header = TRUE)*365/12


### RUN WITH LOADED DATA
beta1vec <- seq(0.1,0.3,0.01)
phivec <- seq(1,6.5,0.5)

sta <- read.csv(paste0(stateid,".csv"))
sta <- sta[,-1]

st <- which(sta == min(sta), arr.ind = TRUE)

controlLength <-12 # months
timestart = 1912
r <- run_rsv_model(b0 = 0.03,
                   b1 = beta1vec[st[[1]]],
                   phi =phivec[st[[2]]],
                   prop_detected_1 = 0.424, # hospitalization rate
                   prop_detected_2 = 0.088,
                   prop_detected_3 = 0.047,
                   prop_detected_4 = 0.020,
                   max_t = max_t,
                   mixing = mixing,
                   timestart = timestart,
                   timeend = timestart + controlLength, 
                   betared = 0.53,
                   init_conds_from_file = 0, maternalimmunity = TRUE, popInput = popset)



all_inc <- r$Incidence
all_inc <- r$DetIncidence
all_inc <- all_inc[1880:2000,]
all_inc_sum <- as.data.frame(all_inc)%>%
  summarize(mon1_3 = V1 + V2 + V3,
            mon4_6 = V4 + V5 + V6,
            mon7_9 = V7 + V8 + V9, 
            mon10_12= V10 + V11 + V12,
            mon13_15 = V13 + V14 + V15, 
            mon16_18 = V16 + V17 + V18,
            mon19_21 = V19 + V20 + V21,
            mon22_24 = V22 + V23 + V24, 
            mon25_60 = V25 + V26 + V27 + V28 + V29 + V30 + V31 + V32 + V33 + V34 + V35 + V36 + V37 + V38 + V39 + V40 + 
              V41 + V42 + V43 + V44 + V45 + V46 + V47 + V48 + V49 + V50 + 
        
                    V51 + V52 + V53 + V54 + V55 + V56 + V57 + V58 + V59 + V60)

obs_dets  <-  argav$mean_cases
mod_dets <- rowSums(all_inc_sum)[1:length(obs_dets)]
mod_rr <- (mean(obs_dets)/mean(mod_dets))
mod_dets <- mod_dets*mod_rr
plot(mod_dets)
lines(obs_dets)


require("plotrix")

pdf(paste0("agestructure_",stateid,"DetI.pdf"), width=10, height=3)
layout(matrix(c(1,1,1,2,2), nrow = 1, ncol = 5, byrow = TRUE))
par(mar=c(3,3,1,1))
par(cex = 0.8)
barplot(t(all_inc_sum)*mod_rr,border = NA, col = paluse, space = 0, xlab="")
axislabs <- seq(2017,2027,1)[1:11]
axis(1, at = seq(1,121,12),label = axislabs )
timeuse <- seq(1,121,1)
polygon(c(timeuse[1914 - 1874],timeuse[1914 - 1874],
          timeuse[timestart + controlLength - 1874],timeuse[timestart + controlLength - 1874]), 
        c(-1,sum(all_inc_sum),sum(all_inc_sum),-1), border = NA,col=rgb(0,0,0,0.3))
obs_dets  <-  argav$mean_cases
argst <- argst[argst$year >= 2017,]
mod_dets <- rowSums(all_inc_sum)[1:length(obs_dets)]
points(seq(0.5,2000.5,1)[1:length(argst$sum)],argst$sum, pch= 16, cex = 0.8, col="grey26")
lines(seq(0.5,2000.5,1)[1:length(argst$sum)],argst$sum, pch= 16, cex = 0.8, col="grey64")
title(xlab="Year", line =2)
title(ylab="Hospitalizations", line = 2)


timeyears <- seq(2017,2027,1/12)[1:length(timeuse)]
all_inc_postcontrol <- all_inc_sum[timeyears >= 2021.1 & timeyears < 2022.1,]
postcontroldis <-  colSums(all_inc_postcontrol)*mod_rr
postcontroldis_norm <- postcontroldis/sum(postcontroldis)
all_inc_nocontrol <- all_inc_sum[floor(timeyears)==2018,]
nocontroldis <- colSums(all_inc_nocontrol)*mod_rr
nocontroldis_norm <- nocontroldis/sum(nocontroldis)

bartable <- as.matrix(rbind(round(nocontroldis[1:8]),round(postcontroldis[1:8])))
colnames(bartable) <- c("0-2","3-5","6-8","9-11","12-14","15-17","18-20","21-23")
barplot(bartable, beside = TRUE, border = "NA", col= c(paluse[1], paluse[9]),xlab="", cex.names = 0.8) 
title(xlab="Age (months)", line = 2)
title(ylab="Total hospitalizations", line = 2)
legend("topright",
       c("Typical","Post-NPI"),
       fill = c(paluse[1], paluse[9]))
dev.off()