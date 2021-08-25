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
require("data.table")
paluse <- brewer.spectral(9)
stateid <- "Argentina"
setwd('')
arg <- read.csv("Argentina_manual.csv") #from https://ais.paho.org/phip/viz/ed_flu.asp

# cleaning
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

##########################################################
# run model fitting beta1 and phi

controlLength <- 14 # at this stage, control length does not matter 
timestart = 1914

beta1vec <- seq(0.1,0.3,0.01)
phivec <- seq(1,6.5,0.5)
matout <- matrix(NA, nrow=length(beta1vec), ncol=length(phivec))
for(i in 1:length(beta1vec)){
  for(j in 1:length(phivec)){
    r <- run_rsv_model(b0 = 0.03,
                   b1 = beta1vec[i],
                   phi =phivec[j],
                   prop_detected_1 = 0.424, # hospitalization rate
                   prop_detected_2 = 0.088,
                   prop_detected_3 = 0.047,
                   prop_detected_4 = 0.020,
                   max_t = max_t,
                   mixing = mixing,
                   timestart = timestart,
                   timeend = timestart + controlLength, 
                   betared = 0.6,
                   init_conds_from_file = 0, maternalimmunity = TRUE, popInput = 44940000)

#all_inc <- r$Incidence
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
mod_dets <- mod_dets*(mean(obs_dets)/mean(mod_dets))
mean(mod_dets)
mean(obs_dets)
matout[i,j] <- sum((obs_dets - mod_dets)^2)
  }
}

save(matout, file = paste0(stateid,".RData"))
write.csv(matout, file = paste0(stateid,".csv"))

