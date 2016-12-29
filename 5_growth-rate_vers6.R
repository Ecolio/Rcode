# last script for extracting the data out the compilation and plotting the results


rm(list = ls())


# define the xy position of the lane 
#endof1treat <- 15
#endof2treat <- 25
#endof3treat <- 40



# needed libraries
library(data.table)
library(utils)
library(reshape)
library(ggplot2)

###############################################

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=TRUE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop = .drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm = T),
                     mean = mean   (xx[[col]], na.rm = T),
                     sd   = sd     (xx[[col]], na.rm = T)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

###############################################

# will set the file location and worfoldering directory (library(utils) needed):
# to install, http://cran.r-project.org/web/pacfolderages/R.utils/index.html, 
# download: r-devel: R.utils_1.34.0.zip, install from zip (not cran)
compdirname <- choose.dir(getwd(), "Choose a suitable folder")
filenames <- list.files(compdirname)
setwd(compdirname)
# read the data in
rawdata1 <- read.csv(filenames[2], header = T)
rawdata2 <- read.csv(filenames[3], header = T)
rawdata3 <- read.csv(filenames[4], header = T)
rawdata4 <- read.csv(filenames[5], header = T)
rawdata5 <- read.csv(filenames[6], header = T)
rawdata6 <- read.csv(filenames[7], header = T)
rawdata7 <- read.csv(filenames[8], header = T)
rawdata8 <- read.csv(filenames[9], header = T)
rawdata9 <- read.csv(filenames[10], header = T)




# name of the treatments
nb1 <- nrow(rawdata1)
nb2 <- nrow(rawdata2)
nb3 <- nrow(rawdata3)
nb4 <- nrow(rawdata4)
nb5 <- nrow(rawdata5)
nb6 <- nrow(rawdata6)
nb7 <- nrow(rawdata7)
nb8 <- nrow(rawdata8)
nb9 <- nrow(rawdata9)
nbtot <- nb1 + nb2 + nb3 + nb4 + nb5 + nb6 + nb7 + nb8 + nb9
nameof1treat <- paste("33°C, N = ", nb1, "cells")
nameof2treat <- paste("35°C, N = ", nb2, "cells")
nameof3treat <- paste("36°C, N = ", nb3, "cells")
nameof4treat <- paste("37°C, N = ", nb4, "cells")
nameof5treat <-  paste("38°C, N = ", nb5, "cells")
nameof6treat <- paste("39°C, N = ", nb6, "cells")
nameof7treat <- paste("40°C, N = ", nb7, "cells")
nameof8treat <-  paste("41°C, N = ", nb8, "cells")
nameof9treat <-  paste("42°C, N = ", nb9, "cells")

##############################################################################################################
 ################################################################################
 # mean by time of the growth rates
# building the long table ####################################################################################

## 33
 #variables
 Nbrow <- nrow(rawdata1)
 longdata1 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata1[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata1[k, 4:750])
   temp$ToD <- rawdata1[k, 3]
   temp$xy <- rawdata1[k, 2]
   temp$Treat <- nameof1treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata1 <- rbind(longdata1,temp)
 }

## 35
 Nbrow <- nrow(rawdata2)
 longdata2 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata2[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata2[k, 4:750])
   temp$ToD <- rawdata2[k, 3]
   temp$xy <- rawdata2[k, 2]
   temp$Treat <- nameof2treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata2 <- rbind(longdata2,temp)
 }

## 36
 Nbrow <- nrow(rawdata3)
 longdata3 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata3[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata3[k, 4:750])
   temp$ToD <- rawdata3[k, 3]
   temp$xy <- rawdata3[k, 2]
   temp$Treat <- nameof3treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata3 <- rbind(longdata3,temp)
 }

## 37
 Nbrow <- nrow(rawdata4)
 longdata4 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata4[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata4[k, 4:750])
   temp$ToD <- rawdata4[k, 3]
   temp$xy <- rawdata4[k, 2]
   temp$Treat <- nameof4treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata4 <- rbind(longdata4,temp)
 }

## 38
 Nbrow <- nrow(rawdata5)
 longdata5 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata5[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata5[k, 4:750])
   temp$ToD <- rawdata5[k, 3]
   temp$xy <- rawdata5[k, 2]
   temp$Treat <- nameof5treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata5 <- rbind(longdata5,temp)
 }

## 39
 Nbrow <- nrow(rawdata6)
 longdata6 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata6[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata6[k, 4:750])
   temp$ToD <- rawdata6[k, 3]
   temp$xy <- rawdata6[k, 2]
   temp$Treat <- nameof6treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata6 <- rbind(longdata6,temp)
 }

## 40
 Nbrow <- nrow(rawdata7)
 longdata7 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata7[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata7[k, 4:750])
   temp$ToD <- rawdata7[k, 3]
   temp$xy <- rawdata7[k, 2]
   temp$Treat <- nameof7treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata7 <- rbind(longdata7,temp)
 }

## 41
 Nbrow <- nrow(rawdata8)
 longdata8 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata8[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata8[k, 4:750])
   temp$ToD <- rawdata8[k, 3]
   temp$xy <- rawdata8[k, 2]
   temp$Treat <- nameof8treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata8 <- rbind(longdata8,temp)
 }

## 42
 Nbrow <- nrow(rawdata9)
 longdata9 <- as.data.frame(c())
 k <- 1
  #data retrieval
  for(k in ( 1 : Nbrow)) {
   dummylength <- ncol(rawdata9[k, 4:750])
   temp <- as.data.frame((matrix(nrow = dummylength, ncol = 1)))
   temp$Mean <-  t(rawdata9[k, 4:750])
   temp$ToD <- rawdata9[k, 3]
   temp$xy <- rawdata9[k, 2]
   temp$Treat <- nameof9treat
   temp <- subset(temp, !is.na(temp$Mean))
   dummylength <- nrow(temp)
   temp$V1 <- c(1: dummylength) * 3
   longdata9 <- rbind(longdata9,temp)
 }

   longdata <- rbind(longdata1, longdata2, longdata3, 
                     longdata4, longdata5, longdata6,
                     longdata7, longdata8, longdata9)

   write.csv(longdata, file = ".\\longdata\\longdata-gr.csv")
   longdata <- read.csv(filenames[3])
   

 datameanc <- summarySE(longdata, 
                        measurevar="Mean", 
                        groupvars= c("V1", "Treat"), 
                        na.rm = TRUE)
 
 datameanc <- subset(datameanc, !is.na(datameanc$sd))
 datameanc$Treat <- factor(datameanc$Treat, 
                           levels = c(nameof1treat, nameof2treat, nameof3treat, 
                                      nameof4treat, nameof5treat, nameof6treat, 
                                      nameof7treat, nameof8treat, nameof9treat))
 datameanc$V1 <- datameanc$V1 * 0.2
 
 
 ###########################################################################
 limits <- aes(ymax = Mean + ci, ymin=Mean - ci)
 ggplot(datameanc[10:3691,], aes(V1,Mean, group = Treat)) + 
   geom_point(shape = 21, colour = "chartreuse3", fill = "chartreuse3", size = 0.1, stroke = 0.1)+
   facet_wrap(~Treat, ncol=3, scales="fixed")+
   geom_errorbar(limits, width=0.1) +
   xlim(0,200) + 
   ylim(0.91, 1.11) +
   labs(x="Time (hour)", y = "Mean growth rate") + 
   ggtitle(expression(atop("Growth rate", atop(italic("in a temperature range"), "")))) +
   geom_ribbon(aes(ymin=Mean-ci, ymax=Mean+ci), alpha=0.5, fill="chartreuse2") +
   stat_smooth(method = "lm", col = "red")
 ############################################################################################################################################################

 