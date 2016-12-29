# last script for extracting the data out the compilation and plotting the results


rm(list = ls())



# needed libraries
library(data.table)
library(utils)
library(reshape)
library(ggplot2)
library(zoo)


#will set the file location and worfoldering directory (library(utils) needed):
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

#estimates lx, qx, for initial cells
T1dummy <-rawdata1$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT1<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT1<-data.frame(survrateT1)
colnames(survrateT1)<-c("dx", "surv", "lx","age","qx","sx")
survrateT1<-survrateT1[1:max(T1dummy)-1,]
survrateT1$Treat<-nameof1treat

#estimates lx, qx, for initial cells
T2dummy <-rawdata2$V2
T2dummy <-subset(T2dummy, !is.na(T2dummy))
T2dummy <- sort(T2dummy)
length(T2dummy)
T2dummy2<-rep(NA,(ceiling(max(T2dummy)/3)))
for (i in 0:(ceiling(max(T2dummy)/3))) {
  T2dummy2[ i + 1 ] <- length((which(T2dummy >= i * 3 + 1 & T2dummy <= i * 3 + 3)))
}
dx<-T2dummy2
sum(dx)
T2dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T2dummy2)+1))
lx[2:length(lx)]<-T2dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T2dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT2<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT2<-data.frame(survrateT2)
colnames(survrateT2)<-c("dx", "surv", "lx","age","qx","sx")
survrateT2<-survrateT2[1:max(T2dummy)-1,]
survrateT2$Treat<-nameof2treat

#estimates lx, qx, for initial cells
T3dummy <-rawdata3$V2
T3dummy <-subset(T3dummy, !is.na(T3dummy))
T3dummy <- sort(T3dummy)
length(T3dummy)
T3dummy2<-rep(NA,(ceiling(max(T3dummy)/3)))
for (i in 0:(ceiling(max(T3dummy)/3))) {
  T3dummy2[ i + 1 ] <- length((which(T3dummy >= i * 3 + 1 & T3dummy <= i * 3 + 3)))
}
dx<-T3dummy2
sum(dx)
T3dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T3dummy2)+1))
lx[2:length(lx)]<-T3dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T3dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT3<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT3<-data.frame(survrateT3)
colnames(survrateT3)<-c("dx", "surv", "lx","age","qx","sx")
survrateT3<-survrateT3[1:max(T3dummy)-1,]
survrateT3$Treat<-nameof3treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata4$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT4<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT4<-data.frame(survrateT4)
colnames(survrateT4)<-c("dx", "surv", "lx","age","qx","sx")
survrateT4<-survrateT4[1:max(T1dummy)-1,]
survrateT4$Treat<-nameof4treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata5$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT5<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT5<-data.frame(survrateT5)
colnames(survrateT5)<-c("dx", "surv", "lx","age","qx","sx")
survrateT5<-survrateT5[1:max(T1dummy)-1,]
survrateT5$Treat<-nameof5treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata6$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT6<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT6<-data.frame(survrateT6)
colnames(survrateT6)<-c("dx", "surv", "lx","age","qx","sx")
survrateT6<-survrateT6[1:max(T1dummy)-1,]
survrateT6$Treat<-nameof6treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata7$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT7<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT7<-data.frame(survrateT7)
colnames(survrateT7)<-c("dx", "surv", "lx","age","qx","sx")
survrateT7<-survrateT7[1:max(T1dummy)-1,]
survrateT7$Treat<-nameof7treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata8$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT8<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT8<-data.frame(survrateT8)
colnames(survrateT8)<-c("dx", "surv", "lx","age","qx","sx")
survrateT8<-survrateT8[1:max(T1dummy)-1,]
survrateT8$Treat<-nameof8treat

#estimates lx, qx, for initial cells
T1dummy <-rawdata9$V2
T1dummy <-subset(T1dummy, !is.na(T1dummy))
T1dummy <- sort(T1dummy)
length(T1dummy)
T1dummy2<-rep(NA,(ceiling(max(T1dummy)/3)))
for (i in 0:(ceiling(max(T1dummy)/3))) {
  T1dummy2[ i + 1 ] <- length((which(T1dummy >= i * 3 + 1 & T1dummy <= i * 3 + 3)))
}
dx<-T1dummy2
sum(dx)
T1dummy3<-sum(dx)-cumsum(dx)
lx<-rep(NA,(length(T1dummy2)+1))
lx[2:length(lx)]<-T1dummy3
lx[1]<-sum(dx)
age<-seq(1, max(T1dummy), 3)
qx<-dx/(lx[1:length(lx)-1])
sx<-cumprod(1-qx)
survrateT9<-cbind(dx,lx[2:length(lx)],lx[1:length(lx)-1],age,qx,sx)
survrateT9<-data.frame(survrateT9)
colnames(survrateT9)<-c("dx", "surv", "lx","age","qx","sx")
survrateT9<-survrateT9[1:max(T1dummy)-1,]
survrateT9$Treat<-nameof9treat


survT1 <- subset(survrateT1, !is.nan(survrateT1[, 5]) & !is.na(survrateT1[, 1]))
nbsurvT1 <- nrow(survT1)
survT1 <- survT1[1 : nbsurvT1 - 1, ]
survT1$m.av <- rollmean(survT1$qx, 15,fill = list(NA, NULL, NA))

survT2 <- subset(survrateT2, !is.nan(survrateT2[, 5]) & !is.na(survrateT2[, 1]))
nbsurvT2 <- nrow(survT2)
survT2 <- survT2[1 : nbsurvT2 - 1, ]
survT2$m.av <- rollmean(survT2$qx, 15,fill = list(NA, NULL, NA))

survT3 <- subset(survrateT3, !is.nan(survrateT3[, 5]) & !is.na(survrateT3[, 1]))
nbsurvT3 <- nrow(survT3)
survT3 <- survT3[1 : nbsurvT3 - 1, ]
survT3$m.av <- rollmean(survT3$qx, 15,fill = list(NA, NULL, NA))

survT4 <- subset(survrateT4, !is.nan(survrateT4[, 5]) & !is.na(survrateT4[, 1]))
nbsurvT4 <- nrow(survT4)
survT4 <- survT4[1 : nbsurvT4 - 1, ]
survT4$m.av <- rollmean(survT4$qx, 15,fill = list(NA, NULL, NA))

survT5 <- subset(survrateT5, !is.nan(survrateT5[, 5]) & !is.na(survrateT5[, 1]))
nbsurvT5 <- nrow(survT5)
survT5 <- survT5[1 : nbsurvT5 - 1, ]
survT5$m.av <- rollmean(survT5$qx, 15,fill = list(NA, NULL, NA))

survT6 <- subset(survrateT6, !is.nan(survrateT6[, 5]) & !is.na(survrateT6[, 1]))
nbsurvT6 <- nrow(survT6)
survT6 <- survT6[1 : nbsurvT6 - 1, ]
survT6$m.av <- rollmean(survT6$qx, 15,fill = list(NA, NULL, NA))

survT7 <- subset(survrateT7, !is.nan(survrateT7[, 5]) & !is.na(survrateT7[, 1]))
nbsurvT7 <- nrow(survT7)
survT7 <- survT7[1 : nbsurvT7 - 1, ]
survT7$m.av <- rollmean(survT7$qx, 15,fill = list(NA, NULL, NA))

survT8 <- subset(survrateT8, !is.nan(survrateT8[, 5]) & !is.na(survrateT8[, 1]))
nbsurvT8 <- nrow(survT8)
survT8 <- survT8[1 : nbsurvT8 - 1, ]
survT8$m.av <- rollmean(survT8$qx, 15,fill = list(NA, NULL, NA))

survT9 <- subset(survrateT9, !is.nan(survrateT9[, 5]) & !is.na(survrateT9[, 1]))
nbsurvT9 <- nrow(survT9)
survT9 <- survT9[1 : nbsurvT9 - 1, ]
survT9$m.av <- rollmean(survT9$qx, 15,fill = list(NA, NULL, NA))

survExp1 <- rbind(survT1, survT2, survT3, survT4, survT5, survT6, survT7, survT8, survT9)
write.csv(survExp1, file = ".\\longdata\\longdata-deathrates.csv")
survExp1 <- read.csv(filenames[1])

survExp1$age <- survExp1$age * 0.2

ggplot(survExp1, aes(age, qx, group = Treat, colour = Treat)) +
  geom_path(alpha = 0.5)

library(RColorBrewer)
library(colorspace)
pal <- choose_palette()

myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(dat$grp)
colScale <- scale_colour_manual(name = "grp",values = myColors)

ggplot(survExp1, aes(x=age, y=m.av, group=Treat)) +
  facet_wrap(~Treat, ncol = 3, scales = "fixed")+
  #geom_line(aes(color=Treat), size = 1) +
  #geom_point(shape = 21, size = 0.1, stroke = 0.1) +
  xlim(0, 300) +
  ylim(0, 0.075) +
  labs(x="Time (hour)", y = "qx") + 
  ggtitle(expression(atop("Mortality rate", atop(italic("total = 10240 Cells"), "")))) +
  stat_smooth(aes(color=Treat), method = "loess")
#+
 # facet_grid(. ~ prog, margins=TRUE)

  
  
  
  
