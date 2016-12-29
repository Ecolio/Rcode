#### DEAD CELL DETECTION ####


rm(list = ls())

# needed libraries
library(data.table)
library(utils)
library(reshape)
library(zoo)


# Begining of the script
# setwd("C:/Users/lionel/Dropbox/Code")


################################################################################
################################################################################
# functions
insertRow          <- function(comp, newrow, radd) {
  
  ##############################################################################
  # Add the data to the compilation matrix   
  # --- FOR THE FIRST LINE  ---
  #
  # Args: - comp = appended matrix
  #       - nrow = new row to be added
  #       - radd = row level
  #
  # Returns: - comp = appended matrix
  ##############################################################################
  
  
  comp[seq(radd+1,nrow(comp)+1),] <- comp[seq(radd,nrow(comp)),]
  comp[radd,] <- newrow
  return(comp)
}

vecIn             <- function(vector,sequence){
  
  ##############################################################################
  # look for a sequence in a vector   
  # --- FOR THE FIRST LINE  ---
  #
  # Args: - vector = matrix row analysed
  #       - sequence = div/no div sequence reference for: - stop dividing
  #                                                       - dead
  #
  # Returns: - position of found sequences (first column of the sequence)
  ##############################################################################
  
  which(
    Reduce('+', lapply(seq_along(y <- lapply(sequence, '==', vector)), function(x){
      y[[x]][x:(length(vector) - length(sequence) +x)]
    }
    )
    ) == length(sequence)
  )
}


################################################################################
################################################################################


# will set the file location and worfoldering directory (library(utils) needed):
# to install, http://cran.r-project.org/web/pacfolderages/R.utils/index.html, 
# download: r-devel: R.utils_1.34.0.zip, install from zip (not cran)
wdpathGeneral <- choose.dir(getwd(), "Choose a suitable folder")
setwd(wdpathGeneral)

# create a list of all the folders in the experiment
csv.data_dirs<-list.dirs(wdpathGeneral, recursive=FALSE)
cat(getwd())

# extract the names of all the the folders
dirnames <- list.dirs(wdpathGeneral, recursive=F, full.names=F)
dirnamesanalysed <- list.dirs(wdpathGeneral, recursive=F, full.names=T)

# contruct the name of the matrix folders
namefoldermatrixb <- substring(dirnames[1], 1,10)
MatrixFolderName <- paste(namefoldermatrixb,"demo-parameters-fine", ".csv", sep = "")
MatrixFolderPath <- paste(wdpathGeneral, MatrixFolderName, sep = "/")  

################################################################################

dead.line.1.div <- 0
folderDiv <- 0.8                                                                # max coef of division
folderUncut <- 1.4                                                              # max coef of normal growth
Dead.cell <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0)

# create output files "cell.data" 
# .stopdiv: the time point at which the treatment stop cell growth (1, 2, 3)
# nb.div.: the number of division in between treatments (1, 2, 3)
# gr../3: the mean growth rate after treatments

cell.data <- as.data.frame((matrix(nrow=750, ncol=1750, NA)))
cell.data.line <- 0

# loop for all folders##########################################################
################################################################################
folders <- length(csv.data_dirs)
folder <- 1

for(folder in seq(along = csv.data_dirs)) {
  # set the target folder
  cat(folder)
  setwd(dirnamesanalysed[folder])
  cat(getwd())
  wdpathLocal <-csv.data_dirs[folder]
  csv.files <- list.files(pattern = ".csv")
  
  # import data
  Div.cells <- read.csv(csv.files[1], header = T)
  Fluo.cells <- read.csv(csv.files[3])
  Gr.cells <- read.csv(csv.files[4])
  L.cells <- read.csv(csv.files[5])
  NbRow <- nrow(Div.cells)
  NbCol <- ncol(Div.cells)  
  
  time.point <- min(which(is.na(Div.cells[7,])))
  time.point1 <- min(which(is.na(Div.cells[13,])))
  time.point2 <- min(which(is.na(Div.cells[19,])))
  
  if (time.point == time.point1 & time.point1 == time.point2) {
    
    time.point <- time.point
  } else {
    time.point <- max(c(time.point, time.point1, time.point2) , na.rm = T)
  }
  
  # remove the first column to get only the data
  Div.cells <- Div.cells[, 2 : time.point]
  Div.cells[is.na(Div.cells)] <- 0
  NbRow <- nrow(Div.cells)
  NbCol <- ncol(Div.cells)
  Div.cells <- Div.cells[seq(1,NbRow,6),]
  Fluo.cells <- Fluo.cells[, 2 : time.point]
  Fluo.cells[is.na(Fluo.cells)] <- 0
  Fluo.cells <- Fluo.cells[seq(1,NbRow,6),]
  Gr.cells <- Gr.cells[, 2 : time.point]
  Gr.cells <- Gr.cells[seq(1,NbRow,6),]
  Gr.cells <- Gr.cells[!(is.na(Gr.cells$V2)) & !(is.na(Gr.cells$V3)),]
  NbRowg <- nrow(Gr.cells)
  L.cells <- L.cells[1 : NbRowg,]
  Fluo.cells <- Fluo.cells[1 : NbRowg, ]
  Div.cells <- Div.cells[1 : NbRowg, ]

  NbRow <- nrow(Div.cells)
  NbCol <- ncol(Div.cells)
  
  
  # loop for analysis ##########################################################
  ##############################################################################
#  options(warn=1)
  channel <- 1
  results<-matrix(nrow=NbRow,ncol=1,NA)
  # loop through the matrix with modulo 6 (6 row of cells per channel)
  for(channel in 1 : NbRow) {
    cell.data.line <- cell.data.line + 1
    
    # test for the time of the last division with first treatment
    firstline <- Div.cells[channel, ]
    if (!is.na(Gr.cells[channel, 2])) {
   ##### if cells die ########################################################    

      # search for fluorescence detection for cell death of the first cell
      if (sum(Fluo.cells[channel,], na.rm = T) > 0) {
        
        fluo.death <- min(which(Fluo.cells[channel,] > 0))
        # define the number of column for the mean
        n <- 3
        gr.means <- c()
        # prepare the data to be meaned
        gr.rate.means <- t(Gr.cells[channel,])
        # replace the division rate by 1
        gr.rate.means <- replace(gr.rate.means, gr.rate.means <= 0.8 , NA)        
        # mean over n column (in theis case 3 column or time points)
        gr.means <- aggregate(gr.rate.means,list(rep(1:(nrow(gr.rate.means)%/%n+1),each=n,len=nrow(gr.rate.means))),mean,na.rm=T)[-1]
        # similar process with the division but with the sum of division in the interval
        div.sum <- c()
        div.sum <- t(Div.cells[channel,])
        div.sum <- aggregate(div.sum,list(rep(1:(nrow(div.sum)%/%n+1),each=n,len=nrow(div.sum))),sum)[-1]
        # contruct the output
        col.transfer <- round(fluo.death / 3) 
        gr.means <- as.data.frame(t(gr.means[1:col.transfer,]))
        div.sum <- as.data.frame(t(div.sum[1:col.transfer,]))
        cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - 2 - col.transfer), NA)))
        gr.means <- cbind(gr.means, cmpl)
        cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - col.transfer), NA)))
        div.sum <- cbind(div.sum, cmpl)
        
        ##### output lane with xy position, time of death and mean values and division rates
        compil <- as.data.frame(t(c(folder, fluo.death)))
        compil <- cbind(compil, gr.means, div.sum)
        
        
      # if no fluorecence de detected, then use the none division detection  
      } else { 
        
        dead.line.1.div <- min(vecIn(firstline, Dead.cell))
        #####  growth rates ###########
             if (dead.line.1.div == Inf ) {
               
        # define the number of column for the mean
                  n <- 3
                  gr.means <- c()
        # prepare the data to be meaned
                  gr.rate.means <- t(Gr.cells[channel,])
        # replace the division rate by 1
                  gr.rate.means <- replace(gr.rate.means, gr.rate.means <= 0.8 , NA)        
        # mean over n column (in theis case 3 column or time points)
                  gr.means <- aggregate(gr.rate.means,list(rep(1:(nrow(gr.rate.means)%/%n+1),each=n,len=nrow(gr.rate.means))),mean,na.rm=T)[-1]
        # similar process with the division but with the sum of division in the interval
                  div.sum <- c()
                  div.sum <- t(Div.cells[channel,])
                  div.sum <- aggregate(div.sum,list(rep(1:(nrow(div.sum)%/%n+1),each=n,len=nrow(div.sum))),sum)[-1]
        # contruct the output
                  col.transfer <- round(NbCol / 3) 
                  gr.means <- as.data.frame(t(gr.means[1:col.transfer,]))
                  div.sum <- as.data.frame(t(div.sum[1:col.transfer,]))
                  cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - 2 - col.transfer), NA)))
                  gr.means <- cbind(gr.means, cmpl)
                  cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - col.transfer), NA)))
                  div.sum <- cbind(div.sum, cmpl)
                  
                  #####  growth rates ###########
                  compil <- as.data.frame(t(c(folder, NbCol)))
                  compil <- cbind(compil, gr.means, div.sum)

                } else {

        # define the number of column for the mean
                  n <- 3
                  gr.means <- c()
        # prepare the data to be meaned
                  gr.rate.means <- t(Gr.cells[channel,])
        # replace the division rate by 1
                  gr.rate.means <- replace(gr.rate.means, gr.rate.means <= 0.8 , NA)        
        # mean over n column (in theis case 3 column or time points)
                  gr.means <- aggregate(gr.rate.means,list(rep(1:(nrow(gr.rate.means)%/%n+1),each=n,len=nrow(gr.rate.means))),mean,na.rm=T)[-1]
        # similar process with the division but with the sum of division in the interval
                  div.sum <- c()
                  div.sum <- t(Div.cells[channel,])
                  div.sum <- aggregate(div.sum,list(rep(1:(nrow(div.sum)%/%n+1),each=n,len=nrow(div.sum))),sum)[-1]
        # contruct the output
                  col.transfer <- round(dead.line.1.div / 3) 
                  gr.means <- as.data.frame(t(gr.means[1:col.transfer,]))
                  div.sum <- as.data.frame(t(div.sum[1:col.transfer,]))
                  cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - 2 - col.transfer), NA)))
                  gr.means <- cbind(gr.means, cmpl)
                  cmpl <- as.data.frame((matrix(nrow = 1, ncol = (750 - col.transfer), NA)))
                  div.sum <- cbind(div.sum, cmpl)
                  
                  #####  growth rates ###########
                  compil <- as.data.frame(t(c(folder, dead.line.1.div)))
                  compil <- cbind(compil, gr.means, div.sum)

                }
      }
    }
    cell.data[cell.data.line, ] <- compil

  }

}
    
    write.csv(cell.data, file = MatrixFolderPath)
      
      
      