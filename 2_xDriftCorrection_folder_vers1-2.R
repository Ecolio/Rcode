# version 1-2 addition of verification prompt allowing to set the parameters for each fields of views
# viosion 1-2 folder addition of loop through whole folder



rm(list = ls())

# call library
library(plyr)
library(utils)

# set working directory and list of file names 
wdpathGeneral <- choose.dir(getwd(), "Choose a suitable folder")
setwd(wdpathGeneral)


csv.data_dirs<-list.dirs(wdpathGeneral, recursive=FALSE)
cat(getwd())

# loop for all folders##########################################################
################################################################################
k = 1

for(k in seq(along = csv.data_dirs)) {
  
  j = 1
  files_full <- list.files(csv.data_dirs[k], full.names=T)
  files_name <- list.files(csv.data_dirs[k], full.names=F)
  nb <- length(files_full)
  files_full <- files_full[1:(nb-2)]
  files_name <- files_name[1:(nb-2)]
  
  # construct of required vectors
  tmp <- vector(mode = "list", length = length(files_full))

  # concanation of all csv files
  for (i in 1 : (nb - 2)) {
      tmp[[i]] <- read.csv(files_full[[i]], header=T, sep="\t")
      }

  # binding of all the data
  output <- do.call(rbind, tmp)
  tmpcor <- c()

  Max <- max(output$x.chanel)
  MaxTp <- as.numeric(max(output$Image.nb))
  MinTp <- as.numeric(min(output$Image.nb))
  
  x <- ""
  while(x == "") {
  ###### value that could be adapted ##########
  #############################################
  Min <- as.numeric(readline("lower limit : "))
  Modulo <- as.numeric(readline("modulo (around 109) : ")) 
  ab <- seq(Min, Max, Modulo) 
  ab <- c(ab, (Max + Modulo), (Max + Modulo), (Max + Modulo), (Max + Modulo))
  hist(output$x.chanel, breaks = 1901, xlim = c(0, (Max + 2 * Modulo)))
  abline(v = ab)
  x <- readline(" To correct, Press enter 
 To continue, Press 1 then Enter")  

  }
  
  # group all channel and export new corrected csv file
    for (j in 1 : (nb - 3)) {
      tmpcor[[j]] <- subset(output, 
                          x.chanel <= ab[[j + 1]] & 
                          x.chanel >= ab[[j]] )
      tmpcor[[j]] <- arrange( tmpcor[[j]], Image.nb, y)
      plot(xlim = c(ab[[j]], ab[[j + 1]]), ylim = c(MinTp,MaxTp),
           tmpcor[[j]]$x.chanel, tmpcor[[j]]$Image.nb, pch = "-")
      abline(v = ((ab[[j]] + ab[[j + 1]]) / 2))
      name <- paste("corrected", sep="_", files_name[[j]])
      namerecur <- substring(files_full[[j]], 1, nchar(files_full[[j]]) - 
                             nchar(files_name[[j]]))
      filename <- paste(namerecur, name, sep="")
      cat("------")
      cat(name)
      cat("------")
      write.csv(tmpcor[[j]], file = filename)
    }


}
################################################################################
#### DONE ####### DONE ####### DONE ####### DONE ####### DONE ####### DONE #####

#plot(output$x.chanel,output$Image.nb)
