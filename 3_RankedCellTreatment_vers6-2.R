# version control and comments
# vers1.x : development phase
# vers2.EvoDemo : multiple files analysed and export of the results, first cell
# vers3.1 : progress bar
# vers3.2 : peaks detection for grey levels, left aside
# vers4.0 : New function for line analysis, Change 0 to NA 
# vers4.4 : intergration of gap removal for cells that are in the middle of the channels
# vers5.0 : vertical analysis of the matrix
# vers5.2 : automatic adjustment of input file format (, or tabbulated)
# vers6.0 : multiple row analysis, not only first cell
# 3_RankedCellTreatment_vers6-2.R

rm(list = ls())

# needed libraries
library(data.table)
library(utils)
library(reshape)
library(miscTools)
library(BBmisc)
library(plyr)



# TO BE ACTIVATED IF RUNNING THROUGH SOURCE 
 setwd(dirname(sys.frame(1)$ofile))

# TO BE ACTIVATED IF RUNNING LINE BY LINE OR THROUGH RUN
# setwd("C:\\Users\\lionel\\Dropbox\\Code")

# functions
source("function\\RankedCellTreatmentfunction_vers6-2.R") 

# will set the file location and working directory (library(utils) needed):
## to install, http://cran.r-project.org/web/packages/R.utils/index.html, 
#download: r-devel: R.utils_1.34.0.zip, install from zip (not cran)
# the supressMessages function allows to avoid to display the error messages due
#to the choose.dir function
suppressMessages( wdpathGeneral <- choose.dir(getwd(), "Choose a suitable folder"))
setwd(wdpathGeneral)
csv.data_dirs <- list.dirs(wdpathGeneral, recursive=FALSE)
cat(getwd())
dir.nb <- length(csv.data_dirs)
dirnames <- list.dirs(wdpathGeneral, recursive=F, full.names=F)

# loop for all folders##########################################################
################################################################################
folder.analysed = 1

# display progression bar
pba <- winProgressBar(title = "Matrices Extraction", min = 0, max = 100, width = 300)

for(folder.analysed in seq(along = csv.data_dirs)) {
  
  # set the target folder
  cat(folder.analysed)                                                          # display cycle number
  setwd(csv.data_dirs[folder.analysed])                                         # set folder treated
  cat(getwd())                                                                  # display folder treated
  wdpathLocal <-csv.data_dirs[folder.analysed]                                  # folder treated
  csv.files <- list.files(pattern = ".csv")                                     # list of .csv files
  files.nb <- length(csv.files)                                                 # length of list of .csv files
  
  # contruct the name of the matrix folders
  # begining of the name
  namefoldermatrixb <- substring(dirnames[folder.analysed], 1, nchar(dirnames[folder.analysed]) - 12)
  # ending of the name
  namefoldermatrixe <- substring(dirnames[folder.analysed], (nchar(dirnames[folder.analysed]) - 3), 
                                 nchar(dirnames[folder.analysed]))
  # contruction of the folder name
  MatrixFolderName <- paste("matrix_", namefoldermatrixb,
                            namefoldermatrixe, sep = "")
  # creation of the matrix destination folder
  MatrixFolderPath <- paste(wdpathGeneral, MatrixFolderName, sep = "/")  
  
  # create output folders "matrix_xy-n"
  dir.create(MatrixFolderPath)
  
  # create output files "parameter-compilation" l:length, div:division, f:fluorescence, gr:growth rate
  l.comp <- as.data.frame((matrix(nrow = 300, ncol = 2000)))
  div.comp <- as.data.frame((matrix(nrow = 300, ncol = 2000)))
  f.comp <- as.data.frame((matrix(nrow = 300, ncol = 2000)))
  gr.comp <- as.data.frame((matrix(nrow = 300, ncol = 2000)))
  empty.file.list <- as.data.frame((matrix(nrow = 14, ncol = 1,0)))
  radd <- 1

  # loop for analysis ############################################################
  ################################################################################
  file.analysed <- 1
  for(file.analysed in 1 : files.nb) {
    
    # contruction of the file analysed
    currentfilename <- basename(csv.files[[file.analysed]])
    # begining of the name
    filenamenoext <- as.character(substring(currentfilename,
                                            1, nchar(currentfilename) - 16))
    # display cycle number in the XY folder
    cat(file.analysed)                            
    # display file treated
    cat(" treated file:", currentfilename, "  -|-  ")
    
    # import the data for analysis
    raw <- read.csv(currentfilename, header=T, sep = "\t")
    testcol <- ncol(raw)
    testraw <- nrow(raw)
    # allows to avoid data misread between comma and tab separated files
    if (testcol == 1) {raw <- read.csv(currentfilename, header=T)}
    if (testraw < 10) { cat("empty file")
    } else {
      first_time_raw <- as.numeric(raw$Image.nb[1])
    
    # correct misnumbering of the time (if no misnumbering, need to be hashed)
    if (first_time_raw > 1){raw$Image.nb <- (raw$Image.nb - first_time_raw) + 1}
    last_time <- as.numeric(last(raw$Image.nb))
    
    # add the limit of the dead end of the side channel if not detected
    y.position <- 1
    cat("missing time point ")
    for (y.position in 1 : last_time) {
      first.y <- min(which(raw$Image.nb == y.position))
      if (first.y == Inf) { 
        cat(y.position)
      }
      else if (y.position == 1) {
        newrowinsert <- raw[first.y,]
        newrowinsert$y <- 3
        raw <- rbind(newrowinsert, 
                     raw[(first.y + 1) : nrow(raw),] )
      }
      else if (raw$y[[first.y]] > 10) {
        newrowinsert <- raw[first.y,]
        newrowinsert$y <- 3
        raw <- rbind(raw[ 1 : (first.y - 1), ], newrowinsert, 
                     raw[first.y : nrow(raw),] )
      }
      
    }
    
    
    raw <- raw[, c("Image.nb",                                                  # define the column name of the raw data
                   "y",
                   "grey.slope",
                   "grey.before",
                   "grey.after",
                   "x.chanel",
                   "chanel.end",
                   "gravity.center.color1",
                   "surface.color1",
                   "surface.mean.color1",
                   "surface.threshold1",
                   "gravity.center.color2",
                   "surface.color2",
                   "surface.mean.color2",
                   "surface.threshold2")]
    
    
    # Organise the data ########################################################

    # Automatic detection of grey levels to set the exclusion levels
    grey.value.gap.remover.int <- quantile(raw$grey.before)                     # determine the spread of grey values at detection points
    grey.value.gap.remover <- mean(tail(grey.value.gap.remover.int, n = 3))     
    grey.value.dark <- (as.numeric(grey.value.gap.remover.int[2]) - 
                          (as.numeric(grey.value.gap.remover.int[2])/10))       # formula to determine cut off point to exclude dark points in bacteria
    raw.clean <- c()
    image.number <- 1
    for (image.number in 1 : last_time) {
      temp <- subset(raw, raw$Image.nb == image.number)
      gap.remover <- 1
      rowtemp <- nrow(temp)
      if (rowtemp < 2 ) { 
        cat(image.number)
      } else {
        value.to.be.substracted <- 0
        for (gap.remover in 1 : (rowtemp)) {
          tempor <- temp[gap.remover,]
          if (gap.remover == 1) {
            temporary <- tempor
          }
          else if (temp$grey.after[gap.remover - 1] > grey.value.gap.remover &
                     temp$grey.before[gap.remover] > grey.value.gap.remover) {
            value.to.be.substracted <- (temp$y[gap.remover] - 
                                        temp$y[gap.remover - 1]) + 
                                        value.to.be.substracted
          }
          else if (temp$grey.after[gap.remover - 1] < grey.value.dark &
                     temp$grey.before[gap.remover] < grey.value.dark) {
            temporary <- temporary[1:(gap.remover - 1), ]
          } else {
            tempor$y <- tempor$y - value.to.be.substracted
            temporary <- rbind(temporary, tempor)
          }
        } 
        raw.clean <- rbind(raw.clean, temporary)
      }
    }
    raw.clean <- na.omit(raw.clean)
    raw.clean$weighting <- 0

    

    
    # bring back all the cells to the bottom of the tube
    slope.cutoff.int <- quantile(abs(raw.clean$grey.slope))
    slope.cutoff <- mean(slope.cutoff.int[2:3])
    #raw.clean <- subset(raw.clean, abs(raw.clean$grey.slope) > slope.cutoff)
    raw.clean <- subset(raw.clean, raw.clean$y > 50)
    
    y.position <- 1
    for (y.position in 1 : last_time) {
      first.y <- min(which(raw.clean$Image.nb == y.position))
      if (first.y == Inf) { 
        cat(y.position)
      }
      else if (first.y == 1) {
        newrowinsert <- raw.clean[first.y,]
        newrowinsert$y <- 3
        raw.clean <- rbind(newrowinsert, 
                           raw.clean[(first.y) : nrow(raw.clean),] )
      }
      else if (raw.clean$y[[first.y]] > 10) {
        newrowinsert <- raw.clean[first.y,]
        newrowinsert$y <- 3
        raw.clean <- rbind(raw.clean[ 1 : (first.y - 1), ], newrowinsert, 
                           raw.clean[first.y : nrow(raw.clean),] )
      }
      
    }  
    
    
    
    # shift the y column to calculate the cell length ######################
    shift <- function(x, n){
      c(x[ - (seq(n))], rep(NA, n))
    }
    
    if (!(is.data.frame(raw.clean) && nrow(raw.clean) == 0)){
      raw.clean$y_y.1 <- shift(raw.clean$y, 1)
      raw.clean$cell.length <- raw.clean$y_y.1 - raw.clean$y
    }

    # remove double/triple points and weighting the result ###########################
    raw.clean.select <- raw.clean
    raw.clean.select[!is.na(raw.clean.select)] <- 0
    no.data.row <- raw.clean.select[1, ]
    no.data.row[!is.na(no.data.row)] <- NA
    row.number <- 1
    for (row.number in 1 : nrow(raw.clean - 2)) {
        if (!is.na(raw.clean.select$y[row.number]) &&
            !is.na(raw.clean.select$cell.length[row.number])  ) {
              if (raw.clean$cell.length[row.number] < 0) {
                raw.clean.select[row.number, ] <- raw.clean[row.number, ]
                raw.clean.select$weighting[row.number] <- 1
              }
              else if (raw.clean$y[row.number] < 10  ) {
                  raw.clean.select[row.number, ] <- raw.clean[row.number, ]
                  raw.clean.select$weighting[row.number] <- 1
              }
              else if ((raw.clean$cell.length[row.number] +
                  raw.clean$cell.length[row.number + 1])  < 45 &&
                  raw.clean$cell.length[row.number] > 0 &&
                  !is.na(raw.clean$cell.length[row.number + 1]) ) {
                  raw.clean.select[row.number, ] <- raw.clean[row.number, ]
                  raw.clean.select[row.number + 1, ] <- no.data.row
                  raw.clean.select$weighting[row.number] <- 1
              }
              else if (raw.clean$cell.length[row.number] < 20 &&
                  raw.clean$cell.length[row.number] > 0 &&
                  !is.na(raw.clean$cell.length[row.number])) {
                  raw.clean.select[row.number, ] <- raw.clean[row.number, ]
                  raw.clean.select[row.number + 1, ] <- no.data.row
                  raw.clean.select$weighting[row.number] <- 1
                
              } else { 
                raw.clean.select[row.number, ] <- raw.clean[row.number, ]
              }
        } else {
          raw.clean.select[row.number, ] <-no.data.row
        }
    }  
    
    raw.clean.select[((nrow(raw.clean) - 1) : nrow(raw.clean)), ] <-
        raw.clean[((nrow(raw.clean) - 1) : nrow(raw.clean)), ]
    raw.clean.select <- na.omit(raw.clean.select)
    
    if (!(is.data.frame(raw.clean.select) && nrow(raw.clean.select) == 0)){
      raw.clean.select$y_y.1 <- shift(raw.clean.select$y, 1)
      raw.clean.select$cell.length <- raw.clean.select$y_y.1 - 
                                      raw.clean.select$y
    }  
    
    # testlowvalues <- subset(raw.clean.select, raw.clean.select$y < 50)
    kCellLenghtLLInt <- quantile(raw.clean.select$cell.length, na.rm = T)
    kCellLenghtLL <- as.numeric(kCellLenghtLLInt[2] + 
                               (kCellLenghtLLInt[2] / 2))
    if (kCellLenghtLL < 30 && !is.na(kCellLenghtLL)) {kCellLenghtLL <- 30}
    cells <- subset(raw.clean.select, 
                    raw.clean.select$cell.length > kCellLenghtLL)
################################################################################
# used mainly for data checking in case of abberant analysis
      # plot(raw$Image.nb, raw$y, xlim = c(0,100), ylim = c(0,600), pch = "-")
      # plot(raw.clean$Image.nb, raw.clean$y, xlim = c(0,100), ylim = c(0,600), pch = "-")
      # plot(raw.clean.select$Image.nb, raw.clean.select$y, xlim = c(0,100), ylim = c(0,600), pch = "-")
      # plot(cells$Image.nb, cells$y, xlim = c(0,100), ylim = c(0,600), pch = "-")
################################################################################
    if ( nrow(cells) <= 10 | is.na(kCellLenghtLL)) { 
      cat("empty tube") 
      EmptyTube <- as.data.frame(matrix(nrow = 6, ncol = 2000))
      l.comp <- rbind(l.comp, EmptyTube)
      div.comp <- rbind(div.comp, EmptyTube)
      f.comp <- rbind(f.comp, EmptyTube)
      gr.comp <- rbind(gr.comp, EmptyTube)
      empty.file.list <- rbind(empty.file.list,currentfilename)
    } else {
      # rank the cells in regard to y
      first_cells <- transform(cells, 
                               cell.rank = ave(y, Image.nb  , 
                                               FUN = function(x) rank(x, ties.method = "first")))
      # reshape in a matrix for cell length, point weighting (.p) and fluorescences
      matrix.length <- cast(first_cells, cell.rank ~ Image.nb, 
                            value = 'cell.length')
      matrix.length.p <- cast(first_cells, cell.rank ~ Image.nb, 
                              value = 'weighting')
      matrix.length.f1 <- cast(first_cells, cell.rank ~ Image.nb, 
                              value = 'surface.mean.color1')
      matrix.length.f2 <- cast(first_cells, cell.rank ~ Image.nb, 
                              value = 'surface.mean.color2')
      # make it a data.frame and give cols the same names as data
      newrow <- matrix(nrow = nrow(matrix.length), ncol = ncol(matrix.length))
      colnames(newrow) <- colnames(matrix.length)
      # rbind the empty row to matrix.length and expend the ranking to last row
      matrix.length <- rbind(matrix.length,newrow)
      matrix.length.p <- rbind(matrix.length.p,newrow)
      matrix.length.f1 <- rbind(matrix.length.f1,newrow)
      matrix.length.f2 <- rbind(matrix.length.f2,newrow)
      r <- as.numeric(nrow(matrix.length))
      matrix.length$cell.rank <- c(1 : r)
      # extract only the Matrix.Data as Matrix.Data-frame and create composition matrix 
      m.dat.l <- matrix.length[, -1]
      m.dat.p <- matrix.length.p[, -1]
      m.dat.f1 <- matrix.length.f1[, -1]
      m.dat.f2 <- matrix.length.f2[, -1]
      
      # variables 
      nbrow <- nrow(m.dat.l)
      tp <- tp.f <- 1                                                       # time variable 
      rp <- rp.f <- 1                                                       # row number
      fp <- 1                                                                   # first time point variable
      shift <- 0                                                                # shift due to cell re-cut
      div <- 0                                                                  # mark div on the lower ranked cell only
      cor.length <- 0
      gr.cor.pred <- 1
      fill.forward <- NA
      row.empty <- NA
      firstposition <- 1
      m.cor.value <- NA      
      
      # constant for cell correction 
      kDiv <- 0.8                                                               # max coef of division
      kUncut <- 1.4                                                             # max coef of normal growth
      kShrink <- 0.9                                                            # mac coef of normal shrinkage
      kGR <- 1.1                                                                # coef of normal growth
      kPixelForLongCell <- 200                                                  # max lenght of normal cell in pixel
      kPixelNewCell <-                                                          # min lenght of normal cell in pixel
      kOrd <- 10000                                                             # value added to allow ordering of the matrix
      
      # preparation of the matrix 
      ColNmall <- as.numeric(colnames(m.dat.l)) #sub("X", "", ))                # collect data matrix column number
      virtual <- seq.default(1, max(ColNmall))    #min(ColNmall)                # generate data matrix column sequence number
      missg <- virtual[!virtual %in% ColNmall] + kOrd                           # find missing column and add kOrd (10000) for later ordering
      missing <- matrix(nrow = nbrow, ncol = length(missg))                            
      missingp <- matrix(nrow = nbrow, ncol = length(missg))                            
      ColNmall <- ColNmall + kOrd
      colnames(m.dat.l) <- ColNmall                                             # rename the columns for ordering
      colnames(m.dat.p) <- ColNmall                                             # rename the columns for ordering
      colnames(m.dat.f1) <- ColNmall                                            # rename the columns for ordering
      colnames(m.dat.f2) <- ColNmall                                            # rename the columns for ordering
      colnames(missing) <- missg
      colnames(missingp) <- missg
      
      missingp[,] <- 0
      missing[,] <- NA
      
      m.dat.l <- cbind(m.dat.l, missing)                                        # reconstruct the completed matrix
      m.dat.l <- m.dat.l[,order(colnames(m.dat.l))]
      m.dat.p <- cbind(m.dat.p, missingp)                                       
      m.dat.p <- m.dat.p[,order(colnames(m.dat.p))]
      m.dat.f1 <- cbind(m.dat.f1, missing)                                      
      m.dat.f1 <- m.dat.f1[,order(colnames(m.dat.f1))]
      m.dat.f2 <- cbind(m.dat.f2, missing)                                      
      m.dat.f2 <- m.dat.f2[,order(colnames(m.dat.f2))]

      m.div <- as.data.frame(apply(m.dat.l, 2, FUN = function(x) x * NA))                                                                 # output of div matrix
      m.shift <- m.div                                                          # keep track of the cell shift
      m.cor <- m.div                                                            # output of m.dat.l corrected
      m.gr.raw <- m.div                                                         # output of m.dat.l grow rate
      m.gr.cor <- m.gr.raw                                                      # output of gr.cor
      
      
      nb.tp <- as.numeric(length(m.dat.l[1, ]))                                 # column number(time)
      
      
      ################################################################################
      ################################################################################
      
      options(warn=1)
      
      incremTooSmall <- 0
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      # Analysis of the 1st row 
      for (tp in fp : (nb.tp - 1)) {                                            # loop through the first line, from t=2 til last time point
        tp.ff <- tp.f <- tp
        
        if (is.na(m.dat.l[rp, tp])) {                                           # Correct for empty time point
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else { 
          m.cor[rp, tp] <- RowEmpty(rp, tp)                                     #--- Function RowEmpty 
        } 
        
        if (is.na(m.dat.f1[rp, tp]) &&
            tp > 2) {
          m.dat.f1[rp, tp] <- m.dat.f1[rp, (tp - 1)]
        }
        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw
        gr.cor.pred <- GRCorPred(rp, tp)
        analysis <- FirstLine(rp, tp, gr.cor.pred)                              #--- Function FirstLine 
        ResultsWritting(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
        
      }
      
      # Analysis of the 2nd row 
      
      rp <- rp.f + 1
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      for (tp in fp : (nb.tp - 1)) {            
        tp.ff <- tp.f <- tp
        shift <- 0
        div <- 0
        cor.length <- 0
        analysis <- c()
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)                          # rank after correction of uncut
        div <- sum(m.div[, tp], na.rm = T)                                      # rank after cell division
        
        if (is.na(m.dat.l[rp, tp])) {                                           # Correct for empty time point
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else { 
          if (m.cor[rp, tp] > 45 && !is.na(m.cor[rp, tp])) {
             m.cor[rp + 1, tp] <- RowEmpty(rp, tp) 
          } else {
             m.cor[rp, tp] <- RowEmpty(rp, tp)                                  #--- Function RowEmpty 
          }
        } 
        
        if (is.na(m.dat.f1[rp, tp]) &&
            tp > 2) {
          m.dat.f1[rp, tp] <- m.dat.f1[rp, (tp - 1)]
        }        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw 
        gr.cor.n <- GRCorPred(rp, tp)                                           #--- Function GRCorPrev 
        analysis <- OtherLines(rp, tp, gr.cor.n)                                #--- Function OtherLines 
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)
        ResultsWrittingN(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
      }
     
      # Analysis of the 3rd row 
      rp <- rp.f + 1
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      for (tp in fp : (nb.tp - 1)) {            
        shift <- 0
        div <- 0
        cor.length <- 0
        analysis <- c()
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)                          # rank after correction of uncut
        div <- sum(m.div[, tp], na.rm = T)                                      # rank after cell division
        
        if (is.na(m.dat.l[rp, tp])) {       
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else {                                                  
           if (!is.na(m.cor[rp, tp])) {
             if (m.cor[rp, tp] > 45 && is.na(m.cor[rp + 1, tp])) {
             m.cor[rp + 1, tp] <- RowEmpty(rp, tp) 
             } 
             else if (m.cor[rp, tp] > 45 
                   && m.cor[rp + 1, tp] > 45 
                   && is.na(m.cor[rp + 2, tp])) { 
             m.cor[rp + 2, tp] <- RowEmpty(rp, tp) 
             }
         } else {
             m.cor[rp, tp] <- RowEmpty(rp, tp)                                  #--- Function RowEmpty 
          }                                    
        }
        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw 
        gr.cor.n <- GRCorPred(rp, tp)                                           #--- Function GRCorPrev 
        analysis <- OtherLines(rp, tp, gr.cor.n)                                #--- Function OtherLines 
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)
        ResultsWrittingN(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
      }
      
      # Analysis of the 4th row 
      rp <- rp.f + 1
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      for (tp in fp : (nb.tp - 1)) {            
        tp.ff <- tp.f <- tp
        shift <- 0
        div <- 0
        cor.length <- 0
        analysis <- c()
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)                                   # rank after correction of uncut
        div <- sum(m.div[, tp], na.rm = T)                                       # rank after cell division
        
        if (is.na(m.dat.l[rp, tp])) {       
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else {                                                  
          if (!is.na(m.cor[rp, tp])) {
             if (m.cor[rp, tp] > 45 && is.na(m.cor[rp + 1, tp])) {
             m.cor[rp + 1, tp] <- RowEmpty(rp, tp) 
             } 
             else if (m.cor[rp, tp] > 45 
                   && m.cor[rp + 1, tp] > 45 
                   && is.na(m.cor[rp + 2, tp])) { 
             m.cor[rp + 2, tp] <- RowEmpty(rp, tp) 
             }
             else if (m.cor[rp, tp] > 45 
                   && m.cor[rp + 1, tp] > 45 
                   && m.cor[rp + 2, tp] > 45 
                   && is.na(m.cor[rp + 3, tp])) { 
             m.cor[rp + 2, tp] <- RowEmpty(rp, tp) 
             }
         } else {
             m.cor[rp, tp] <- RowEmpty(rp, tp)                                  #--- Function RowEmpty 
          }
        }
        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw 
        gr.cor.n <- GRCorPred(rp, tp)                                           #--- Function GRCorPrev 
        analysis <- OtherLines(rp, tp, gr.cor.n)                                #--- Function OtherLines 
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)
        ResultsWrittingN(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
      }
      
      
      # Analysis of the 5th row 
      rp <- rp.f + 1
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      for (tp in fp : (nb.tp - 1)) {            
        # reset row variables
        #tp<-tp+1
        tp.ff <- tp.f <- tp
        shift <- 0
        div <- 0
        cor.length <- 0
        analysis <- c()
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)                                   # rank after correction of uncut
        div <- sum(m.div[, tp], na.rm = T)                                       # rank after cell division
        
        if (is.na(m.dat.l[rp, tp])) {       
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else {                                                  
          m.cor[rp, tp] <- RowEmpty(rp, tp)                                     #--- Function FillForward
        }
        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw 
        gr.cor.n <- GRCorPred(rp, tp)                                           #--- Function GRCorPrev 
        analysis <- OtherLines(rp, tp, gr.cor.n)                                #--- Function OtherLines 
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)
        ResultsWrittingN(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
      }
      
      
      # Analysis of the 6th row 
      rp <- rp.f + 1
      rp.ff <- rp.f <- rp
      fp <- tp.ff <- tp.f <- tp <- 1
      for (tp in fp : (nb.tp - 1)) {            
        # reset row variables
        #tp<-tp+1
        tp.ff <- tp.f <- tp
        shift <- 0
        div <- 0
        cor.length <- 0
        analysis <- c()
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)                                   # rank after correction of uncut
        div <- sum(m.div[, tp], na.rm = T)                                       # rank after cell division
        
        if (is.na(m.dat.l[rp, tp])) {       
          m.cor[rp, tp] <- FillForward(rp, tp, gr.cor.pred)                     #--- Function FillForward
          m.dat.p[rp, tp] <- 0
        } else {                                                  
          m.cor[rp, tp] <- RowEmpty(rp, tp)                                     #--- Function FillForward
        }
        
        m.gr.raw[rp, tp] <- GRRawF(rp, tp)                                      #--- Function gr.raw 
        gr.cor.n <- GRCorPred(rp, tp)                                           #--- Function GRCorPrev 
        analysis <- OtherLines(rp, tp, gr.cor.n)                                #--- Function OtherLines 
        shift <- sum(m.shift[(rp - 1), tp], na.rm = T)
        ResultsWrittingN(rp, tp)
        m.gr.cor[rp, tp] <- CorrectGRF(rp, tp)                                  # Corrected growth rate 
      }
      
      nbcol <- ncol(m.cor)
      add.comp <- as.data.frame(matrix(nrow = 6, ncol = (2000 - nbcol)))
      newrow1 <- m.cor[1 : 6, ]
      newrow1 <- suppressWarnings(cbind(newrow1, add.comp))
      newrow1 <- as.data.frame(newrow1)
      l.comp[(radd : (radd + 5)), ] <- newrow1
      newrow2 <- m.div[1 : 6, ] 
      newrow2 <- suppressWarnings(cbind(newrow2, add.comp))
      newrow2 <- as.data.frame(newrow2)
      div.comp[(radd : (radd + 5)), ] <- newrow2
      newrow3 <- m.dat.f1[1 : 6, ]
      newrow3 <- suppressWarnings(cbind(newrow3, add.comp))
      newrow3 <- as.data.frame(newrow3)
      f.comp[(radd : (radd + 5)), ] <- newrow3
      newrow4 <- m.gr.cor[1 : 6, ]
      newrow4 <- suppressWarnings(cbind(newrow4, add.comp))
      newrow4 <- as.data.frame(newrow4)
      gr.comp[(radd : (radd + 5)), ] <- newrow4

#      newrow1 <- m.cor[1, ]
#      l.comp <- suppressWarnings(insertRow(l.comp, radd, newrow1))
#      newrow2 <- m.div[1, ]
#      div.comp <- suppressWarnings(insertRow(div.comp, radd, newrow2))
#      newrow3 <- m.dat.f1[1, ]
#      f.comp <- suppressWarnings(insertRow(f.comp, radd, newrow3))
#      newrow4 <- m.gr.cor[1, ]
#      gr.comp <- suppressWarnings(insertRow(gr.comp, radd, newrow4))
      
      
    }
    radd <- radd + 6
    }
  }
  # end loop for analysis ########################################################
  
  fname <- paste("fluo1_", MatrixFolderName, ".csv", sep = "")
  lname <- paste("length_", MatrixFolderName, ".csv", sep = "")
  divname <- paste("div_", MatrixFolderName, ".csv", sep = "")
  grname <- paste("gr_", MatrixFolderName, ".csv", sep = "")
  empties <- paste("emptyfiles_", MatrixFolderName, ".csv", sep = "")
  
  setwd(MatrixFolderPath)
  write.csv(l.comp, file = lname)
  write.csv(div.comp, file = divname)
  write.csv(f.comp, file = fname)
  write.csv(gr.comp, file = grname)
  write.csv(empty.file.list, file = empties)
  
  
  
  
  
  # incrementation of the progression bar
  info <- sprintf("%d%% done", round(folder.analysed*(100/dir.nb)))
  setWinProgressBar(pba, round(folder.analysed*(100/dir.nb)), sprintf("Matrices Extraction (%s)", info), info)
  
  
} 

# close the progession bar
close(pba)

################################################################################
#### DONE ####### DONE ####### DONE ####### DONE ####### DONE ####### DONE #####
