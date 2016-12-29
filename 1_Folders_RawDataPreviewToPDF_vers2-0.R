rm(list = ls())


# needed software for the pdf generation 
# https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/
# make sur to change the path the program in line 129 and 153

# needed libraries
library(data.table)
library(utils)

# will set the file location and working directory (library(utils) needed):
## to install, http://cran.r-project.org/web/packages/R.utils/index.html, 
#download: r-devel: R.utils_1.34.0.zip, install from zip (not cran)
wdpath <- choose.dir(getwd(), "Choose a suitable folder")
setwd(wdpath)
csv.data_dirs<-list.dirs(wdpath, recursive=FALSE)
dir.nb <- length(csv.data_dirs)




# loop for all folders##########################################################
################################################################################
k=1

# display progression bar
pba <- winProgressBar(title = "Quality Control: Cells and X Drift", min = 0, max = 100, width = 300)


for(k in seq(along = csv.data_dirs)) {

  # set the target folder
  cat(k)
  setwd(csv.data_dirs[k])
  cat(getwd())
  wdpath <-csv.data_dirs[k]
  csv.files <- list.files(pattern = ".csv")
  files.nb <- length(csv.files)
  
  
  # create output folders "PDF_xdrift", "PDF_cells" and "DataOnly"
  dir.create("PDF_xdrift")
  dir.create("PDF_cells")
  
  # path output folders "PDF_xdrift", "PDF_cells" and "DataOnly"
  PDF_xdrift_Path <- paste(dirname("PDF_xdrift"),"/","PDF_xdrift/",sep="")
  PDF_cells_Path <- paste(dirname("PDF_cells"),"/","PDF_cells/",sep="")





# loop for analysis ############################################################
################################################################################
  i <- 1
  for(i in 1 : files.nb) {

      # contruction of the file analysed
    currentfilename <- basename(csv.files[[i]])
    filenamenoext <- as.character(substring(currentfilename,
                                            1, nchar(currentfilename) - 4))
    cat(i)
    cat(" treated file:", currentfilename, "  -|-  ")
    
      # import the data for analysis
    
    # import the data for analysis
    raw <- read.csv(currentfilename, header=T, sep = "\t")
    testcol <- ncol(raw)
    # allows to avoid data misread between comma and tab separated files
    if (testcol == 1) {raw <- read.csv(currentfilename, header=T)}
    last_time_raw <- as.numeric(last(raw$Image.nb))
    # if  (k > 1){raw$Image.nb <- raw$Image.nb - last_time_raw +1}
    last_time <- as.numeric(last(raw$Image.nb))
    ###PLOTS####################################################################  
 
          # channel drift
    xdrift <- paste("xdrift-", filenamenoext, ".pdf", sep="")
    pdf(file=xdrift, width = 50, height = 10)
    mytitle = paste("Image name:", currentfilename)
    plot(raw$Image.nb , raw$x.chanel ,
         ylim = c(0,2560),
         main = mytitle,
         axes = FALSE,
         pch = ".")
    axis(side = 1, at = seq(0, last_time_raw, by = 50))
    axis(side = 2, at = seq(0, 2560, by = 100))
    abline(h = seq(0, 2560, by = 50), 
           v = seq(0, last_time_raw, by = 25), 
           col = "gray", lty = 3)
    dev.off()
    file.rename(from = paste(wdpath, "/", xdrift, sep = ""), 
                to = paste(wdpath, "/PDF_xdrift/", xdrift, sep = ""))
    
          # detected bacteria
    csv_cells_file <- paste("cells-", filenamenoext, ".pdf", sep = "")
    pdf(file=csv_cells_file, width = 50, height = 10)
    mytitle = paste("Image name:", currentfilename)
    plot(raw$Image.nb, raw$y, 
         main = mytitle,
         axes = FALSE,
         pch = "'")
    axis(side = 1, at = seq(0, last_time_raw, by = 50))
    axis(side = 2, at = seq(0, 800, by = 100))
    abline(h = seq(0, 800, by = 50), v = seq(0, last_time_raw, by = 25),
           col = "gray", lty = 3)
    dev.off()
    file.rename(from = paste(wdpath, "/", csv_cells_file, sep=""), 
                to = paste(wdpath, "/PDF_cells/", csv_cells_file, sep=""))
    
    ###END PLOTS################################################################  
}  
# end loop for analysis ########################################################

# list the files .pdf in the target folder
wdpath1<-paste(wdpath,"/PDF_cells/", sep="")
setwd(wdpath1)

## Collect the names of the figures to be glued together
csv.files.c <- list.files(pattern = ".pdf")

# contruction of the file analysed
namerecur <- substring(filenamenoext, 1, nchar(filenamenoext) - 9)
nameend <- substring(filenamenoext, (nchar(filenamenoext) - 3), 
                     nchar(filenamenoext))
PdfName <- paste("cells-", namerecur, nameend, sep = "")

## The name of the pdf doc that will contain all the figures
FileName <- paste(PdfName,"pdf",sep = ".")


## Make a system call to pdftk
## enter the path to the program pdftk
system2(command = "C:/Program Files (x86)/PDFtk/bin/pdftk.exe",
        args = c(shQuote(csv.files.c), "cat output", shQuote(FileName)))

## The command above is equiv. to typing the following at the system command line
## pdftk "*xy*.pdf" cat output "outFileName.pdf"
################################################################################

# list the files .pdf in the target folder
wdpath2<-paste(wdpath,"/PDF_xdrift/", sep="")
setwd(wdpath2)

## Collect the names of the figures to be glued together
csv.files.x <- list.files(pattern = ".pdf")

# contruction of the file analysed
namerecur <- substring(filenamenoext, 1, nchar(filenamenoext) - 9)
nameend <- substring(filenamenoext, (nchar(filenamenoext) - 3), 
                     nchar(filenamenoext))
PdfName <- paste("xdrift-", namerecur, nameend, sep = "")
FileName <- paste(PdfName,"pdf",sep=".")


## Make a system call to pdftk
## enter the path to the program pdftk
system2(command = "C:/Program Files (x86)/PDFtk/bin/pdftk.exe",
        args = c(shQuote(csv.files.x), "cat output", shQuote(FileName)))

# incrementation of the progression bar
info <- sprintf("%d%% done", round(k*(100/dir.nb)))
setWinProgressBar(pba, round(k*(100/dir.nb)), sprintf("Quality Control: Cells and X Drift (%s)", info), info)



## The command above is equiv. to typing the following at the system command line
## pdftk "*xy*.pdf" cat output "outFileName.pdf"
}

# close the progession bar
close(pba)

################################################################################
#### DONE ####### DONE ####### DONE ####### DONE ####### DONE ####### DONE #####


