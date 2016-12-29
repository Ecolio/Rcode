# Version 1.2 -> inclusion of forward testing:

#    - Time control
#    - Time insertion
#    - Forward probing
#    - GRcorrection on t-5 -> t-1
# all analysis are done in the m.cor, not the m.dat.l

# Version 2.0 -> foward completion of missing time points

# Version 3.0 -> exportation of results to files
# Version 3.1 -> removed unsed fonctions
# Version 3.2 -> integration of NA and corrected bugs of growth
# Version 3.4 -> growth predictor simplified
# Version 4.0 -> function for the 1+lines for easy reading
# Version 5.0 -> reading of t+1 for artefact correction
# Version 6.0 -> weighting of the point



################################################################################

# COMMON FUNCTIONS 

# corrected growth rates

GRRawF             <- function(rp.f, tp.f) {
  

  div <- sum(m.div[, tp.f], na.rm = T) 
  if (tp.f <= 2) {
    gr.raw <- kGR }
  else if (!is.na(m.cor[rp.f, tp.f])) {
      if (rp.f >= 2) {
          if ((m.div[(rp.f - 1), tp.f] == 1)  && 
              !is.na(m.div[(rp.f - 1), tp.f] == 1)){
            gr.raw <- m.gr.cor[(rp.f - 1), tp.f]
          } else {    
            gr.raw <- round(as.numeric((m.cor[rp.f, tp.f] /
                                         m.cor[rp.f-div, tp.f - 1])), 2)
          }
      } else if (rp.f < 2){    
        gr.raw <- round(as.numeric((m.cor[rp.f, tp.f] /
                                     m.cor[rp.f, tp.f - 1])), 2)                                   
      }
  } else {
    gr.raw <- kGR
  }
  
  if ((gr.raw >= 0.1) & (gr.raw <= 10^10) & !is.na(gr.raw)) {
    gr.raw <- gr.raw 
  } else { gr.raw <- kGR }
  
  return(gr.raw)
}

#_______________________________________________________________________________


#_______________________________________________________________________________

CorrectGRF         <- function(rp.f, tp.f) {

  div <- sum(m.div[, tp.f], na.rm = T) 
  
  if (rp.f == 1) {
    if (tp.f == 1) {
      gr.cor <- kGR
      return(gr.cor)
    } else {
      gr.cor <- round(as.numeric((m.cor[rp.f, tp.f] / 
                                m.cor[rp.f, (tp.f - 1)])), 2)
      return(gr.cor)
    }
  } 
  else if (rp.f > 1 & tp.f <= 3) {
    gr.cor <- kGR
    return(gr.cor)
  } 
  else if (is.na(m.cor[rp.f, tp.f])) {
    gr.cor <- NA
    return(gr.cor)
  }
  else if (is.na(m.cor[(rp.f - div), (tp.f - 1)])) {
    gr.cor <- NA
    return(gr.cor)
  } else {
    gr.cor <- round(as.numeric(m.cor[rp.f, tp.f] / 
                               m.cor[(rp.f - div), (tp.f - 1)]), 2)
  return(gr.cor)  
  }
}


#_______________________________________________________________________________
#|______________________________________________________________________________|



# FORWARD FUNCTIONS 

# corrected growth rates form t-3 -> t-1

GRCorPred          <- function(rp.f, tp.f) {

  ##############################################################################
  # constant of growth by using gr.cor values from 3 previous steps to fit forward projection
  # excluding division in the calculation
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #
  # Returns: - gr.cor.pred = average growth of previous time points
  #          
  ##############################################################################
  
  if (tp.f <= 8){ gr.cor.pred <<- kGR

  } else { 
  predcor <- c()
  predcor <- m.gr.cor[rp.f, (tp.f-8):tp.f]
  predcor <- predcor[!is.na(predcor)]
  predcor <- subset(predcor, (predcor <= kUncut) & (predcor >= kDiv))
  gr.cor.pred <<- round(mean(predcor), 2)
      if (is.na(gr.cor.pred) | gr.cor.pred <= 1) { 
        gr.cor.pred <<- kGR 
      } 
  }
  return(gr.cor.pred)
  
}


Plus              <- function(x) {
    
  ##############################################################################
  # Allows to deal with sum function limitation
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #
  # Returns: - fill.forward = value corrected 
  #          
  ##############################################################################
  
  
 if(all(is.na(x))){
   c(x[0], NA)} else {
   sum(x, na.rm = TRUE)}

}




# complete missing m.dat.l values in the m.cor matrix

FillForward       <- function(rp.f, tp.f, gr.cor.pred.f) {
  
  ##############################################################################
  # Check and complete empty forward time position
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #
  # Returns: - fill.forward = value corrected 
  #          
  ##############################################################################
  
  if (is.na(m.cor[rp.f, tp.f])) {
    if (tp.f < 3) {
      fill.forward <- 45
      return(fill.forward) 
    }  
    else if (tp < (nb.tp - 4)) {                                                # Correct for empty time point
      if (!is.na(Plus(m.dat.l[rp, (tp + 1):(tp + 3)])) && 
          Plus(m.dat.l[rp, (tp + 1):(tp + 3)]) > 0) {
        fill.forward <- round(m.cor[rp.f, (tp.f - 1)] * gr.cor.pred)
        return(fill.forward) 
      } else {
        fill.forward <- NA
        return(fill.forward) 
      }  

    } else {
      fill.forward <- NA
      return(fill.forward) 
    }
  }
  else if (m.cor[rp.f, tp.f] < 0){
    fill.forward <- abs(m.cor[rp.f, tp.f])
    return(fill.forward) 
  } else {
    fill.forward <- round(m.cor[rp.f, (tp.f - 1)] * gr.cor.pred)
  } 
}





#_______________________________________________________________________________



RowEmpty           <- function(rp.f, tp.f) {
  
  ##############################################################################
  # test for corrected values in matrix.cell.correct to avoid data overwritting
  # --- FOR THE REMAINING LINES (need to take into account the shift matrix) ---
  #
  # Args: - row.position.f = line analysed specific to the loop, incremented outside of the loop
  #       - time.position.f = time analysed incremented through the loop
  #
  # Returns: - row.empty value corrected with accounted cell shift
  #          
  ##############################################################################

  if (is.na(m.cor[rp.f, tp.f])) {
    row.empty <- m.dat.l[rp.f, tp.f]
    return(row.empty) 
  }
  else if (m.cor[rp.f, tp.f] > 45) {                                            # avoid over-writting of corrected values
    row.empty <- m.cor[rp.f, tp.f]
    return(row.empty)
  }
  else if (m.cor[rp.f, tp.f] < 45) {
    row.empty <- (m.dat.l[rp.f, tp.f] + m.cor[rp.f, tp.f])
    return(row.empty)
  } else {                                                                      # insert raw data
    row.empty <- m.dat.l[rp.f, tp.f]
    return(row.empty)
  }
}



#_______________________________________________________________________________
#|______________________________________________________________________________|


#_______________________________________________________________________________

# line analysis based on cell growth and rank


FirstLine          <- function(rp.f, tp.f, gr.cor.pred.f) {
  # function file
 
  
  ##############################################################################
  # test growth by using gr.raw variable for cell   
  # --- FOR THE FIRST LINE  ---
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #
  # Returns: - m.cor = matrix corrected (corrected for uncut cell, shrinkage)
  #          - m.div = matrix marking cell div, keeping trace of daughter cells
  #          - m.shift = matrix keeping track of uncut cell needed to trace cell lineage
  #
  ##############################################################################
  


  if (is.na(m.gr.raw[rp.f, tp.f])) {                                            # NA m.gr.raw
    analysis <- c(NA,                                                           # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  NA,                                                           # value for m.div
                  NA)                                                           # value for m.shift
  }
  # normal growth
  else if ((m.gr.raw[rp.f, tp.f] >= kShrink) & (m.gr.raw[rp.f, tp.f] < kUncut)) {
    if (is.na(m.cor[rp.f, tp.f])) {                                             # NA value
      analysis <- c(NA,                                                         # value for m.cor
                    NA,                                                         # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)                                                          # value for m.shift
    } 
    else if (m.cor[rp.f, tp.f] >= 600) {    # too long
      analysis <- c(600,                                                        # value for m.cor
                    NA,                                                         # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)                                                          # value for m.shift
    } else {
    analysis <- c(m.cor[rp.f, tp.f],                                            # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  0,                                                            # value for m.div
                  0)                                                            # value for m.shift
    }
  }
  # cell div  
  else if (m.gr.raw[rp.f, tp.f] <= kDiv) {                                      
    if (!is.na(m.dat.l[rp.f, (tp.f + 1)]) &&
         !is.na(m.dat.l[rp.f, (tp.f - 1)]) &&
         !is.na(m.dat.l[rp.f, tp.f]) &&
          m.dat.l[rp.f, (tp.f + 1)] <= (m.dat.l[rp.f, (tp.f - 1)] * 
                                        gr.cor.pred) * kUncut && 
        (m.dat.l[rp.f, (tp.f + 1)] > (m.dat.l[rp.f, tp.f] * kUncut))) { 
      analysis <- c(round(m.cor[rp.f, (tp.f - 1)] * gr.cor.pred),               # value for m.cor
                    NA,                                                         # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)   
    } else {
      analysis <- c(m.cor[rp.f, tp.f],                                          # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  1,                                                            # value for m.div
                  0)                                                            # value for m.shift
    }
  } 
  # uncut cell  
  else if (m.gr.raw[rp.f, tp.f] >= kUncut) {                                    
    if (!is.na(m.dat.l[rp.f, (tp.f + 1)]) &&
        !is.na(m.dat.l[rp.f, (tp.f - 1)]) &&
        !is.na(m.dat.l[rp.f, tp.f]) &&
        m.dat.l[rp.f, (tp.f + 1)] <= (m.dat.l[rp.f, (tp.f - 1)] * 
                                      gr.cor.pred) * gr.cor.pred) { 
      corr <- round(m.cor[rp.f, tp.f - 1] * gr.cor.pred)
      difference <- round(m.dat.l[rp.f, tp.f] - corr)
      if (difference < kPixelNewCell) {
        m.dat.l[(rp.f + 1), tp.f] <<- (m.dat.l[(rp.f + 1), tp.f] + difference)
        incremTooSmall <<- incremTooSmall+1
               analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred),        # value for m.cor
                             m.dat.l[(rp.f + 1), tp.f],                         # value for m.cor after cut
                             0,                                                 # value for m.div
                             0)                                                 # value for m.shift
      } else { analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred),        # value for m.cor
                             difference,                                        # value for m.cor after cut
                             0,                                                 # value for m.div
                             1)                                                 # value for m.shift 
      }  
    } else {
    corr <- round(m.cor[rp.f, tp.f - 1] * gr.cor.pred)
    difference <- round(m.dat.l[rp.f, tp.f] - corr)
      if (difference < kPixelNewCell) {m.dat.l[(rp.f + 1), tp.f] <<- 
                                (m.dat.l[(rp.f + 1), tp.f] + difference)
        incremTooSmall <<- incremTooSmall+1
        analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred),               # value for m.cor
                      m.dat.l[(rp.f + 1), tp.f],                                # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)                                                        # value for m.shift
      } else { analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred),        # value for m.cor
                             difference,                                        # value for m.cor after cut
                             0,          # extra condition if division          # value for m.div
                             1)                                                 # value for m.shift 
      }
    }
  } else {                                                                      # cells too shrunk
    GRCorPred(rp.f, tp.f)
    analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred),                   # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  0,                                                            # value for m.div
                  0)                                                            # value for m.shift
  }
  return(analysis)
}

#_______________________________________________________________________________
  

OtherLines         <- function(rp.f, tp.f, gr.cor.pred.f) {
  
  ##############################################################################
  # test growth by using gr.raw variable for cell less than variable 'kPixelForLongCell'  
  # --- FOR THE REMAINING LINES (need to take into account the shift matrix) ---
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #
  # Returns: - m.cor = matrix corrected (corrected for uncut cell, shrinkage)
  #          - m.div = matrix marking cell division, keeping trace of daughter cells
  #          - matrix.cell.shift = matrix keeping track of uncut cell needed to trace cell lineage
  ##############################################################################

if (is.na(m.cor[rp.f, tp.f])){
  analysis <- c(NA,                                                           # value for m.cor
                NA,                                                           # value for m.cor after cut
                NA,                                                           # value for m.div
                NA)                                                           # value for m.shift
}
else if (is.na(m.cor[(rp.f - 1), tp.f])){
  analysis <- c(NA,                                                           # value for m.cor
                NA,                                                           # value for m.cor after cut
                NA,                                                           # value for m.div
                NA)                                                           # value for m.shift
}
else if ( tp.f > 2 && is.na(m.cor[(rp.f - 1), (tp.f - 1)])){
  analysis <- c(NA,                                                           # value for m.cor
                NA,                                                           # value for m.cor after cut
                NA,                                                           # value for m.div
                NA)                                                           # value for m.shift
}
    
else if (m.dat.p[rp.f, tp.f] == 1) {  
  # normal growth
  if (is.na(m.gr.raw[rp.f, tp.f])) {                                            # NA value for m.gr.raw
    analysis <- c(NA,                                                           # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  NA,                                                           # value for m.div
                  NA)                                                           # value for m.shift
  }
  else if ((m.gr.raw[rp.f, tp.f] >= kShrink) & (m.gr.raw[rp.f, tp.f] < kUncut)) {   
    if (m.cor[rp.f, tp.f] >= 600) {    # too long
      analysis <- c(600,                                                        # value for m.cor
                    NA,                                                         # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)                                                          # value for m.shift
    } else {
      analysis <- c(m.cor[rp.f, tp.f],                                            # value for m.cor
                    NA,                                                           # value for m.cor after cut
                    0,                                                            # value for m.div
                    0)                                                            # value for m.shift
    }
  }
  # cell division line - 1
  else if ((m.div[rp.f - 1, tp.f] == 1) &&
           !is.na(m.div[rp.f - 1, tp.f])) {                                       
    if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) <= 
         m.cor[(rp.f - 1), (tp.f - 1)] * kUncut) {
        analysis <- c(m.cor[rp.f, tp.f],                                        # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)                                                        # value for m.shift
    }
    else if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) > 
               m.cor[(rp.f - 1), (tp.f - 1)] * kUncut)        {                 # uncut cell from dividing cells
      div <- sum(m.div[, tp.f], na.rm = T) 
      corr <- round(m.cor[rp.f, tp.f] * m.gr.raw[(rp.f - div), tp.f])
      difference <- m.cor[rp.f, tp.f] - corr
      if (difference > 45) {
        analysis <- c(corr,                                                     # value for m.cor
                      difference,                                               # value for m.cor after cut
                      0,                                                        # value for m.div
                      1)                                                        # value for m.shift
      } else {
        analysis <- c(corr,                                                       # value for m.cor
                      NA,                                                         # value for m.cor after cut
                      0,                                                          # value for m.div
                      0)                                                          # value for m.shift
      }
    }
  }
  # cell division
  else if ((m.div[rp.f - 1, tp.f] == 0) && (m.gr.raw[rp.f, tp.f] < kDiv) &&
           !is.na(m.div[rp.f - 1, tp.f])) {     
    if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) <= 
          m.cor[(rp.f - 1), (tp.f - 1)] * kUncut) {
      if (!is.na(m.dat.l[rp.f, (tp.f + 1)]) &&
          !is.na(m.dat.l[rp.f, (tp.f - 1)]) &&
          !is.na(m.dat.l[rp.f, tp.f]) &&
          (m.dat.l[rp.f, (tp.f + 1)] <= (m.dat.l[rp.f, (tp.f - 1)] * 
                                           gr.cor.pred) * kUncut) && 
            (m.dat.l[rp.f, (tp.f + 1)] > (m.dat.l[rp.f, tp.f] * 
                                            gr.cor.pred))) { 
        analysis <- c(round(m.cor[rp.f, (tp.f - 1)] * gr.cor.pred),             # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)   
      } else {
        analysis <- c(m.cor[rp.f, tp.f],                                        # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      1,                                                        # value for m.div
                      0)                                                        # value for m.shift
      }                     
    } else {
      analysis <- c(m.cor[rp.f, tp.f],                                        # value for m.cor
                    NA,                                                       # value for m.cor after cut
                    1,                                                        # value for m.div
                    0)                                                        # value for m.shift
    }
    
  }
  # uncut cell
  else if (m.gr.raw[rp.f, tp.f] >= kUncut) {                                    
    div <- sum(m.div[, tp.f], na.rm = T) 
    corr <- round(m.cor[(rp.f - div), tp.f - 1] * gr.cor.pred.f)
    difference <- m.cor[rp.f, tp.f] - corr
      if (difference > 45) {
        analysis <- c(corr,                                                     # value for m.cor
                      difference,                                               # value for m.cor after cut
                      0,                                                        # value for m.div
                      1)                                                        # value for m.shift
      } else {
        analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred.f),             # value for m.cor
                      NA,                                                        # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)                                                        # value for m.shift
          }
  }
   else {                                                                      # too shrunk cells
    analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred.f),                 # value for m.cor
                  NA,                                                            # value for m.cor after cut
                  0,                                                            # value for m.div
                  0)                                                            # value for m.shift  
  }
  return(analysis)
} else {
  # normal growth  
  if ((m.gr.raw[rp.f, tp.f] >= kShrink) & (m.gr.raw[rp.f, tp.f] < kUncut)) {    
    if (m.cor[rp.f, tp.f] >= 600) {    # too long
      analysis <- c(600,                                                        # value for m.cor
                    NA,                                                         # value for m.cor after cut
                      0,                                                          # value for m.div
                      0)                                                          # value for m.shift
    } else {
      analysis <- c(m.cor[rp.f, tp.f],                                            # value for m.cor
                    NA,                                                           # value for m.cor after cut
                    0,                                                            # value for m.div
                    0)                                                            # value for m.shift
    }
  }
  # cell division line - 1
  else if (m.div[rp.f - 1, tp.f] == 1 &&
           !is.na(m.div[rp.f - 1, tp.f]))  {                                      
    if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) <= 
          m.cor[(rp.f - 1), (tp.f - 1)] * kUncut) {
      analysis <- c(m.cor[rp.f, tp.f],                                          # value for m.cor
                    NA,                                                         # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)                                                          # value for m.shift
    }
    else if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) > 
               m.cor[(rp.f - 1), (tp.f - 1)] * kUncut)        {                 # uncut cell from dividing cells
      div <- sum(m.div[, tp.f], na.rm = T) 
      corr <- round(m.cor[rp.f, tp.f] * m.gr.raw[(rp.f - div), tp.f])
      difference <- m.cor[rp.f, tp.f] - corr
      if (difference > 45) { 
        analysis <- c(corr,                                                     # value for m.cor
                      difference,                                               # value for m.cor after cut
                      0,                                                        # value for m.div
                      1)                                                        # value for m.shift
      } else {
        analysis <- c(corr,                                                     # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)                                                        # value for m.shift
      }
    }
  }
  # cell division
  else if ((m.div[rp.f - 1, tp.f] == 0) && (m.gr.raw[rp.f, tp.f] < kDiv) &&
           !is.na(m.div[rp.f - 1, tp.f])) {    
    if ((m.cor[rp.f, tp.f] + m.cor[(rp.f - 1), tp.f]) <= 
          m.cor[(rp.f - 1), (tp.f - 1)] * kUncut) {
      if (!is.na(m.dat.l[rp.f, (tp.f + 1)]) &&
          !is.na(m.dat.l[rp.f, (tp.f - 1)]) &&
          !is.na(m.dat.l[rp.f, tp.f]) &&
          (m.dat.l[rp.f, (tp.f + 1)] <= (m.dat.l[rp.f, (tp.f - 1)] * 
                                           gr.cor.pred) * kUncut) && 
            (m.dat.l[rp.f, (tp.f + 1)] > (m.dat.l[rp.f, tp.f] * 
                                            gr.cor.pred))) { 
        analysis <- c(round(m.cor[rp.f, (tp.f - 1)] * gr.cor.pred),             # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      0,                                                        # value for m.div
                      0)   
      } else {
        analysis <- c(m.cor[rp.f, tp.f],                                        # value for m.cor
                      NA,                                                       # value for m.cor after cut
                      1,                                                        # value for m.div
                      0)                                                        # value for m.shift
      }                     
    } else {
      analysis <- c(m.cor[rp.f, tp.f],                                        # value for m.cor
                    NA,                                                       # value for m.cor after cut
                    1,                                                        # value for m.div
                    0)                                                        # value for m.shift
    }
  }
  # uncut cell
  else if (m.gr.raw[rp.f, tp.f] >= kUncut) {                                   
    div <- sum(m.div[, tp.f], na.rm = T) 
    corr <- round(m.cor[(rp.f - div), tp.f - 1] * gr.cor.pred.f)
    difference <- m.cor[rp.f, tp.f] - corr
    if (difference > 45) {
      analysis <- c(corr,                                                       # value for m.cor
                    difference,                                                 # value for m.cor after cut
                    0,                                                          # value for m.div
                    1)                                                          # value for m.shift
    } else {
      analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred.f),               # value for m.cor
                    NA,                                                          # value for m.cor after cut
                    0,                                                          # value for m.div
                    0)                                                          # value for m.shift
    }
  }
  else if (is.na(m.gr.raw[rp.f, tp.f])) {                                       
    analysis <- c(NA,                                                           # value for m.cor
                  NA,                                                           # value for m.cor after cut
                  NA,                                                           # value for m.div
                  NA)                                                           # value for m.shift
  } else {                                                                      # too shrunk cells
    analysis <- c(round(m.cor[rp.f, tp.f - 1] * gr.cor.pred.f),                 # value for m.cor
                  NA,                                                            # value for m.cor after cut
                  0,                                                            # value for m.div
                  0)                                                            # value for m.shift  
  }
  return(analysis)
}
} 


#_______________________________________________________________________________
ResultsWritting    <- function(rp.ff, tp.ff) {
  # function file
  
  
  ##############################################################################
  # transfer the data to the matices   
  # --- FOR THE FIRST LINE  ---
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #       - analysis = output of the FirstLine analysis
  #
  # Returns: - m.cor = matrix corrected (corrected for uncut cell, shrinkage)
  #          - m.div = matrix marking cell div, keeping trace of daughter cells
  #          - m.shift = matrix keeping track of uncut cell needed to trace cell lineage
  #
  ##############################################################################
  
  
  m.cor[rp.ff, tp.ff] <<- analysis[1]
  m.cor[rp.ff + 1, tp.ff] <<- analysis[2]
  m.div[rp.ff, tp.ff] <<- analysis[3]
  m.gr.raw[rp.ff, tp.ff] <<- m.gr.raw[rp.f, tp.f]
  m.shift[rp.ff, tp.ff] <<- analysis[4]
  
}

#_______________________________________________________________________________

ResultsWrittingN   <- function(rp.ff, tp.ff) {
  
  
  ##############################################################################
  # transfer the data to the matices   
  # --- FOR THE FIRST LINE  ---
  #
  # Args: - rp.f = line analysed specific to the loop, incremented outside of the loop
  #       - tp.f = time analysed incremented through the loop
  #       - analysis = output of the FirstLine analysis
  #
  # Returns: - m.cor = matrix corrected (corrected for uncut cell, shrinkage)
  #          - m.div = matrix marking cell div, keeping trace of daughter cells
  #          - m.shift = matrix keeping track of uncut cell needed to trace cell lineage
  ##############################################################################

  shift <- sum(m.shift[(rp.ff - 1), tp.ff], na.rm = T)
  if ((m.gr.raw[rp.ff, tp.ff] >= kUncut) &
      (m.cor[rp.ff, tp.ff] > 0)) {
    m.cor[((rp.ff + 2):(rp.ff + 6)), tp.ff] <<- 
      m.cor[(rp.ff + 1) :(rp.ff + 5), tp.ff]
    m.cor[rp.ff, tp.ff] <<- analysis[1]
    m.cor[(rp.ff + 1), tp.ff] <<- analysis[2]
    m.div[rp.ff, tp.ff] <<- analysis[3]
    m.shift[rp.ff, tp.ff] <<- analysis[4]
  
  } else {
    m.cor[(rp.ff + shift), tp.ff] <<- analysis[1]
    m.cor[(rp.ff + shift + 1), tp.ff] <<- analysis[2]
    m.div[rp.ff, tp.ff] <<- analysis[3]
    m.shift[rp.ff, tp.ff] <<- analysis[4]
    
  }
  
}




#################################################################################
