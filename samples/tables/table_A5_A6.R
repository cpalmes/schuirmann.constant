# Calculation of Table A.5 and A.6. The output is written in table_A5_A6.txt of 
# the current working directory and contains two tables of nabla constants (fixed
# mean-difference and weighted mean-difference) in LateX format using the xtable 
# package. The content of this txt-file might be copy-paste into a LaTeX-file for 
# compilation to a pdf document.
#
rm(list = ls())
library(xtable)
library(schuirmann.constant)

# Parameter settings for the subsequent nabla calculations
is.weighted <- c('F', 'T')
alphas <- c(0.025, 0.05, 0.1)
pwrs <- c(0.8, 0.85, 0.9, 0.95)
densities.str <- c('U', 'T', 'HT')
rs <- c(0.5, 0.55, 0.6)

# The nabla constants are first stored in the matrices A.nabla.point and 
# A.nabla.weighted 
nrow <- length(pwrs)
ncol.point <- length(alphas)*length(rs)
ncol.weighted <- length(alphas)*length(densities.str)
A.nabla.point <- matrix(nrow = nrow, ncol = ncol.point)
A.nabla.weighted <- matrix(nrow = nrow, ncol = ncol.weighted)

# Calculation of all matrices entries by looping over all parameter combinations
idx.alpha <- 0
for (alpha in alphas) {
  idx.alpha <- idx.alpha + 1
  idx.pwr <- 0
  for (pwr in pwrs) {
    idx.pwr <- idx.pwr + 1
    idx.density <- 0
    idx.r <- 0
    for (weight.it in is.weighted) {
      if (weight.it) {
        for (density in densities.str) {
          idx.density <- idx.density + 1
          idx.row <- idx.pwr
          idx.col <- (idx.density - 1)*length(alphas) + idx.alpha
          A.nabla.weighted[idx.row, idx.col] <- 
            nabla.weighted(alpha, pwr, apriori.density(density))
        }    
      } else {
        for (r in rs) {
          idx.r <- idx.r + 1
          idx.row <- idx.pwr
          idx.col <- (idx.r - 1)*length(alphas) + idx.alpha
          A.nabla.point[idx.row, idx.col] <- nabla.point(alpha, pwr, r)
        }
      }
    }  #end of is.weighted
  }
}

# Restructure the calculated matrices to corresponding data.frames for the subsequent 
# transformations to the LaTeX tables
df.nabla.point <- as.data.frame(A.nabla.point)
df.nabla.weighted <- as.data.frame(A.nabla.weighted)
row.names(df.nabla.point) <- pwrs
row.names(df.nabla.weighted) <- pwrs
col.names.1.point <- c(rep(rs[1],3), rep(rs[2],3), rep(rs[3],3))
col.names.1.weighted <- c(rep('uniform',3), rep('triangle',3), 
                          rep('half-triangle',3))
col.names.2 <- rep(alphas,3)
col.names.point <- paste0(col.names.1.point, '_', col.names.2)
col.names.weighted <- paste0(col.names.1.weighted, '_', col.names.2)
colnames(df.nabla.point) <- col.names.point
colnames(df.nabla.weighted) <- col.names.weighted

# Write the both data.frames to the table_A5_to_A6.txt file
capture.output( {
  chr.nabla.point <- print(xtable(df.nabla.point, digits=3, 
                                  caption = 'Nabla-point'))
  chr.nabla.weighted <- print(xtable(df.nabla.weighted, digits=3, 
                                     caption = 'Nabla-weighted'))
  }, file='NUL')

fileConn<-file('table_A5_A6.txt')
writeLines(c(chr.nabla.point, chr.nabla.weighted), fileConn)
close(fileConn)
