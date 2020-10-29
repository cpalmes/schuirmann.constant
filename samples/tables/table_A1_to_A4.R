# Calculation of Table A.1, A.2, A.3 and A.4. The output is written in 
# table_A1_to_A4.txt of the current working directory and contains four tables
# with weighted TOST sample sizes in LateX format using the xtable package. The 
# content of this txt-file might be copy-paste into a LaTeX-file for compilation 
# to a pdf document.
# This code might calculate a while. But on a modern standard PC it should be far 
# less than 5 minutes.
#
rm(list = ls())
library(xtable)
library(schuirmann.constant)

# Setting the alpha, power and sigma parameters for the subsequent sample size 
# calculations
alphas <- c(0.025, 0.05, 0.1)
pwrs <- c(0.8, 0.85, 0.9, 0.95)
sigmas <- seq(0.05, 0.5, 0.05)

# The sample sizes are first stored in the matrices A.zero, A.uni, A.tri and
# A.Htri with respect to their a-priori weighting.
nrow <- length(sigmas)
ncol <- length(alphas)*length(pwrs)
A.zero <- matrix(nrow = nrow, ncol = ncol)
A.uni <- matrix(nrow = nrow, ncol = ncol)
A.tri <- matrix(nrow = nrow, ncol = ncol)
A.Htri <- matrix(nrow = nrow, ncol = ncol)

# Calculation of all matrices entries by looping over all parameter combinations
idx.alpha <- 0
for (alpha in alphas) {
  idx.alpha <- idx.alpha + 1
  idx.pwr <- 0
  for (pwr in pwrs) {
    idx.pwr <- idx.pwr + 1
    idx.sigma <- 0
    for (sigma in sigmas) {
      idx.sigma <- idx.sigma + 1
      idx.row <- idx.sigma
      idx.col <- (idx.alpha - 1)*4 + idx.pwr
      
      A.zero[idx.row, idx.col] <- n.point(alpha, -0.5, 0.5, 0, sigma, pwr)
      A.uni[idx.row, idx.col] <-n.weighted(alpha, -0.5, 0.5, apriori.density('U'), 
                                           sigma, pwr) 
      A.tri[idx.row, idx.col] <- n.weighted(alpha, -0.5, 0.5, apriori.density('T'), 
                                            sigma, pwr)
      A.Htri[idx.row, idx.col] <- n.weighted(alpha, -0.5, 0.5, apriori.density('HT'), 
                                             sigma, pwr) 
    }
  }
}

# Restructure the calculated matrices to corresponding data.frames for the subsequent 
# transformations to the LaTeX tables
df.zero <- as.data.frame(A.zero)
df.uni <- as.data.frame(A.uni)
df.tri <- as.data.frame(A.tri)
df.Htri <- as.data.frame(A.Htri)
row.names(df.zero) <- sigmas
row.names(df.uni) <- sigmas
row.names(df.tri) <- sigmas
row.names(df.Htri) <- sigmas
col.names.1 <- c(rep(0.025,4), rep(0.05,4), rep(0.1,4))
col.names.2 <- rep(c(0.8,0.85,0.9,0.95),3)
col.names <- paste0(col.names.1, '_', col.names.2)
colnames(df.zero) <- col.names
colnames(df.uni) <- col.names
colnames(df.tri) <- col.names
colnames(df.Htri) <- col.names

# Write the all four data.frames to the table_A1_to_A4.txt file
capture.output( {
  chr.zero <- print(xtable(df.zero, digits=0, caption = 'zero'))
  chr.uni <- print(xtable(df.uni, digits=0, caption = 'uniform'))
  chr.tri <- print(xtable(df.tri, digits=0, caption = 'triangle'))
  chr.Htri <- print(xtable(df.Htri, digits=0, caption = 'half-triangle'))
  }, file='NUL')
fileConn<-file('table_A1_to_A4.txt')
writeLines(c(chr.zero, chr.uni, chr.tri, chr.Htri), fileConn)
close(fileConn)
