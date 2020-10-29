# Calculation of Table 1. The output is written in table_1.txt of the current
# working directory and contains a table of Schuirmann constants in LateX format 
# using the xtable package. The content of this txt-file might be copy-paste into 
# a LaTeX-file for compilation to a pdf document.
#
rm(list = ls())
library(xtable)
library(schuirmann.constant)

# Setting the alpha, power and densities for the subsequent Schuirmann constant
# calculations
alphas <- c(0.025, 0.05, 0.1)
pwrs <- c(0.8, 0.85, 0.9, 0.95)
densities.str <- c('U', 'T', 'HT')

# The Schuirmann constants are first stored in the matrix A.schuirmann
nrow <- length(pwrs)
ncol <- length(alphas)*length(densities.str)
A.schuirmann <- matrix(nrow = nrow, ncol = ncol)

# Calculation of the matrix entries by looping over all parameter combinations
idx.alpha <- 0
for (alpha in alphas) {
  idx.alpha <- idx.alpha + 1
  idx.pwr <- 0
  for (pwr in pwrs) {
    idx.pwr <- idx.pwr + 1
    idx.density <- 0
    for (density in densities.str) {
      idx.density <- idx.density + 1
      idx.row <- idx.pwr
      idx.col <- (idx.density - 1)*length(alphas) + idx.alpha
      A.schuirmann[idx.row, idx.col] <- 
        schuirmann.constant(alpha, pwr, apriori.density(density))
    }
  }
}

# Restructure the calculated matrix to a data.frame for subsequent transformation
# to a LaTeX table
df.schuirmann <- as.data.frame(A.schuirmann)
row.names(df.schuirmann) <- pwrs
col.names.1 <- c(rep('uniform',3), rep('triangle',3), rep('half-triangle',3))
col.names.2 <- rep(alphas,3)
col.names <- paste0(col.names.1, '_', col.names.2)
colnames(df.schuirmann) <- col.names

# Write the df.schuirmann data.frame to the table_1.txt file
capture.output( {
  chr.schuirmann <- print(xtable(df.schuirmann, digits = 3, 
                                 caption = 'Schuirmann-Constants'))
  }, file='NUL')
fileConn <- file('table_1.txt')
writeLines(chr.schuirmann, fileConn)
close(fileConn)
