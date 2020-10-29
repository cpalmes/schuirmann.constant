# Calculation of Figure 6. The result is written to the file
# figure_6.pdf in the current working directory
#
rm(list = ls())
library(schuirmann.constant)

# The uniform (U), triangle (T) and half-triangle (HT) distribution are provided
# by the apriori.density function
density.U <- apriori.density('U')
density.T <- apriori.density('T')
density.HT <- apriori.density('HT')

# All three apriori density functions are evaluated on an equidistant grid of
# r values for the subsequent plot
rs <- seq(0,1,0.01)
dens.U <- c()
dens.T <- c()
dens.HT <- c()
for (r in rs) {
  dens.U <- c(dens.U, density.U(r))
  dens.T <- c(dens.T, density.T(r))
  dens.HT <- c(dens.HT, density.HT(r))
}
# The three apriori densities are plotted in a 1 x 3 grid representation
pdf('figure_6.pdf', height=5, width=14)
par(mfrow = c(1,3), cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2)
plot(rs, dens.U, type = 'l', ylim = c(0,4), xlab = 'r', 
     ylab = 'density', main = 'Uniform weighting density (U)')
plot(rs, dens.T, type = 'l', ylim = c(0,4), xlab = 'r', 
     ylab = '', yaxt = 'n', main = 'Triangle weighting density (T)')
plot(rs, dens.HT, type = 'l', ylim = c(0,4), xlab = 'r', 
     ylab = '', yaxt = 'n', main = 'Half-triangle weighting density (HT)')

# close the pdf file
dev.off()
