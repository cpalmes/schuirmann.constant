# Calculation of Figure 3. The result is written to the file
# figure_3.pdf in the current working directory
#
rm(list = ls())
# Our new schuirmann.constant R-package 
library(schuirmann.constant)
# Package for latex representations in the some of the figure labels
library(latex2exp)

# Definition of some constants for convenience 
alpha <- 0.05
nablas <- c(8, 12)
y.labels <- c('power', '')
y.axis <- c('l', 'n')
tex.nabla <- c(TeX('$\\nabla = 8$'), TeX('$\\nabla = 12$'))

# the output file consists of a 1 x 2 grid. The left-hand side depicts the
# nabla = 8 case and the right-hand side the lambda = 12 case.
pdf('figure_3.pdf', height=7, width=14)
par(mfrow = c(1,2), cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
# Each of the two loops generates on of the two grids
for (k in 1:2) {
  nabla <- nablas[k]
  
  # Asymptotic power calculation
  rs <- seq(0.5,1,0.01)
  pwr.approx <- c()
  for (r in rs) {
    pwr <- power.point.asymp(alpha, r, nabla)
    pwr.approx <- c(pwr.approx, pwr)
  }
  plot(x = rs, y = pwr.approx, ylim = c(0,1), type = 'l', lty = 2, lwd = 2, 
       xlab = 'r', ylab = y.labels[k], yaxt = y.axis[k], main = tex.nabla[k])
  
  # Exact power calculation
  for (n in c(2:5, 10, 20)) {
    pwr.n <- c()
    for (r in rs) {
      pwr <- power.point.nabla(alpha, r, nabla, n)
      pwr.n <- c(pwr.n, pwr)
    }
    lines(x = rs, y = pwr.n)
    
  }
  
  # Insert a legend in the left grid, i.e. for k = 1
  if (k == 1) {
    legend('topright',
           legend = c(TeX('$power_{asymp}$'),
                      TeX('$power\\,(n = 2,3,4,5, 10, 20)$')),
           lty = c(2,1), 
           lwd = c(2,1), cex = 1.3)
  }
}
# close the pdf file
dev.off()

