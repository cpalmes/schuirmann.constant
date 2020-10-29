# Calculation of Figure 1 and Figure 2. The results are written in the files
# figure_1.pdf respectively figure_2.pdf of the working directory
#
rm(list = ls())

# Setting some constants for the later calculation
nablas <- c(8,12)
n <- 15
alpha <- 0.05
rs <- c(0.2, 0.3, 0.4, 0.48)

# Constants for the subsequent graphical representation
x.axis <- c('n', 'n', 'l', 'l')
y.axis <- c('l', 'n', 'l', 'n')
x.lab <- c('', '', 'z', 'z')

# Figure 1 and Figure 2 differ from each other by a different nabla value. Thus,
# the first loop calculates Figure 1 and the second one Figure 2.
for (nabla in nablas) {
  # A pdf output format is chosen
  if (nabla == 8) {
    pdf('figure_1.pdf', height=10, width=14)
    y.range <- c(0,4)
  } else {
    pdf('figure_2.pdf', height=10, width=14)
    y.range <- c(0,6)
  }
  
  # Graphic parameters for the respective 2 x 2 grid representations
  par(mfrow = c(2,2), cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
  # Loop over each subfigure in the 2 x 2 grid
  for (k in 1:4) {
    r <- rs[k]
    
    # The notation A, B, C coincides with the one in the paper
    A <- function(z) {
      pchisq((2*n-2)*nabla^2*z^2/qt(1-alpha,2*n-2)^2,2*n-2)
    }
    B <- function(z) { dnorm(z, mean = r, sd = nabla^(-1))}
    C <- function(z) { dnorm(z, mean = 1-r, sd = nabla^(-1))}
  
    # Plot of the three curves
    curve(A, from=0, to=1, lwd = 2, lty = 4, ylim = y.range, n = 500, 
          xaxt = x.axis[k], yaxt = y.axis[k], xlab = x.lab[k], ylab = '', 
          main = paste('r = ', r))
    curve(B, from=0, to=1, lty = 1, add = T, n = 500)
    curve(C, from=0, to=1, lty = 2, add = T, n = 500)
    
    
    axis(1, at = seq(0,1,0.1))
    
    # Plot of the legend
    if (k == 3) {
      legend('topright',
             legend = c('A(z)','B(z)','C(z)'),
             lty = c(4,1,2), 
             lwd = c(2,1,1), cex = 1.3)
    }
    
    # Depiction of the two vertical lines
    abline(v = c(qnorm(1-alpha)/nabla, 0.5), col = 'grey', lwd = 2)
  }
  # close the device, i.e. close the pdf file
  dev.off()
}
