# Calculation of Figure 4. The result is written to the file
# figure_4.pdf in the current working directory
#
rm(list = ls())
library(latex2exp)
library(schuirmann.constant)

# Our unbiased approximation formula of the TOST sample size
approx.n.corrected <- function(sigma, theta, alpha, pwr, r) {
  nabla <- nabla.point(alpha, pwr, r)
  n <- (sigma/theta)^2 * nabla^2/2
  n <- round(max(n,2))
}

# Implementation of the TOST sample size approximation by Chow and Wang (2001)
approx.n.literature.CW <- function(sigma, theta, alpha, pwr, r) {
  beta <- 1 - pwr
  # Calculation of the right-hand side of the formula (12) in our paper
  get.rs <- function(n) {
    if (r == 0.5) {
      x <- (qt(1-alpha, df = 2*n-2) + qt(1-beta/2, df = 2*n-2))/(1-r)
    } else {
      x <- (qt(1-alpha, df = 2*n-2) + qt(1-beta, df = 2*n-2))/(1-r)
    }
    n.rs <- (sigma/theta)^2 * x^2/2
  }
  
  # Linear search for the smallest n that fulfills the condition in formula (12)
  n <- 2
  while (n < get.rs(n)) {
    n <- n + 1
  }
    
  return(n)    
}

# The calculated parameter setting is saved by the subsequent constants
alpha <- 0.05
theta <- 0.5
sigmas <- seq(0.01, 1, 0.01)
pwr <- 0.8
r <- 0.55
theta.0 <- (1-r)*(-theta) + r*theta

n.corrected <- c()
n.literature.CW <- c()
pwr.corrected <- c()
pwr.literature.CW <- c()

# The approximated sample size (left-hand plot) and TOST power (right-hand plot)
# are calculated via a loop over sigma, i.e. over the x-axis
for (sigma in sigmas) {
  # Calculation of the approximated sample size
  n.cor <- approx.n.corrected(sigma, theta, alpha, pwr, r)
  n.corrected <- c(n.corrected, n.cor)
  n.lit.CW <- approx.n.literature.CW(sigma, theta, alpha, pwr, r)
  n.literature.CW <- c(n.literature.CW, n.lit.CW)
  
  # Calculation of the TOST power approximation
  pwr.cor <- power.point(alpha, -theta, theta, theta.0, sigma, n.cor)
  pwr.corrected <- c(pwr.corrected, pwr.cor)
  pwr.lit.CW <- power.point(alpha, -theta, theta, theta.0, sigma, n.lit.CW)
  pwr.literature.CW <- c(pwr.literature.CW, pwr.lit.CW)
}

# Plot of the approximated sample size (left-hand side) of Figure 4
pdf('figure_4.pdf', height=7, width=14)
par(mfrow = c(1,2), cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
y.range <- c(0, max(n.literature.CW, n.corrected))
plot(x = sigmas, y = n.corrected, ylim = y.range, 
     xlab = TeX('$\\sigma / (\\theta_2 - \\theta_1)$'), pch = 16, ylab = 'n', 
     main = 'Sample size approximation')
points(x = sigmas, y = n.literature.CW, pch = 1)

# Plot of the approximated TOST power (right-hand side) of Figure 4
plot(x = sigmas, y = pwr.corrected, ylim = c(0,1), 
     xlab = TeX('$\\sigma / (\\theta_2 - \\theta_1)$'), pch = 16, 
     ylab = 'power', main = 'Exact TOST power')
abline(h = pwr, lwd = 2, lty = 2)
points(x = sigmas, y = pwr.literature.CW, ylim = c(0,1), pch = 1)
legend('topright',
       legend = c('our improvement',TeX('literature: $n_{CW}$')),
       pch = c(16,1), cex = 1.3)

# close the pdf file
dev.off()
