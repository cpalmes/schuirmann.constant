# Calculation of Figure 7. The result is written to the file
# figure_6.pdf in the current working directory
#
rm(list = ls())
library(latex2exp)
library(schuirmann.constant)

alpha <- 0.05
pwr <- 0.8

# Calculation of the Schuirmann constant for the triangle a-priori distribution
s <- schuirmann.constant(alpha, pwr, apriori.density('T'))

theta.1 <- -0.5
theta.2 <- 0.5
theta.0 <- (1-s)*theta.1 + s*theta.2
sigmas <- seq(0.01,0.5,0.005)

# Calculate the sample size using the Schuirmann constant as a fixed mean-difference
ns <- c()
for (sigma in sigmas) {
  n <- n.point(alpha, theta.1, theta.2, theta.0, sigma, pwr)
  ns <- c(ns, n)
}

# Calculate the (exact) power with the triangle density as a-priori distribution
# for the true mean-difference and the former calculated sample sizes using the
# Schuirmann constant as fixed mean-difference
pwrs <- c()
for (k in 1:length(sigmas)) {
  n <- ns[k]
  pwr <- power.weighted(alpha, theta.1, theta.2, apriori.density('T'), sigmas[k], n) 
  pwrs <- c(pwrs, pwr)
}


# Plot the sample sizes on the left-hand side and the power calculations on the
# right-hand side within a 1 x 2 grid presentation
pdf('figure_7.pdf', height=7, width=14)
par(mfrow = c(1,2), cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
plot(sigmas, ns, xlab = TeX('$\\sigma / (\\theta_2 - \\theta_1)$'), ylab = 'n',
     main = 'Sample size calculation with mean-difference S', pch = 16)
plot(sigmas, pwrs, ylim = c(0,1), 
     xlab = TeX('$\\sigma / (\\theta_2 - \\theta_1)$'), ylab = 'power',
     main = 'Triangle weighted power', pch = 16)
abline(h = 0.8)

# close the pdf file
dev.off()

