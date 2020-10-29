# Calculation of Figure 5. The result is written to the file
# figure_5.pdf in the current working directory
#
rm(list = ls())
library(latex2exp)
library(schuirmann.constant)

alpha <- 0.05
pwr <- 0.8
beta <- 1 - pwr
rs <- seq(0.5, 0.7, 0.001)

# Calculation of the Chow and Wang (2001) approximation terms
# The r > 0.5 case, cf. (13) in the paper
y.CW <- (qnorm(1-alpha) + qnorm(1-beta))/(1-rs)
# The r = 0.5 case, cf. (13)
y.CW.2 <- (qnorm(1-alpha) + qnorm(1-beta/2))/(1-rs)

# Calculation of our unbiased approximation term
nablas <- c()
for (r in rs) {
  nablas <- c(nablas, nabla.point(alpha, pwr, r))
}

# Plot of the three curves
pdf('figure_5.pdf', height=8, width=14)
par(cex.lab = 1.5, cex.axis = 1.2, cex.main = 1.5)
y.range <- range(y.CW, y.CW.2, nablas)
plot(rs, y.CW, type = 'l', lty = 4, lwd = 2, ylim = y.range, 
     ylab = 'Values of the constants', 
     main = 'Three constants in sample size approximation', xlab = 'r')
lines(rs, y.CW.2, lty = 2, lwd = 2)
lines(rs, nablas, lty = 1)

# Adding of a legend
legend('topleft',
       legend = c(TeX('$(\\phi_{1-\\alpha} + \\phi_{1-(\\beta / 2)}) / (1-r)$'), 
                  TeX('$\\nabla_{\\alpha, \\beta, r}$'), 
                  TeX('$(\\phi_{1-\\alpha} + \\phi_{1-\\beta}) / (1-r)$')),
       lty = c(2,1,4), 
       lwd = c(2,1,2), cex = 1.3)

# close the pdf file
dev.off()

