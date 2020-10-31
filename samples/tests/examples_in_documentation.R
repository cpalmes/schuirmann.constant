# Contains and runs all examples in the schuirmann.constant R-package documentation

library(schuirmann.constant)

rm(list = ls())
rho <- apriori.density('T')
r <- seq(0, 1, length = 100)
plot(r, rho(r), type = 'l')

rm(list = ls())
pwr <- power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
n.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, pwr = pwr)

rm(list = ls())
rho <- apriori.density('U')
pwr <- power.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, n = 22)
n.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, pwr = pwr)

rm(list = ls())
nabla <- nabla.point(alpha = 0.05, pwr = 0.8, r = 0.55)
power.point.asymp(alpha = 0.05, r = 0.55, nabla = nabla)

rm(list = ls())
l <- nabla.to.standard(r = 0.75, nabla = 10, theta1 = 0, theta2 = 4, n = 20)
standard.to.nabla(l$theta1,l$theta2,l$theta0, l$sigma, l$n)

rm(list = ls())
rho <- apriori.density('U')
nabla <- nabla.weighted(alpha = 0.05, pwr = 0.8, density = rho)
power.weighted.asymp(alpha = 0.05, density = rho, nabla = nabla)

rm(list = ls())
pwr <- power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
n.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, pwr = pwr)

rm(list = ls())
power.point.nabla(alpha = 0.05, r = 2/3, nabla = 10, n = 10)
power.point.nabla(alpha = 0.05, r = 2/3, nabla = 10, n = 50)
power.point.asymp(alpha = 0.05, r = 2/3, nabla = 10)

rm(list = ls())
power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
l <- standard.to.nabla(theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
power.point.nabla(alpha = 0.05, r = l$r, nabla = l$nabla, n = 15)

rm(list = ls())
rho <- apriori.density('U')
pwr <- power.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, n = 22)
n.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, pwr = pwr)

rm(list = ls())
rho <- apriori.density('HT')
power.weighted.nabla(alpha = 0.05, density = rho, nabla = 7, n = 10)
power.weighted.nabla(alpha = 0.05, density = rho, nabla = 7, n = 50)
power.weighted.asymp(alpha = 0.05, density = rho, nabla = 7)

rm(list = ls())
rho <- apriori.density('U')
power.weighted.nabla(alpha = 0.05, density = rho, nabla = 8, n = 22)
## Since we are performing a weighted power calculation, the r respectively
## theta0 value is not needed. Thus, the following calculation does not depend
## on r.
l <- nabla.to.standard(r = 0.99, nabla = 8, theta1 = -1, theta2 = 4, n = 22)
power.weighted(alpha = 0.05, theta1 = -1, theta2 = 4, density = rho, sigma = l$sigma, n = 22)

rm(list = ls())
rho <- apriori.density('U')
theta1 <- 0
theta2 <- 3
s <- schuirmann.constant(alpha = 0.05, pwr = 0.8, density = rho)
theta0 <- (1-s)*theta1 + s*theta2
n <- n.point(alpha = 0.05, theta1 = theta1, theta2 = theta2, theta0 = theta0, sigma = 1, pwr = 0.8)
## Note that the weighted power is calculated with a sample size that was found for a
## point-difference. Due to the newly introduced duality concept it is -nevertheless-
## very close to the target power 0.8.
power.weighted(alpha = 0.05, theta1 = theta1, theta2 = theta2, density = rho, sigma = 1, n = n)


rm(list = ls())
rho <- apriori.density('U')
theta1 <- 0
theta2 <- 3
s <- schuirmann.constant(alpha = 0.05, pwr = 0.8, density = rho)
(1-s)*theta1 + s*theta2
schuirmann.constant.uniform(alpha = 0.05, pwr = 0.8, theta1 = theta1, theta2 = theta2)

rm(list = ls())
l <- nabla.to.standard(r = 0.75, nabla = 10, theta1 = 0, theta2 = 4, n = 20)
standard.to.nabla(l$theta1,l$theta2,l$theta0, l$sigma, l$n)
