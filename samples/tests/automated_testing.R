# Automated testing functions for the schuirmann.constant R-package
rm(list = ls())
library('schuirmann.constant')

# Automated tests for the functions nabla.point, nabla.weighted and schuirmann.constant. 
# The Schuirmann constant is calculated for a randomly chosen beta distribution.
# Then, the accuracy of its defining property nabla.point = nabla.weight is checked
# together with the nabla defining properties power.asym(nabla) = pwr.
# alpha resp. pwr are drawn randomly from the interval ]0, 0.5[ respectively ]0.5, 1[
# and all calculations are repeated in a loop M = 10^3 times.
#
check.schuirmann.nabla <- function(M = 10^3) {
  pwrs <- c()
  nabla.ps <- c()
  nabla.ws <- c()
  pwr.points <- c()
  pwr.weights <- c()
  ss <- c()
  delta.nabla.rels <- c()
  delta.pwrs <- c()
  # Repeat all simulations M times
  for (m in 1:M) {
    # Simulate alpha, pwr and pick randomly a beta function
    alpha <- runif(n = 1, min = 0, max = 0.5)
    pwr <- runif(n = 1, min = 0.5, max = 1)
    pwrs <- c(pwrs, pwr)
    shapes <- runif(n = 2, min = 0.5, max = 5)
    rho <- function(x) { dbeta(x, shape1 = shapes[1], shape2 = shapes[2]) }
    
    # Calculate the Schuirmann constant
    s <- schuirmann.constant(alpha, pwr, density = rho)
    ss <- c(ss, s)
    # Calculate also the corresponding point nabla and the weighted nabla as well
    # in order to check for their equality since this is the defining property of the
    # Schuirmann constant.
    nabla.p <- nabla.point(alpha, pwr, r = s)
    nabla.w <- nabla.weighted(alpha, pwr, density = rho)
    nabla.ps <- c(nabla.ps, nabla.p)
    nabla.ws <- c(nabla.ws, nabla.w)
    
    # Calculate also both asymptotic powers, which should be close to the target
    # power pwr.
    pwr.point <- power.point.asymp(alpha, r = s, nabla.p)
    pwr.weight <- power.weighted.asymp(alpha, density = rho, nabla.p)
    pwr.points <- c(pwr.points, pwr.point)
    pwr.weights <- c(pwr.weights, pwr.weight)
    
    # Calculate the relative nabla differences and, more importantly, the absolute
    # asymptotic power differences.
    delta.nabla.rel <- (nabla.w - nabla.p) / (nabla.w + nabla.p)
    delta.nabla.rels <- c(delta.nabla.rels, delta.nabla.rel)
    delta.pwr <- pwr.weight - pwr.point
    delta.pwrs <- c(delta.pwrs, delta.pwr)
  }
  
  # Summarize all results in a common data.frame and return it.
  df <- data.frame(pwr = pwrs, 
                   pwr.point = pwr.points,
                   pwr.weight = pwr.weights,
                   nabla.p = nabla.ps,
                   nabla.w = nabla.ws,
                   s = ss,
                   delta.pwr = delta.pwrs,
                   delta.nabla.rel = delta.nabla.rels)
  return(df)
}


# Automated test for the functions n.weighted and power.weighted.
# For randomly chosen beta distributions and 0 < alpha < 0.5, 0.5 < pwr < 1, 
# -10 < theta1 < theta2 < 10 and 0 < sigma < 20, the minimum sample size n for 
# exceeding the given power pwr is calculated using the n.weighted function. 
# For verification, the function power.weighted is applied afterwards. All calculations 
# are repeated in a loop M times.
#
check.n.power.weighted <- function(M = 10^2) {
  ns <- c()
  pwrs <- c()
  pwr.lows <- c()
  pwr.highs <- c()
  # Repeat all simulations M times
  for (m in 1:M) {
    # Print the actual simulation step for awareness 
    print(m)
    # Simulate alpha, pwr, theta1, theta2 and rho
    alpha <- runif(n = 1, min = 0, max = 0.5)
    pwr <- runif(n = 1, min = 0.5, max = 1)
    pwrs <- c(pwrs, pwr)
    thetas <- runif(n = 2, min = -10, max = 10)
    theta1 <- min(thetas)
    theta2 <- max(thetas)
    sigma <- runif(n = 1, min = 0, max = 20)
    shapes <- runif(n = 2, min = 0.5, max = 5)
    rho <- function(x) { dbeta(x, shape1 = shapes[1], shape2 = shapes[2]) }
    
    # Calculate the minimum sample size for exceeding the target power
    n <- n.weighted(alpha, theta1, theta2, density = rho, sigma, pwr)
    ns <- c(ns, n)
    pwr.high <- power.weighted(alpha, theta1, theta2, density = rho, sigma, n)
    if (n > 2) {
      pwr.low <- power.weighted(alpha, theta1, theta2, density = rho, sigma, n-1)
    } else {
      pwr.low <- 0
    }
    
    # Save the established power as well as the power for the reduced sample size
    # n-1 since this has to be strictly lower than the given target power.
    pwr.lows <- c(pwr.lows, pwr.low)
    pwr.highs <- c(pwr.highs, pwr.high)
  }
  
  # Summarize all results in a common data.frame and return it.
  df <- data.frame(n = ns, 
                   pwr.low = pwr.lows, 
                   pwr = pwrs, 
                   pwr.high = pwr.highs)
  
  # The function pair (n.weighted, power.weighted) works as expected if and only
  # if the following condition is true. A later summation over df$passed will be
  # calculated for verification.
  df$passed <- (pwr.high >= pwr) & (pwr.low < pwr)

  return(df)
}


# Automated testing for the function power.point.asymp.
# The goodness of approximation of the power.point.asymp function is checked for
# the sample size n = 20. For this, 0 < alpha < 0.5, 0.5 < r < 1 and  0 < nabla < 20
# are simulated and the exact point power in this setting is calculated for n = 20.
# Next, the asymptotic point power is calculated in the same parameter setting.
# Finally, the absolute difference between both powers is considered. All calculations 
# are repeated in a loop M = 10^3 times. 
#
check.power.point.asymp <- function(M = 10^3) {
  pwr.exacts <- c()
  pwr.asymps <- c()
  
  # Repeat all simulations M times
  for (m in 1:M) {
    # Simulate alpha, nabla and r
    alpha <- runif(n = 1, min = 0, max = 0.5)
    nabla <- runif(n = 1, min = 0, max = 20)
    r <- runif(n = 1, min = 0.5, max = 1)
    
    # Calculate the exact and asymptotic power
    pwr.exact <- power.point.nabla(alpha, r, nabla, n = 20)
    pwr.asymp <- power.point.asymp(alpha, r, nabla)
    pwr.exacts <- c(pwr.exacts, pwr.exact)
    pwr.asymps <- c(pwr.asymps, pwr.asymp)
  }
  
  # Store the exact and asymptotic power in a data.frame together with their 
  # difference
  df <- data.frame(pwr.exact = pwr.exacts, pwr.asymp = pwr.asymps)
  df$diff <- df$pwr.asymp - df$pwr.exact
  return(df)
}


# Automated testing for the function power.weighted.asymp
# The goodness of approximation of the power.weighted.asymp function is checked for
# the sample size n = 20. For this, 0 < alpha < 0.5, 0 < nabla < 20 and a beta distribution
# for the true mean-difference are simulated and the exact point power in this setting 
# is calculated for n = 20. Next, the asymptotic point power is calculated in the 
# same parameter setting. Finally, the absolute difference between both powers is 
# considered. All calculations are repeated in a loop M = 10^3 times. 
#
check.power.weighted.asymp <- function(M = 10^3) {
  pwr.exacts <- c()
  pwr.asymps <- c()
  
  # Repeat all simulations M times
  for (m in 1:M) {
    # Simulate alpha, nabla and beta
    alpha <- runif(n = 1, min = 0, max = 0.5)
    nabla <- runif(n = 1, min = 0, max = 20)
    shapes <- runif(n = 2, min = 0.5, max = 5)
    rho <- function(x) { dbeta(x, shape1 = shapes[1], shape2 = shapes[2]) }
    
    # Calculate the exact and asymptotic weighted power
    pwr.exact <- power.weighted.nabla(alpha, rho, nabla, n = 20)
    pwr.asymp <- power.weighted.asymp(alpha, rho, nabla)
    pwr.exacts <- c(pwr.exacts, pwr.exact)
    pwr.asymps <- c(pwr.asymps, pwr.asymp)
  }
  
  # Store the exact and asymptotic weighted power in a data.frame together with 
  # their difference
  df <- data.frame(pwr.exact = pwr.exacts, pwr.asymp = pwr.asymps)
  df$diff <- df$pwr.asymp - df$pwr.exact
  
  return(df)
}

###############################################################################

# Set a seed for reproducibility 
set.seed(123)
# Run automated tests for the  nabla.point nabla.weighted and schuirmann.constant
# functions and plot relevant differences. See the comments in the check.schuirmann.nabla
# function for further information.
df <- check.schuirmann.nabla(M = 10^3)
plot(df$delta.nabla.rel)
plot(df$delta.pwr)

# Run automated test for the n.weighted and power.weighted function. See the comments
# in the check.n.power.weighted function for further information.
df <- check.n.power.weighted(M = 10^2)
sum(df$passed) # Should be equal to M

# Automated testing for the function power.point.asymp, cf. the comments in the 
# function check.power.point.asymp for further information.
df <- check.power.point.asymp(M = 10^3)
plot(df$diff)

# Automated testing for the function power.weighted.asymp, cf. the comments in the 
# function check.power.weighted.asymp for further information.
df <- check.power.weighted.asymp(M = 10^3)
plot(df$diff)
