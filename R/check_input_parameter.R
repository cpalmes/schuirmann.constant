# Description of intern package function check.input.parameter:
# This function is the first call in all schuirmann.const functions that are
# available to the package user. It checks for the validity of each parameter delivered
# by the input list l.para. If any parameter in l.para has an invalid value, the
# function ends with an informative error message.
#
check.input.parameter <- function(l.para) {
  
  # For performance reasons a parameter check should be done only in case a function
  # is invoked externally by the package-user. If a function is invoked internally, e.g.
  # as part of a calculation algorithm, no parameter check is needed.
  # 
  if (!identical(parent.frame(n = 2), globalenv())) {
    return()
  } 
  
  # Checks for a valid apriori density rho. Three properties are checked:
  # (i) rho has to be a function, i.e. class(rho) == 'function'
  # (ii) rho(r) >= 0 on the interval [0, 1]
  # (iii) rho integrates to one on the interval [0, 1]
  #
  check.density <- function(density, eps = 10^(-6)) {
    if (class(density) != 'function') {
      stop('The a-priori density has to be a function on [0, 1]!')
    } else {
      density.min <- optimize(density, interval = c(0, 1), maximum = F)
      if (density.min$objective < 0) {
        stop('The a-priori density has to be non-negative!')
      } else {
        density.int <- integrate(density, lower = 0, upper = 1)$value
        if (abs(density.int - 1) > eps) {
          stop('The a-priori density needs to integrate to 1 over [0, 1]!')
        }
      }
    }
  }
  
  # Checks for 0 < alpha < 0.5
  check.alpha <- function(alpha) {
    if (alpha <= 0 || alpha >= 0.5) {
      stop('alpha needs to be in the open interval ]0, 0.5[!')
    }
  }
  
  # Checks for nabla > 0
  check.nabla <- function(nabla) {
    if (nabla <= 0) {
      stop('nabla > 0 is not fulfilled!')
    }
  }
  
  # Checks for 0 < pwr < 1
  check.pwr <- function(pwr) {
    if (pwr <= 0.5 || pwr >= 1) {
      stop('pwr needs to be in the open interval ]0.5, 1[!')
    } 
  }
  
  # Checks for :
  # (i) n needs to be an integer
  # (ii) n >= 2
  check.n <- function(n) {
    if (abs(round(n) - n) != 0 || n <= 1) {
      stop('n needs to be an integer >= 2!')
    }
  }
  
  # Checks for theta1 < theta2 and, if theta0 is available additionally for
  # theta0 in ]theta1, theta2[
  check.theta <- function(theta1, theta2, theta0 = NA) {
    if (theta1 >= theta2) {
      stop('theta1 < theta2 is not fulfilled!')
    }
    if (!is.na(theta0)) {
      if (theta0 <= theta1 || theta0 >= theta2) {
        stop('theta1 < theta0 < theta2 is not fulfilled!')
      }
    } 
  }
  
  # Checks for sigma > 0
  check.sigma <- function(sigma) {
    if (sigma <= 0) {
      stop('sigma > 0 is not fulfilled!')
    }
  }
  
  # Checks for 0.5 <= r <= 1
  check.r <- function(r) {
    if (r < 0.5 || r > 1) {
      stop('0.5 <= r <= 1 is not fulfilled!')
    }
  }

  # check all parameters in l.para for validity by invoking their respective
  # check function as defined above.
  l.para <- lapply(l.para, eval)
  if ('density' %in% names(l.para)) {
    check.density(l.para$density)
  }
  
  if ('alpha' %in% names(l.para)) {
    check.alpha(l.para$alpha)
  }
  
  if ('pwr' %in% names(l.para)) {
    check.pwr(l.para$pwr)
  }
  
  if ('n' %in% names(l.para)) {
    check.n(l.para$n)
  }
  
  if ('theta1' %in% names(l.para)) {
    if ('theta0' %in% names(l.para)) {
      check.theta(l.para$theta1, l.para$theta2, l.para$theta0)
    } else {
      check.theta(l.para$theta1, l.para$theta2)
    }
  }
  
  if ('sigma' %in% names(l.para)) {
    check.sigma(l.para$sigma)
  }
  
  if ('r' %in% names(l.para)) {
    check.r(l.para$r)
  }
  
  if ('nabla' %in% names(l.para)) {
    check.nabla(l.para$nabla)
  }
  
}

