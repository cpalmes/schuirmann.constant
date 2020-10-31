#' Returns a standard a-priori density function
#' 
#' Let U denote the uniform distribution on \eqn{[0,1]}, T the triangle distribution 
#' with center point \eqn{0.5} and full support on \eqn{[0,1]} and HT the half-triangle distribution
#' with center point \eqn{0.5} and support on \eqn{[0.25, 0.75]}. This function returns one of
#' these a-priori densities in dependence of its input parameter. These three
#' a-priori distributions are deemed as helpful in experimenting with
#' the Schuirmann constant. 
#' 
#' @param density.id 'U', 'T' or 'HT', depending on which a-priori density should
#' be returned
#' @return Uniform (U), triangle (T) or half triangle (HT) a-priori density function
#' @examples  
#' rho <- apriori.density('T')
#' r <- seq(0, 1, length = 100)
#' plot(r, rho(r), type = 'l')
#' @export
#' 
apriori.density <- function(density.id) {
  # Uniform density
  density.U <- function(r) {
    rho.U <- ((r >= 0) & (r <= 1))*1.
    return(rho.U)
  }
  # Triangle density
  density.T <- function(r) {
    rho.T <- ((r >= 0) & (r <= 1))*(2-4*abs(r-0.5))
    return(rho.T)
  }
  # Half-triangle density
  density.HT <- function(r) {
    rho.HT <- ((r >= 0.25) & (r <= 0.75))*(4-16*abs(r-0.5))
    return(rho.HT)
  }  
  
  density.id <- tolower(density.id)
  if (density.id == 'u') {
    return (density.U)
  } else if (density.id == 't') {
    return (density.T)
  } else if (density.id == 'ht') {
    return (density.HT)
  } else {
    stop('Invalid input! density.id must be \'U\', \'T\' or \'HT\'')
  }
}


#' Exact weighted power; standard parameterization
#' 
#' Calculates the exact power of an unpaired, balanced TOST with 
#' n samples in each group, where the true mean-difference is modeled with an
#' a-priori density rather than a fixed mean-difference.
#' 
#' @param alpha type I error
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param density a-priori density for the mean-difference
#' @param sigma standard deviation of the measurements
#' @param n sample size in each group
#'
#' @return Exact power of the TOST for a weighted true mean-difference
#' @examples
#' rho <- apriori.density('U')
#' pwr <- power.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, n = 22)
#' n.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, pwr = pwr)
#' @export
#' 
power.weighted <- function(alpha, theta1, theta2, density, sigma, n) {
  check.input.parameter(as.list(match.call())[-1])
  integrand <- function(r) {
    theta0 <- (1-r)*theta1 + r*theta2
    power.point(alpha, theta1, theta2, theta0, sigma, n) * density(r)
  }
  pwr <- integrate(integrand, lower = 0, upper = 1)$value
  return(pwr)
}



#' Exact weighted power; sensitivity parameterization
#' 
#' Calculates the exact power of an unpaired, balanced TOST with 
#' n samples in each group, where the true mean-difference is modeled with an
#' a-priori density rather than a fixed mean-difference. The calculation is performed 
#' in the TOST sensitivity parameterization.
#' 
#' @param alpha type I error
#' @param density a-priori density for the mean-difference
#' @param nabla TOST sensitivity index
#' @param n sample size in each group
#'
#' @return Exact power of the TOST for a weighted mean-difference
#' @examples
#' rho <- apriori.density('U')
#' power.weighted.nabla(alpha = 0.05, density = rho, nabla = 8, n = 22)
#' ## Since we are performing a weighted power calculation, the r respectively 
#' ## theta0 value is not needed. Thus, the following calculation does not depend
#' ## on r.
#' l <- nabla.to.standard(r = 0.99, nabla = 8, theta1 = -1, theta2 = 4, n = 22)
#' power.weighted(alpha = 0.05, theta1 = -1, theta2 = 4, density = rho, sigma = l$sigma, n = 22)
#' @export
#' 
power.weighted.nabla <- function(alpha, density, nabla, n) {
  check.input.parameter(as.list(match.call())[-1])
  integrand <- function(r) {
    power.point.nabla(alpha, r, nabla, n) * density(r)
  }
  pwr <- integrate(integrand, lower = 0, upper = 1)$value
  return(pwr)
}


#' Exact sample size; weighted mean-difference
#' 
#' Exact sample size calculation of a balanced, unpaired TOST, where the true
#' mean-difference follows a prescribed a-priori distribution. The lowest sample 
#' size n that exceeds the target power pwd is calculated using a standard binary 
#' search algorithm. The result n is the sample size in one group.
#'
#' @param alpha type I error
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param density a-priori density for the mean-difference
#' @param sigma standard deviation of the measurements
#' @param pwr Exact power of the TOST
#'
#' @return sample size in each group
#' @examples
#' rho <- apriori.density('U')
#' pwr <- power.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, n = 22)
#' n.weighted(alpha = 0.05, theta1 = -2, theta2 = 2, density = rho, sigma = 1, pwr = pwr)
#' @export
#' 
n.weighted <- function(alpha, theta1, theta2, density, sigma, pwr) {
  check.input.parameter(as.list(match.call())[-1])
  
  n.low <- 2
  n.high <- 64
  
  # Check whether two samples are enough. If not, then we have
  # power(n.low) < pwr
  if (power.weighted(alpha, theta1, theta2, density, sigma, n.low) >= pwr) {
    return(n.low)
  }
  
  # Find n.high, such that
  # power(n.low) < pwr <= power(n.high)
  pwr.high <- power.weighted(alpha, theta1, theta2, density, sigma, n.high)
  while (pwr.high < pwr) {
    n.low <- n.high
    n.high <- 2*n.high
    pwr.high <- power.weighted(alpha, theta1, theta2, density, sigma, n.high)
  }
  
  # Repeat until n.high = n.low + 1 and
  # power(n.heigh-1) < pwr <= power(n.heigh)
  while (n.high - n.low > 1) {
    n.inter <- floor((n.low + n.high)/2)
    pwr.inter <- power.weighted(alpha, theta1, theta2, density, sigma, n.inter)
    if (pwr.inter < pwr) {
      n.low <- n.inter  
    } else {
      n.high <- n.inter
    }
  }
  
  return(n.high)
}


#' Asymptotical weighted power
#' 
#' Calculates the asymptotical power of an unpaired, balanced TOST 
#' for a weighted a-priori mean-difference. 
#' 
#' @param alpha type I error
#' @param density a-priori density for the mean-difference
#' @param nabla TOST sensitivity index
#'
#' @return Asymptotic power of the TOST for a weighted mean-difference
#' @examples 
#' rho <- apriori.density('HT')
#' power.weighted.nabla(alpha = 0.05, density = rho, nabla = 7, n = 10)
#' power.weighted.nabla(alpha = 0.05, density = rho, nabla = 7, n = 50)
#' power.weighted.asymp(alpha = 0.05, density = rho, nabla = 7)
#' @export
#' 
power.weighted.asymp <- function(alpha, density, nabla) {
  check.input.parameter(as.list(match.call())[-1])
  integrand <- function(r) { 
    (pnorm(nabla*r - qnorm(1-alpha)) - pnorm(qnorm(1-alpha)-nabla*(1-r))) * density(r)  
  }
  if (nabla > 2*qnorm(1-alpha)) {
    pwr <- integrate(integrand, lower = 0, upper = 1)$value
  } else {
    pwr <- 0
  }
  return(pwr)
}

#' nabla; weighted mean-difference
#'
#' The equality 
#' \deqn{power_{asymp, density}(\nabla, \alpha) = pwr}{
#' power_asymp_density(nabla, alpha) = pwr}
#' is solved for \eqn{\nabla}{nabla} in dependence of \eqn{\alpha}{alpha}, 
#' \eqn{pwr}{pwr} and \eqn{density}{density}.
#' 
#' @param alpha type I error
#' @param pwr asymptotic weighted power of the TOST
#' @param density a-priori density for the mean-difference
#'
#' @return \eqn{\nabla}{nabla} for a weighted mean-difference
#' @examples
#' rho <- apriori.density('U')
#' nabla <- nabla.weighted(alpha = 0.05, pwr = 0.8, density = rho)
#' power.weighted.asymp(alpha = 0.05, density = rho, nabla = nabla)
#' @export
#' 
nabla.weighted <- function(alpha, pwr, density) {
  check.input.parameter(as.list(match.call())[-1])
  fct <- function(x) {
    power.weighted.asymp(alpha, density, x) - pwr
  }
  nabla <- uniroot2(fct, lower = 2*qnorm(1-alpha), upper = 10^5, to.one = F)$root
  return(nabla)
}
