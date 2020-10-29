#' Exact power; point mean-difference; standard parameterization
#'
#' Calculates the exact power of an unpaired, balanced TOST with n samples in 
#' each group. The true mean-difference is assumed to be a point, i.e. fixed value.
#' @param alpha type I error
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param theta0 mean-difference
#' @param sigma standard deviation of the measurements
#' @param n sample size in each group
#'
#' @return Exact power of the TOST for a fixed point mean-difference
#' @importFrom stats qnorm pnorm uniroot integrate optimize
#' @importFrom Rdpack reprompt
#' @examples 
#' pwr <- power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
#' n.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, pwr = pwr)
#' @export
#' 
power.point <- function(alpha, theta1, theta2, theta0, sigma, n) {
  check.input.parameter(as.list(match.call())[-1])
  pwr <- PowerTOST::power.TOST(alpha = alpha, 
                               logscale = F, 
                               theta1 = theta1, 
                               theta2 = theta2, 
                               theta0 = theta0, 
                               CV = sigma, 
                               n = 2*n, 
                               design = 'parallel')
  return(pwr)
}

#' Exact power; point mean-difference; sensitivity parameterization
#' 
#' Calculates the exact power of an unpaired, balanced TOST with n samples in 
#' each group in the TOST sensitivity parameterization. The true mean-difference 
#' is assumed to be a point, i.e. fixed value.
#' 
#' @param alpha type I error
#' @param r relative mean-difference
#' @param nabla TOST sensitivity index
#' @param n sample size in each group
#'
#' @return Exact power of the TOST for a fixed point mean-difference
#' @examples 
#' power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
#' l <- standard.to.nabla(theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
#' power.point.nabla(alpha = 0.05, r = l$r, nabla = l$nabla, n = 15)
#' @export
#' 
power.point.nabla <- function(alpha, r, nabla, n) {
  check.input.parameter(as.list(match.call())[-1])
  l.standard <- nabla.to.standard(r,nabla, theta1 = -1, theta2 = 1, n)
  pwr <- power.point(alpha, theta1 = -1, theta2 = 1, 
                     l.standard$theta0, l.standard$sigma, n)
  return(pwr)
}

#' Exact sample size; point mean-difference
#' 
#' Exact sample size calculation of a balanced, unpaired TOST. The result n are 
#' the samples in one group. The true mean-difference is assumed to be a point, 
#' i.e. fixed value.
#' 
#' @param alpha type I error
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param theta0 mean-difference
#' @param sigma standard deviation of the measurements
#' @param pwr Exact power of the TOST
#' 
#' @return sample size in each group
#' @examples
#' pwr <- power.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, n = 15)
#' n.point(alpha = 0.05, theta1 = -1, theta2 = 2, theta0 = 1, sigma = 1, pwr = pwr)
#' @export
#' 
n.point <- function(alpha, theta1, theta2, theta0, sigma, pwr) {
  check.input.parameter(as.list(match.call())[-1])
  n <- PowerTOST::sampleN.TOST(alpha = alpha, 
                               targetpower = pwr, 
                               logscale = F, 
                               theta1 = theta1, 
                               theta2 = theta2, 
                               theta0 = theta0, 
                               CV = sigma, 
                               design = 'parallel', 
                               print = F)$'Sample size'
  return(n/2)
}



#' Asymptotical TOST power; point mean-difference
#' 
#' Calculates the asymptotical power of an unpaired, balanced TOST for a fixed 
#' point mean-difference. 
#' 
#' @param nabla TOST sensitivity index
#' @param r relative mean-difference
#' @param alpha type I error
#'
#' @return Asymptotic power of the TOST for a fixed point mean-difference
#' @examples 
#' power.point.nabla(alpha = 0.05, r = 2/3, nabla = 10, n = 10)
#' power.point.nabla(alpha = 0.05, r = 2/3, nabla = 10, n = 50)
#' power.point.asymp(alpha = 0.05, r = 2/3, nabla = 10)
#' @export
#' 
power.point.asymp <- function(alpha, r, nabla) {
  check.input.parameter(as.list(match.call())[-1])
  if (nabla > 2*qnorm(1-alpha)) {
    pwr <- pnorm(nabla*r - qnorm(1-alpha)) - pnorm(qnorm(1-alpha)-nabla*(1-r))
  } else {
    pwr <- 0
  }
  return(pwr)
}


#' nabla; point mean-difference
#' 
#' The equality 
#' \deqn{{power}_{asymp}(\nabla, r, \alpha) = pwr}{
#' power_asymp(nabla, r, alpha) = pwr}
#' is solved for \eqn{\nabla}{nabla} in dependence of \eqn{\alpha}{alpha}, \eqn{pwr} 
#' and \eqn{r}.
#' 
#' @param alpha type I error
#' @param pwr asymptotic power of the TOST
#' @param r relative mean-difference
#'
#' @return \eqn{\nabla}{nabla} for a fixed point mean-difference
#' @examples
#' nabla <- nabla.point(alpha = 0.05, pwr = 0.8, r = 0.55)
#' power.point.asymp(alpha = 0.05, r = 0.55, nabla = nabla)
#' @export
#' 
nabla.point <- function(alpha, pwr, r) {
  check.input.parameter(as.list(match.call())[-1])
  fct <- function(x) {
    power.point.asymp(alpha, r, x) - pwr
  }
  nabla <- uniroot2(fct, lower = 2 * qnorm(1- alpha), upper = 10^5, to.one = F)$root
  return(nabla)
}
