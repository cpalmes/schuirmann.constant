#' TOST sensitivity to standard parameterization
#' 
#' The TOST sensitivity parameterization is transformed to the standard parameterization.
#' @param r relative mean-difference
#' @param nabla TOST sensitivity index
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param n sample size in each group
#' 
#' @return A list with the transformed parameters in standard parameterization, i.e. 
#' \eqn{(\theta_1, \theta_2, \theta_0, \sigma, n)}{(theta1, theta2, theta0, sigma, n)}
#' @examples
#' l <- nabla.to.standard(r = 0.75, nabla = 10, theta1 = 0, theta2 = 4, n = 20)
#' standard.to.nabla(l$theta1,l$theta2,l$theta0, l$sigma, l$n)
#' @export
#' 
nabla.to.standard <- function(r, nabla, theta1, theta2, n) {
  check.input.parameter(as.list(match.call())[-1])
  theta0 <- (1-r)*theta1 + r*theta2
  sigma <- (theta2-theta1) / (nabla*sqrt(2/n))
  list(theta1 = theta1, 
       theta2 = theta2, 
       theta0 = theta0, 
       sigma = sigma,
       n = n)
}


#' standard to TOST sensitivity parameterization
#' 
#' The standard parameterization is transformed to the TOST sensitivity parameterization.
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#' @param theta0 mean-difference
#' @param sigma standard deviation of the measurements
#' @param n sample size in each group
#' 
#' @return A list with the transformed parameters in TOST sensitivity parameterization, 
#' i.e. \eqn{(r, \nabla, \theta_1, \theta_2, n)}{(r, nabla, theta1, theta2, n)}
#' @examples
#' l <- nabla.to.standard(r = 0.75, nabla = 10, theta1 = 0, theta2 = 4, n = 20)
#' standard.to.nabla(l$theta1,l$theta2,l$theta0, l$sigma, l$n)
#' @export
#' 
standard.to.nabla <- function(theta1, theta2, theta0, sigma, n) {
  check.input.parameter(as.list(match.call())[-1])
  nabla <- (theta2-theta1) / (sigma*sqrt(2/n))
  r <- (theta0-theta1) / (theta2-theta1)
  list(r = r, 
       nabla = nabla, 
       theta1 = theta1, 
       theta2 = theta2,
       n = n)
}