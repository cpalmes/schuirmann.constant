#' Calculation of the schuirmann-constant
#' 
#' Calculation of the schuirmann-constant, cf. \insertCite{palmes2020}{schuirmann.constant} 
#'
#' @param alpha type I error
#' @param pwr asymptotic weighted power of the TOST
#' @param density a-priori density for the mean-difference
#'
#' @return The Schuirmann constant
#' @examples 
#' rho <- apriori.density('U')
#' theta1 <- 0
#' theta2 <- 3
#' s <- schuirmann.constant(alpha = 0.05, pwr = 0.8, density = rho)
#' theta0 <- (1-s)*theta1 + s*theta2 
#' n <- n.point(alpha = 0.05, theta1 = theta1, theta2 = theta2, theta0 = theta0, sigma = 1, pwr = 0.8)
#' ## Note that the weighted power is calculated with a sample size that was found for a
#' ## point-difference. Due to the newly introduced duality concept it is -nevertheless-
#' ## very close to the target power 0.8. 
#' power.weighted(alpha = 0.05, theta1 = theta1, theta2 = theta2, density = rho, sigma = 1, n = n)
#' @references
#' \insertAllCited{}
#' @export
#' 
schuirmann.constant <- function(alpha, pwr, density) {
  check.input.parameter(as.list(match.call())[-1])
  nabla <- nabla.weighted(alpha, pwr, density)
  fct <- function(r) {nabla.point(alpha, pwr,r) - nabla}
  S <- uniroot2(fct, lower = 0.5, upper = 1-10^(-4), to.one = T)$root
  return(S)
}


#' Calculation of the uniform schuirmann-constant
#' 
#' Calculation of the uniform schuirmann-constant, cf. \insertCite{palmes2020}{schuirmann.constant} 
#'
#' The uniform density is assumed as a-priori distribution. This density is
#' unique with the property that each possible mean-difference is equally weighted.
#' It represents the lack of any information about the true mean-difference. Due
#' to its importance, this special setting is included as a separate function for 
#' convenience. Technical, this function is merely a wrapper that calls the 
#' schuirmann.constant function with the uniform distribution as a-priori density.    
#' 
#' @param alpha type I error
#' @param pwr asymptotic weighted power of the TOST
#' @param theta1 lower limit of equivalence interval
#' @param theta2 upper limit of equivalence interval
#'
#' @return The Schuirmann constant of the uniform a-priori distribtion is returned. 
#' If \eqn{[\theta_1,  \theta_2]}{[theta1,  theta2]} differs from \eqn{[0,1]}, then 
#' the appropriately scaled true-mean difference
#' \deqn{\theta_0 = (S-1) \cdot \theta_1 + S \cdot \theta_2}{
#' theta0 = (S-1)*theta1 + S*theta2}
#' is returned.
#' @examples
#' rho <- apriori.density('U')
#' theta1 <- 0
#' theta2 <- 3
#' s <- schuirmann.constant(alpha = 0.05, pwr = 0.8, density = rho)
#' (1-s)*theta1 + s*theta2
#' schuirmann.constant.uniform(alpha = 0.05, pwr = 0.8, theta1 = theta1, theta2 = theta2) 
#' @references
#' \insertAllCited{}
#' @export
#' 
schuirmann.constant.uniform <- function(alpha, pwr, theta1 = 0, theta2 = 1) {
  check.input.parameter(as.list(match.call())[-1])
  s <- schuirmann.constant(alpha, pwr, apriori.density('U'))
  (s-1)*theta1 + s*theta2
}