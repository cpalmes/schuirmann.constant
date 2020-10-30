# Description of intern package function uniroot2:
# Calculates the root of an univariate function f using the R build-in stats::uniroot
# function. The uniroot function needs a lower and upper interval limit where
# f has different signs, i.e. f(lower) * f(upper) < 0. If this assumption is not
# valid, uniroot cannot calculate the root and breaks with an error message. The
# function uniroot2 generalizes this restriction by searching for a change of sign
# in the following way:
# (i) case (to.one = TRUE):
# The change of sign is searched for the intervals [lower, 1 - (1-upper)/10^k]
# for k = 0, .., k.max
# (ii) case (to.one = FALSE):
# The change of sign is searched for the intervals [lower, upper * 10^k]
# for k = 0, .., k.max
# If such a change of sign cannot be found, then the function ends with an 
# informative error message.
#
uniroot2 <- function(fct, lower, upper, to.one, k.max = 15) {
  upper.modified <- upper
  k <- 1
  while(k <= k.max && fct(lower) * fct(upper.modified) > 0) {
    upper.modified <- ifelse(to.one, 1 - (1-upper)/(10^k), upper*(10^k))
    k <- k + 1
  }
  if (k > k.max) {
    stop(paste0('Calculation not possible due to ill-conditioned numerics! \n Please choose meaningful',
         ' parameters such as 0.001 <= alpha <= 0.200 and 0.800 <= pwr <= 0.999!'))
  } else {
    uniroot(f = fct, lower = lower, upper = upper.modified)
  }
}
