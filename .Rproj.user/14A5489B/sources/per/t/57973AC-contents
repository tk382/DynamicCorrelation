#' Computes the small sample correction
#'
#' @param y main vector to test against
#' @param Y matrix of target vectors
#' @param d degree statistic
#' @param numsim number of simulations for null distribution
#' @return p-value
#' @export
get_p_from_degree = function(y, Y, d, numsim = 5000) {
  bigy = cbind(y, Y)
  K = ncol(bigy)
  est_Sigma = stats::cor(bigy)
  est_H = get_H(est_Sigma)
  lambda = eigen(est_H)$values
  U      = eigen(est_H)$vectors
  null_d = matrix(NA, numsim, K-1)
  for (k in 1:(K-1)){
    null_d[,k] = rgamma(numsim, 1/2, 1/(2*lambda[k]))
  }
  p = sum(rowSums(null_d) > d)/numsim
  return(p)
}

#' Computes the small sample correction
#'
#' @param score score statistic computed from get_score
#' @param coef coefficients for the cubic function
#' @return small-sample adjusted score statistic
#' @export
post_score = function(score, coef) {
  roots = polyroot(c(-score, coef))
  return(Re(roots)[abs(Im(roots)) < 1e-06][1])
}


#' Helper function for get_est_H
#' Returns each element of matrix H
#'
#' @param rho12 correlation between variable 1 and variable 2
#' @param rho23 correlation between variable 2 and variable 3
#' @param rho13 correlation between variable 1 and variable 3
#' @return element (j,k) of H from elements (i,j), (j,k), (i,k) of Sigma
#' @export
get_cor_r = function(rho12, rho23, rho13) {
  num = (rho23 + 2 * rho12 * rho23) *
    (rho12^2 + 1) * (rho13^2 + 1)
  num = num + rho12 * rho13 * (6 + 2 *rho12 + 2 * rho13 + 2 * rho23)
  num = num - rho12 * (rho13^2 + 1) * (3 * rho13 + rho13 + 2 * rho12 *rho23)
  num = num - rho13 * (rho12^2 + 1) * (3 * rho12 + rho12 + 2 * rho13 * rho23)
  denom = (1 - rho12^2) * (1 - rho13^2) * sqrt(1 + rho12^2) * sqrt(1 + rho13^2)
  return(num/denom)
}


#' Computes the correlation matrix of the combined test statistics r
#' Takes input of the estimated or true covariance matrix
#'
#' @param Sigma estimate or true covariance matrix of $K$ variables
#' @return estimated H from the covariance matrix Sigma
#' @export
get_est_H = function(Sigma) {
  K = nrow(Sigma)
  est_H = matrix(0, K - 1, K - 1)
  diag(est_H) = 1
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      eta = get_cor_r(Sigma[1, i],
                      Sigma[i, j], Sigma[j, 1])
      est_H[i - 1, j - 1] = est_H[j -
                                    1, i - 1] = eta
    }
  }
  return(est_H)
}




#' #' Returns rho
#' #'
#' #' @param X scaled covariate
#' #' @param alpha length 2 vector for intercept + slope for X
#' #' @return simulated rho
#' #' @export
#' h1_fisher = function(X, alpha) {
#'   tmp = X %*% alpha
#'   return((exp(tmp) - 1)/(exp(tmp) + 1))
#' }
#'
#' #' Returns rho
#' #'
#' #' @param X scaled covariate
#' #' @param alpha length 2 vector for intercept + slope for X
#' #' @return simulated rho
#' #' @export
#' h2_sqrt = function(X, alpha) {
#'   tmp = (X %*% alpha)/sqrt(1 + (X %*%
#'                                   alpha)^2)
#'   return(tmp * 2 - 1)
#' }
#'
#'
#' #' Returns rho
#' #'
#' #' @param X scaled covariate
#' #' @param alpha length 2 vector for intercept + slope for X
#' #'
#' #' @export
#' h3_cdf = function(X, alpha) {
#'   tmp = pnorm(X %*% alpha, 0, 10) * 2 -
#'     1
#'   return(tmp)
#' }
#'
#' #' Returns rho
#' #'
#' #' @param X scaled covariate
#' #' @param alpha length 2 vector for intercept + slope for X
#' #' @return simulated rho
#' #' @export
#' h4_sin = function(X, alpha) {
#'   tmp = sin(X %*% alpha * 2)
#'   return(tmp)
#' }
#'
#' #' #' Returns rho
#' #' #'
#' #' #' @param X scaled covariate
#' #' #' @param alpha length 2 vector for intercept + slope for X
#' #' #' @return simulated rho
#' #' #' @export
#' #' h5_gumbel = function(X, alpha) {
#' #'   tmp = pgumbel(X %*% alpha, loc = 1,
#' #'                 scale = 2)
#' #'   return(tmp)
#' #' }
#'
#'
#' #' Returns rho
#' #'
#' #' @param X scaled covariate
#' #' @param alpha length 2 vector for intercept + slope for X
#' #' @return simulated rho
#' #' @export
#' h6_quadratic = function(X, alpha) {
#'   sigma1 = (X %*% alpha - 0.1)^2 - 0.99
#'   # sigma1 = pmin(sigma1, 0.9) sigma1 =
#'   # pmax(sigma1, -0.9)
#'   return(sigma1)
#' }
