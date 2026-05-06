#' Generate Simulation Predictors
#'
#' This function generates simulated microbiome data with correlated covariates
#' for Bayesian sparse regression models, using a multivariate normal
#' distribution to introduce correlation structure among features.
#'
#' @importFrom stats rnorm
#' @importFrom MASS mvrnorm
#'
#' @param n Integer. Number of observations.
#' @param gamma Numeric vector of length p. Binary inclusion indicators
#'   for each covariate (1 = included, 0 = excluded).
#' @param b Numeric vector of length p. Coefficients for included covariates.
#' @param theta Numeric vector of length p. Column-wise mean for each
#'   microbiome feature. Used to mimic real microbiome data where different
#'   features have different abundance magnitudes.
#' @param sigmaX Numeric vector of length p. Standard deviation for each
#'   covariate, used to scale the covariance matrix.
#' @param Xcor Matrix of dimensions p x p. Correlation matrix for the
#'   covariates, specifying the correlation structure among microbiome features.
#' @param sigma Numeric. Standard deviation for the noise term epsilon.
#'
#' @return A list containing:
#' \describe{
#'   \item{X}{An n x p matrix of generated correlated covariates.}
#'   \item{beta}{A numeric vector of length p, the element-wise product
#'     of gamma and b.}
#'   \item{epsilon}{A numeric vector of length n, the generated noise.}
#' }
#'
#' @export
#'
#' @examples
#' ## Generate simulation data with correlated covariates
#' ## n = 100 observations, p = 10 microbiome features
#'
#' n <- 50
#' p <- 100
#' gamma <- rbinom(p, 1, 0.3)
#' b <- rnorm(p, mean = 0, sd = 1)
#' theta <- rep(0, p)
#' sigmaX <- rep(1, p)
#' rho <- 0.5
#' Xcor <- rho^abs(outer(1:p, 1:p, "-"))
#' sigma <- 0.5
#' result <- gen_simdata(n, gamma, b, theta, sigmaX, Xcor, sigma)
#' str(result)
gen_simdata <- function(n, gamma, b, theta, sigmaX, Xcor, sigma) {
  beta <- gamma * b
  X <- MASS::mvrnorm(n = n, mu = theta,
                     Sigma = diag(sigmaX) %*% Xcor %*% diag(sigmaX))
  epsilon <- sigma * rnorm(n)
  return(list(X = X, beta = beta, epsilon = epsilon))
}
