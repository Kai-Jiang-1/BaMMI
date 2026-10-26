#' Generate Simulation Data with Missing Values
#'
#' Generates simulated microbiome compositional data with correlated
#' covariate structure and missing values in the outcome, mimicking
#' real-world microbiome data with block correlation structure and
#' multiple missing data mechanisms.
#'
#' @importFrom stats rnorm quantile
#' @importFrom Matrix nearPD
#'
#' @param n Integer. Number of observations. Must be at least 10.
#' @param p Integer. Number of microbiome features. Must be at least 10.
#' @param snr Numeric. Signal-to-noise ratio controlling the noise level.
#' @param proportion Numeric. Proportion of missing values (between 0 and 1,
#'   exclusive).
#' @param mnar_ratio Numeric between 0 and 1. Fraction of the missing
#'   values that are Missing Not At Random (MNAR, i.e. below the limit of
#'   detection); the remaining \code{1 - mnar_ratio} are Missing At Random
#'   (MAR, above the LOD but missing due to technical issues). For example,
#'   \code{mnar_ratio = 1} gives pure MNAR, \code{0} gives pure MAR and
#'   \code{0.25} gives 25\% MNAR + 75\% MAR. If supplied, it overrides
#'   \code{type}.
#' @param type Character. Legacy shortcut for \code{mnar_ratio}, one of
#'   \code{"mnar"} (1), \code{"mar"} (0), \code{"mnar0.33"} (1/3) or
#'   \code{"mnar0.67"} (2/3). Ignored when \code{mnar_ratio} is given.
#'
#' @return A list containing:
#' \describe{
#'   \item{coeffs}{A matrix with two columns: \code{gammatrue} (binary
#'     inclusion indicators) and \code{beta} (true coefficients).}
#'   \item{Y_true}{A numeric vector of length n. True outcome without
#'     missing values.}
#'   \item{X}{An n x p matrix of log-ratio transformed compositional
#'     microbiome features.}
#'   \item{Y_miss}{A numeric vector of length n. Outcome with missing
#'     values set to 0.}
#'   \item{data_type}{Character. Label of the missing data mechanism
#'     used, e.g. \code{"mnar0.25"}.}
#'   \item{mnar_ratio}{Numeric. Fraction of missing values that are MNAR.}
#'   \item{mnar_idx}{Integer vector. Indices of the MNAR missing values.}
#'   \item{mar_idx}{Integer vector. Indices of the MAR missing values.}
#'   \item{Q}{A p x p matrix encoding prior correlation information for
#'     the Bayesian sparse regression model.}
#'   \item{sigma}{Numeric. Noise standard deviation used for the outcome.}
#' }
#'
#' @export
#'
#' @examples
#' ## Generate correlated microbiome data with 1/3 MNAR + 2/3 MAR missingness
#' set.seed(123)
#' mydata <- gen_missing_value(
#'   n = 200, p = 200, snr = 1,
#'   proportion = 0.30,
#'   mnar_ratio = 1/3
#' )
#' str(mydata)
#' which(mydata$coeffs[, 1] == 1)
#'
#' ## Pure MAR missingness
#' set.seed(123)
#' mydata_mar <- gen_missing_value(
#'   n = 200, p = 200, snr = 1,
#'   proportion = 0.30,
#'   mnar_ratio = 0
#' )
#' str(mydata_mar)
gen_missing_value <- function(
    n,
    p,
    snr,
    proportion,
    mnar_ratio = NULL,
    type       = NULL
) {

  ## Input validation ----
  if (p < 10) {
    stop("p must be at least 10 to generate meaningful simulation data.")
  }
  if (n < 10) {
    stop("n must be at least 10 observations.")
  }
  if (snr <= 0) {
    stop("snr must be positive.")
  }
  if (proportion <= 0 || proportion >= 1) {
    stop("proportion must be strictly between 0 and 1.")
  }
  mech <- .resolve_mnar_ratio(type, mnar_ratio)

  ## Block indices (proportional to p) ----
  step         <- max(1, round(20 * p / 1000))
  block1_start <- round(0.18 * p)
  block1_end   <- round(0.40 * p)
  block2_start <- round(0.58 * p)
  block2_end   <- round(0.80 * p)

  block1_idx <- seq(block1_start, block1_end, by = step)
  block2_idx <- seq(block2_start, block2_end, by = step)
  true_index <- c(block1_idx, block2_idx)

  ## Noise block indices ----
  n_noise      <- max(2, round(0.016 * p))
  noise1_start <- round(0.444 * p)
  noise1_end   <- min(noise1_start + n_noise - 1, p)
  noise2_start <- round(0.944 * p)
  noise2_end   <- min(noise2_start + n_noise - 1, p)

  noise1_idx <- setdiff(noise1_start:noise1_end, true_index)
  noise2_idx <- setdiff(noise2_start:noise2_end, true_index)
  noise_idx  <- c(noise1_idx, noise2_idx)

  ## gamma and b ----
  gammatrue             <- rep(0, p)
  gammatrue[true_index] <- 1

  b             <- rep(0, p)
  b[true_index] <- rnorm(length(true_index), mean = 0, sd = 1)
  b[true_index] <- ifelse(
    abs(b[true_index]) < 0.3,
    sign(b[true_index]) * 0.3,
    b[true_index]
  )
  b[true_index] <- b[true_index] - mean(b[true_index])

  sigma <- mean(abs(b[b != 0])) / snr

  b <- 10 * b

  ## theta (feature means) ----
  theta             <- rep(0, p)
  theta[true_index] <- log(0.5 * p)
  theta[noise_idx]  <- log(0.25 * p)

  ## Generate correlated X ----
  sigmaX <- rep(1, p)
  Xcor   <- diag(0, p)

  # Signal block 1 - actual index distance
  for (i in block1_idx)
    for (j in block1_idx[block1_idx > i])
      Xcor[i, j] <- max(0, 0.75 - 0.0015 * abs(i - j))

  # Signal block 2 - actual index distance
  for (i in block2_idx)
    for (j in block2_idx[block2_idx > i])
      Xcor[i, j] <- max(0, 0.75 - 0.0015 * abs(i - j))

  # Noise block 1 - actual index distance
  for (i in noise1_idx)
    for (j in noise1_idx[noise1_idx > i])
      Xcor[i, j] <- max(0, 0.4 - 0.02 * abs(i - j))

  # Noise block 2 - actual index distance
  for (i in noise2_idx)
    for (j in noise2_idx[noise2_idx > i])
      Xcor[i, j] <- max(0, 0.4 - 0.02 * abs(i - j))

  Xcor <- Xcor + t(Xcor) + diag(1, p)

  # Fix positive definiteness using nearPD
  Xcor <- as.matrix(Matrix::nearPD(Xcor, corr = TRUE)$mat)

  sim_data <- gen_simdata(n, gammatrue, b, theta, sigmaX, Xcor, sigma)

  ## Compositional transformation ----
  Xorg    <- sim_data$X
  temp    <- exp(2 * Xorg)
  Z       <- t(apply(temp, 1, function(x) x / sum(x)))
  X       <- log(Z)
  beta    <- sim_data$beta
  epsilon <- sim_data$epsilon
  Y_true  <- X %*% b + epsilon

  ## Missing value generation ----
  Y_miss <- Y_true
  miss   <- .gen_missing_idx(Y_true, proportion, mech$mnar_ratio)
  Y_miss[miss$missing_idx] <- 0

  ## Q matrix ----
  Q <- 0.002 * (matrix(1, nrow = p, ncol = p) - diag(1, p))

  # Signal block 1
  for (i in block1_idx)
    for (j in block1_idx[block1_idx > i]) {
      Q[i, j] <- 4
      Q[j, i] <- 4
    }

  # Signal block 2
  for (i in block2_idx)
    for (j in block2_idx[block2_idx > i]) {
      Q[i, j] <- 4
      Q[j, i] <- 4
    }

  # Noise block 1
  for (i in noise1_idx)
    for (j in noise1_idx[noise1_idx > i]) {
      Q[i, j] <- 4
      Q[j, i] <- 4
    }

  # Noise block 2
  for (i in noise2_idx)
    for (j in noise2_idx[noise2_idx > i]) {
      Q[i, j] <- 4
      Q[j, i] <- 4
    }

  ## Return ----
  return(list(
    coeffs    = cbind(gammatrue, beta),
    Y_true    = Y_true,
    X         = X,
    Y_miss    = Y_miss,
    data_type  = mech$label,
    mnar_ratio = mech$mnar_ratio,
    mnar_idx   = miss$mnar_idx,
    mar_idx    = miss$mar_idx,
    Q         = Q,
    sigma     = sigma
  ))
}
