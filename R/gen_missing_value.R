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
#'   Must be positive.
#' @param proportion Numeric. Proportion of missing values (between 0 and 1,
#'   exclusive).
#' @param type Character. Missing data mechanism, one of:
#'   \itemize{
#'     \item \code{"mnar"} - Missing Not At Random (missing due to
#'       concentration values lower than the limit of detection.)
#'     \item \code{"mar"} - Missing At Random (value is greater than LOD
#'       but missing due to technical issues.)
#'     \item \code{"mnar0.33"} - Mixed: 1/3 MNAR + 2/3 MAR
#'     \item \code{"mnar0.67"} - Mixed: 2/3 MNAR + 1/3 MAR
#'   }
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
#'   \item{data_type}{Character. The missing data mechanism used.}
#'   \item{Q}{A p x p matrix encoding prior correlation information for
#'     the Bayesian sparse regression model.}
#' }
#'
#' @export
#'
#' @examples
#' ## Generate correlated microbiome data with MNAR0.33 missing mechanism
#' set.seed(123)
#' mydata <- gen_missing_value(
#'   n = 100, p = 200, snr = 1,
#'   proportion = 0.30,
#'   type = "mnar0.33"
#' )
#' str(mydata)
#' which(mydata$coeffs[, 1] == 1)
#'
#' ## MAR missing mechanism
#' set.seed(123)
#' mydata_mar <- gen_missing_value(
#'   n = 100, p = 200, snr = 1,
#'   proportion = 0.30,
#'   type = "mar"
#' )
#' str(mydata_mar)
gen_missing_value <- function(
    n,
    p,
    snr,
    proportion,
    type = c("mnar", "mar", "mnar0.33", "mnar0.67")
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
  if (!type %in% c("mnar", "mar", "mnar0.33", "mnar0.67")) {
    stop('type must be one of "mnar", "mar", "mnar0.33", "mnar0.67".')
  }

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

  if (type == "mnar") {
    lod         <- quantile(Y_true, probs = proportion)
    missing_idx <- which(Y_miss < lod)

  } else if (type == "mar") {
    lod         <- min(Y_true)
    missing_idx <- sample(n, floor(n * proportion))

  } else if (type == "mnar0.33") {
    cutoff      <- proportion / 3
    lod         <- quantile(Y_true, probs = cutoff)
    mnar_idx    <- which(Y_true < lod)
    remaining   <- setdiff(1:n, mnar_idx)
    mar_idx     <- sample(remaining, floor((proportion - cutoff) * n))
    missing_idx <- c(mnar_idx, mar_idx)

  } else if (type == "mnar0.67") {
    cutoff      <- (2 * proportion) / 3
    lod         <- quantile(Y_true, probs = cutoff)
    mnar_idx    <- which(Y_true < lod)
    remaining   <- setdiff(1:n, mnar_idx)
    mar_idx     <- sample(remaining, floor((proportion - cutoff) * n))
    missing_idx <- c(mnar_idx, mar_idx)
  }

  Y_miss[missing_idx] <- 0

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
    data_type = type,
    Q         = Q
  ))
}
