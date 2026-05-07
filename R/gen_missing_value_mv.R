#' Generate Multivariate Simulation Data with Missing Values
#'
#' Generates simulated microbiome data with multivariate metabolite
#' outcomes and missing values, using a block correlation structure
#' to mimic real microbiome data. This function is designed for the
#' Microbial Association with Left-censored Missing Metabolic Outputs
#' (MALMO) model.
#'
#' @importFrom stats rnorm quantile runif
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#'
#' @param n Integer. Number of observations. Must be at least 10.
#' @param p Integer. Number of microbiome features. Must be at least 10.
#' @param q Integer. Number of metabolite outcomes. Must be at least 2.
#' @param snr Numeric. Signal-to-noise ratio controlling the noise level.
#'   Must be positive.
#' @param proportion Numeric. Proportion of missing values in the
#'   outcome matrix (between 0 and 1, exclusive).
#' @param type Character. Missing data mechanism, one of:
#'   \itemize{
#'     \item \code{"mnar"} - Missing Not At Random (missing due to
#'       concentration values lower than the limit of detection.)
#'     \item \code{"mar"} - Missing At Random (value is greater than
#'       LOD but missing due to technical issues.)
#'     \item \code{"mnar0.33"} - Mixed: 1/3 MNAR + 2/3 MAR
#'     \item \code{"mnar0.67"} - Mixed: 2/3 MNAR + 1/3 MAR
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{Y_miss}{An n x q matrix of outcomes with missing values
#'     set to 0.}
#'   \item{Y_true}{An n x q matrix of true outcomes without missing
#'     values.}
#'   \item{C}{A p x q coefficient matrix. Non-zero entries correspond
#'     to true associations between microbiome features and outcomes.}
#'   \item{C_idx}{A p x q binary matrix. 1 indicates a true non-zero
#'     coefficient, 0 otherwise.}
#'   \item{X}{An n x p standardized microbiome feature matrix.}
#'   \item{p_true}{Integer vector of true signal feature indices.}
#'   \item{data_type}{Character. The missing data mechanism used.}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' mydata <- gen_missing_value_mv(
#'   n          = 50,
#'   p          = 100,
#'   q          = 3,
#'   snr        = 1,
#'   proportion = 0.20,
#'   type       = "mnar0.33"
#' )
#' str(mydata)
gen_missing_value_mv <- function(
    n,
    p,
    q,
    snr,
    proportion,
    type = c("mnar", "mar", "mnar0.33", "mnar0.67")
) {

  ## Input validation ----
  if (n < 10) {
    stop("n must be at least 10 observations.")
  }
  if (p < 10) {
    stop("p must be at least 10 features.")
  }
  if (q < 2) {
    stop("q must be at least 2 outcomes.")
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
  p_true     <- c(block1_idx, block2_idx)

  ## Noise block indices ----
  n_noise      <- max(2, round(0.016 * p))
  noise1_start <- round(0.444 * p)
  noise1_end   <- min(noise1_start + n_noise - 1, p)
  noise2_start <- round(0.944 * p)
  noise2_end   <- min(noise2_start + n_noise - 1, p)

  noise1_idx <- setdiff(noise1_start:noise1_end, p_true)
  noise2_idx <- setdiff(noise2_start:noise2_end, p_true)
  noise_idx  <- c(noise1_idx, noise2_idx)

  ## Coefficient matrix C ----
  C <- matrix(0, nrow = p, ncol = q)
  C[p_true, ] <- matrix(
    round(
      sample(c(1, -1), length(p_true) * q, replace = TRUE) *
        runif(length(p_true) * q, 1, 2),
      2
    ),
    nrow = length(p_true),
    ncol = q
  )


  # Binary indicator matrix
  C_idx           <- matrix(0, nrow = p, ncol = q)
  C_idx[p_true, ] <- 1

  ## Covariance matrix for X ----
  Xcor <- diag(0, p)

  # Signal block 1
  for (i in block1_idx)
    for (j in block1_idx[block1_idx > i])
      Xcor[i, j] <- max(0, 0.75 - 0.0015 * abs(i - j))

  # Signal block 2
  for (i in block2_idx)
    for (j in block2_idx[block2_idx > i])
      Xcor[i, j] <- max(0, 0.75 - 0.0015 * abs(i - j))

  # Noise block 1
  for (i in noise1_idx)
    for (j in noise1_idx[noise1_idx > i])
      Xcor[i, j] <- max(0, 0.4 - 0.02 * abs(i - j))

  # Noise block 2
  for (i in noise2_idx)
    for (j in noise2_idx[noise2_idx > i])
      Xcor[i, j] <- max(0, 0.4 - 0.02 * abs(i - j))

  Xcor <- Xcor + t(Xcor) + diag(1, p)

  # Fix positive definiteness
  Xcor <- as.matrix(Matrix::nearPD(Xcor, corr = TRUE)$mat)

  ## Mean vector theta ----
  theta            <- rep(0, p)
  theta[p_true]    <- log(0.5 * p)
  theta[noise_idx] <- log(0.25 * p)

  ## Generate X ----
  W <- MASS::mvrnorm(n, mu = theta, Sigma = Xcor)
  X <- scale(W, center = TRUE, scale = TRUE)

  ## Error matrix ----
  sigma <- mean(abs(C[p_true, ])) / snr
  E     <- matrix(rnorm(n * q, mean = 0, sd = sigma),
                  nrow = n, ncol = q)

  ## True outcome matrix ----
  Y_true <- X %*% C + E

  ## Missing value generation ----
  Y_miss <- Y_true

  if (type == "mnar") {
    lod         <- quantile(Y_true, probs = proportion)
    missing_idx <- which(Y_miss < lod)
    Y_miss[missing_idx] <- 0

  } else if (type == "mar") {
    total_obs   <- n * q
    missing_idx <- sample(total_obs,
                          size = floor(total_obs * proportion))
    Y_miss[missing_idx] <- 0

  } else if (type == "mnar0.33") {
    # 1/3 MNAR
    lod              <- quantile(Y_true, probs = proportion * (1/3))
    Y_miss[Y_miss < lod] <- 0

    # 2/3 MAR from remaining
    remain_idx <- which(Y_miss != 0)
    mar_size   <- ceiling(n * q * proportion * (2/3))
    mar_size   <- min(mar_size, length(remain_idx))
    miss_idx   <- sample(remain_idx, size = mar_size)
    Y_miss[miss_idx] <- 0

  } else if (type == "mnar0.67") {
    # 2/3 MNAR
    lod              <- quantile(Y_true, probs = proportion * (2/3))
    Y_miss[Y_miss < lod] <- 0

    # 1/3 MAR from remaining
    remain_idx <- which(Y_miss != 0)
    mar_size   <- ceiling(n * q * proportion * (1/3))
    mar_size   <- min(mar_size, length(remain_idx))
    miss_idx   <- sample(remain_idx, size = mar_size)
    Y_miss[miss_idx] <- 0
  }

  cat(sprintf("Missing proportion: %.3f (target: %.3f)\n",
              mean(Y_miss == 0), proportion))

  ## Return ----
  return(list(
    Y_miss    = Y_miss,
    Y_true    = Y_true,
    C         = C,
    C_idx     = C_idx,
    X         = X,
    p_true    = p_true,
    data_type = type
  ))
}
