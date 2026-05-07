#' Gibbs Sampler for Microbial Association with Left-censored
#' Missing Metabolic Outputs (MALMO)
#'
#' Implements a Gibbs sampler for Bayesian reduced rank regression
#' for microbiome-metabolite data integration with missing value
#' imputation. The model assumes Y = X * B * A^T + E where B and A
#' are low-rank factor matrices.
#'
#' @importFrom stats rnorm rbeta rbinom pnorm runif
#' @importFrom truncnorm rtruncnorm
#' @importFrom MCMCpack rinvgamma
#'
#' @param Y Numeric matrix of dimension n x q. Outcome matrix of
#'   metabolite measurements with missing values set to 0.
#' @param X Numeric matrix of dimension n x p. Standardized predictor
#'   matrix of microbiome features.
#' @param nburn Integer. Number of burn-in iterations. Default is 1000.
#' @param niter Integer. Number of posterior sampling iterations after
#'   burn-in. Default is 2000.
#' @param display Logical. If \code{TRUE}, prints progress every 1000
#'   iterations. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{B_result}{List of posterior samples of matrix B (p x q).}
#'   \item{A_result}{List of posterior samples of matrix A (q x q).}
#'   \item{C_result}{List of posterior samples of coefficient matrix
#'     C = B * A^T (p x q).}
#'   \item{Sigma_result}{List of posterior samples of diagonal
#'     covariance matrix Sigma (q x q).}
#'   \item{beta0_result}{List of posterior samples of intercept
#'     vector beta0 (length q).}
#'   \item{lambda2_result}{List of posterior samples of local
#'     shrinkage parameters (p x q).}
#'   \item{tau2_result}{List of posterior samples of global
#'     shrinkage parameters (length q).}
#'   \item{psi_result}{List of posterior samples of hyperparameters
#'     for tau2 (length q).}
#'   \item{nu_result}{List of posterior samples of hyperparameters
#'     for lambda2 (p x q).}
#'   \item{phi_result}{List of posterior probability matrices of
#'     each missing value being below the LOD (n x q).}
#'   \item{Z_result}{List of missing value classification matrices
#'     at each iteration: 0 = MAR, 1 = MNAR, 2 = observed (n x q).}
#'   \item{Y_mvi_result}{List of imputed outcome matrices at each
#'     iteration (n x q).}
#'   \item{theta_result}{Numeric vector of posterior samples of the
#'     MAR probability parameter theta.}
#' }
#'
#' @seealso \code{\link{malmo_summarize}}, \code{\link{gen_missing_value_mv}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## Generate simulation data
#' set.seed(123)
#' mydata <- gen_missing_value_mv(
#'   n = 100, p = 200, q = 3, snr = 1,
#'   proportion = 0.20, type = "mnar0.33"
#' )
#'
#' ## Prepare inputs
#' X <- mydata$X
#' Y <- mydata$Y_miss
#'
#' ## Run MALMO Gibbs sampler
#' result <- malmo_gibbs(
#'   Y      = Y,
#'   X      = X,
#'   nburn  = 1000,
#'   niter  = 2000,
#'   display = FALSE
#' )
#'
#' ## Summarize results
#' C_summary <- malmo_summarize(result$C_result, nburn = 1000, niter = 2000)
#' C_est     <- C_summary$mean
#'
#' ## Compare with truth
#' cat("True non-zero features:", mydata$p_true, "\n")
#' cat("Estimated non-zero rows:", which(rowSums(abs(C_est)) > 0), "\n")
#' }
malmo_gibbs <- function(Y, X, nburn = 1000, niter = 2000,
                        display = FALSE) {

  ## Dimensions ----
  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  ## Storage ----
  nu_result     <- vector("list", nburn + niter)
  psi_result    <- vector("list", nburn + niter)
  lambda2_result <- vector("list", nburn + niter)
  tau2_result   <- vector("list", nburn + niter)
  B_result      <- vector("list", nburn + niter)
  A_result      <- vector("list", nburn + niter)
  C_result      <- vector("list", nburn + niter)
  Sigma_result  <- vector("list", nburn + niter)
  beta0_result  <- vector("list", nburn + niter)
  phi_result    <- vector("list", nburn + niter)
  Z_result      <- vector("list", nburn + niter)
  Y_mvi_result  <- vector("list", nburn + niter)
  theta_result  <- numeric(nburn + niter + 1)

  ## Initialize parameters ----

  # nu: IG(1/2, 1)
  nu <- matrix(rinvgamma(p * q, shape = 0.5, scale = 1),
               nrow = p, ncol = q, byrow = FALSE)

  # lambda2: IG(1/2, 1/nu)
  lambda2 <- matrix(rinvgamma(p * q, shape = 0.5, scale = c(1/nu)),
                    nrow = p, ncol = q, byrow = FALSE)

  # psi: IG(1/2, 1)
  psi <- rinvgamma(q, shape = 0.5, scale = 1)

  # tau2: IG(1/2, 1/psi)
  tau2 <- rinvgamma(q, shape = 0.5, scale = 1/psi)

  # Sigma: diagonal matrix
  Sigma <- diag(runif(q, min = 0.5, max = 1), nrow = q, ncol = q)

  # A: q x q matrix
  A <- matrix(rnorm(q * q, mean = 0, sd = 1),
              nrow = q, ncol = q, byrow = FALSE)

  # beta0: intercept vector
  beta0 <- rnorm(q, mean = 0, sd = sqrt(diag(Sigma)))

  ## Missingness setup ----
  R   <- matrix(0, nrow = n, ncol = q)
  R[Y != 0] <- 1
  LOD <- min(Y[Y != 0])

  theta          <- runif(1, 0, 1)
  theta_result[1] <- theta

  ## Main Gibbs loop ----
  for (i in 1:(nburn + niter)) {

    ## Sample B ----
    Y_tilde       <- Y - matrix(1, nrow = n, ncol = 1) %*% t(beta0)
    y             <- as.vector(t(Y_tilde))
    Sigma_tilde   <- rep(diag(Sigma), times = n)
    Sigma_tilde_inv <- 1 / Sigma_tilde
    sqrt_Sigma_inv  <- sqrt(Sigma_tilde_inv)
    y_tilde         <- sqrt_Sigma_inv * y
    Lambda          <- as.vector(t(sweep(lambda2, 2, tau2, `*`)))

    # Sample u and delta
    u     <- rnorm(p * q, mean = 0, sd = sqrt(Lambda))
    delta <- rnorm(n * q)
    U     <- matrix(u, nrow = q, ncol = p, byrow = FALSE)
    v     <- sqrt_Sigma_inv * as.vector(A %*% U %*% t(X)) + delta

    # Solve for w
    M0           <- fast_kron_diag_kronT(X, A, Lambda)
    M1           <- sweep(M0, 1, sqrt_Sigma_inv, FUN = "*")
    middle_mat   <- sweep(M1, 2, sqrt_Sigma_inv, FUN = "*") + diag(n * q)
    middle_chol  <- chol(middle_mat)
    w            <- backsolve(middle_chol,
                              forwardsolve(t(middle_chol), (y_tilde - v)))

    W    <- matrix(sqrt_Sigma_inv * w, nrow = q, ncol = n, byrow = FALSE)
    beta <- u + Lambda * as.vector(t(A) %*% W %*% X)
    B    <- matrix(beta, nrow = p, ncol = q, byrow = TRUE)
    B_result[[i]] <- B

    ## Sample A ----
    XtX      <- crossprod(X)
    BtXtXB   <- crossprod(B, XtX %*% B)
    Sigma_inv <- diag(1/diag(Sigma))
    omega_a   <- kron_function(BtXtXB, Sigma_inv) + diag(q^2)
    omega_chol <- t(chol(omega_a))
    v_a  <- forwardsolve(omega_chol,
                         as.vector(Sigma_inv %*% t(Y_tilde) %*% (X %*% B)))
    m_a  <- backsolve(t(omega_chol), v_a)
    z_a  <- rnorm(q^2)
    w_a  <- backsolve(t(omega_chol), z_a)
    a    <- m_a + w_a
    A    <- matrix(a, nrow = q, ncol = q, byrow = FALSE)
    A_result[[i]] <- A

    ## Save C = B * A^T ----
    C             <- tcrossprod(B, A)
    C_result[[i]] <- C

    ## Sample shrinkage parameters ----
    # lambda2
    scale_lambda2 <- (B^2) / (2 * matrix(tau2, nrow = p, ncol = q,
                                         byrow = TRUE)) + 1/nu
    lambda2 <- matrix(rinvgamma(p * q, shape = 1,
                                scale = c(scale_lambda2)),
                      nrow = p, ncol = q)
    lambda2_result[[i]] <- lambda2

    # nu
    scale_nu <- 1 + 1/lambda2
    nu <- matrix(rinvgamma(p * q, shape = 1, scale = c(scale_nu)),
                 nrow = p, ncol = q)
    nu_result[[i]] <- nu

    # tau2
    scale_tau2 <- 0.5 * colSums(B^2 / lambda2) + 1/psi
    tau2 <- rinvgamma(q, shape = (p + 1)/2, scale = scale_tau2)
    tau2_result[[i]] <- tau2

    # psi
    scale_psi <- 1 + 1/tau2
    psi <- rinvgamma(q, shape = 1, scale = scale_psi)
    psi_result[[i]] <- psi

    ## Sample beta0 ----
    Y_star     <- Y - X %*% C
    mean_beta0 <- colSums(Y_star) / (n + 1)
    cov_beta0  <- diag(Sigma) / (n + 1)
    beta0      <- rnorm(q, mean = mean_beta0, sd = sqrt(cov_beta0))
    beta0_result[[i]] <- beta0

    ## Sample Sigma ----
    Residuals      <- Y - matrix(1, nrow = n, ncol = 1) %*% t(beta0) -
      X %*% C
    scale_sigma2   <- 0.5 * colSums(Residuals^2)
    sigma2_values  <- rinvgamma(q, shape = n/2, scale = scale_sigma2)
    Sigma          <- diag(sigma2_values, nrow = q, ncol = q)
    Sigma_result[[i]] <- Sigma

    ## Missing value imputation ----
    mu_matrix  <- matrix(1, nrow = n, ncol = 1) %*% t(beta0) + X %*% C
    sigma_sds  <- matrix(sqrt(diag(Sigma)), nrow = n, ncol = q,
                         byrow = TRUE)
    is_miss    <- (R == 0)

    # phi: P(Y < LOD | observed)
    phi <- matrix(NA_real_, n, q)
    phi[is_miss] <- pnorm(LOD,
                          mean       = mu_matrix[is_miss],
                          sd         = sigma_sds[is_miss],
                          lower.tail = TRUE)
    phi_result[[i]] <- phi

    # phi_value: P(Z = 1 | missing)
    phi_value <- matrix(NA_real_, n, q)
    phi_value[is_miss] <- 1 / (1 + theta *
                                 ((1 - phi[is_miss]) / phi[is_miss]))

    # Z: 2 = observed, 1 = MNAR, 0 = MAR
    Z <- matrix(2L, n, q)
    Z[is_miss] <- rbinom(sum(is_miss), 1, phi_value[is_miss])
    Z_result[[i]] <- Z

    # Impute Y_mvi
    Y_mvi <- matrix(NA_real_, n, q)

    mask_low <- (Z == 1)
    if (any(mask_low)) {
      Y_mvi[mask_low] <- rtruncnorm(
        sum(mask_low), a = -Inf, b = LOD,
        mean = mu_matrix[mask_low],
        sd   = sigma_sds[mask_low]
      )
    }

    mask_high <- (Z == 0)
    if (any(mask_high)) {
      Y_mvi[mask_high] <- rtruncnorm(
        sum(mask_high), a = LOD, b = Inf,
        mean = mu_matrix[mask_high],
        sd   = sigma_sds[mask_high]
      )
    }

    Y_mvi[Z == 2]     <- Y[Z == 2]
    Y_mvi_result[[i]] <- Y_mvi

    # Update theta
    theta            <- rbeta(1,
                              shape1 = sum(Z == 0) + 1,
                              shape2 = sum(Z == 2) + 1)
    theta_result[i]  <- theta

    # Update Y
    Y <- Y_mvi

    if (display && i %% 1000 == 0) {
      cat("Iteration:", i, "\n")
    }
  }

  if (display) cat("MALMO Gibbs Sampling complete.\n")

  ## Return ----
  return(list(
    B_result      = B_result,
    A_result      = A_result,
    C_result      = C_result,
    Sigma_result  = Sigma_result,
    beta0_result  = beta0_result,
    lambda2_result = lambda2_result,
    tau2_result   = tau2_result,
    psi_result    = psi_result,
    nu_result     = nu_result,
    phi_result    = phi_result,
    Z_result      = Z_result,
    Y_mvi_result  = Y_mvi_result,
    theta_result  = theta_result
  ))
}
