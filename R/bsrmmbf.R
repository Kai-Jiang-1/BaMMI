#' Bayes Factor Computation for Bayesian Sparse Regression
#'
#' Computes the Bayes Factor (BF) for variable inclusion or exclusion
#' in the Bayesian sparse regression model. This function is called
#' internally by \code{\link{bsrmmgibbs}} at each iteration of the
#' Gibbs sampler to evaluate whether a proposed variable should be
#' included or excluded from the model.
#'
#' @importFrom stats pnorm
#'
#' @param Y Numeric matrix of dimension n x 1. Outcome vector.
#' @param Xri Numeric matrix of dimension n x q. Current predictor matrix
#'   excluding the proposed variable.
#' @param Xi Numeric vector of length n. The proposed predictor variable.
#' @param XIi Numeric matrix of dimension n x (q+1). Combined predictor
#'   matrix including the proposed variable.
#' @param Tri Numeric matrix. Submatrix of the constraint matrix N
#'   excluding the proposed variable.
#' @param Ti Numeric vector. Column of the constraint matrix N
#'   corresponding to the proposed variable.
#' @param TIi Numeric matrix. Combined constraint matrix including the
#'   proposed variable.
#' @param invAri Numeric matrix of dimension q x q. Inverse of the
#'   current precision matrix excluding the proposed variable.
#' @param n Integer. Number of observations.
#' @param tau Numeric. Prior scale parameter for the regression
#'   coefficients.
#' @param nu Numeric. Prior degrees of freedom for the inverse-gamma
#'   prior on sigma squared.
#' @param omega Numeric. Prior scale parameter for the inverse-gamma
#'   prior on sigma squared.
#' @param flag Logical. If \code{TRUE}, computes BF for variable
#'   exclusion (proposed variable is currently in the model). If
#'   \code{FALSE}, computes BF for variable inclusion (proposed
#'   variable is currently outside the model).
#' @param keep List. Cached quantities from the previous iteration
#'   including \code{resAi} (residual sum of squares) and
#'   \code{sqrtdetinvAi} (log determinant of inverse precision matrix).
#'
#' @return A list containing:
#' \describe{
#'   \item{BF or BF1}{Numeric. The computed Bayes Factor. \code{BF} is
#'     returned when \code{flag = FALSE} (inclusion), \code{BF1} when
#'     \code{flag = TRUE} (exclusion).}
#'   \item{keep}{Updated list of cached quantities for the next
#'     iteration.}
#'   \item{invAi}{Numeric matrix. Updated inverse precision matrix
#'     (only returned when \code{flag = FALSE}).}
#'   \item{lri1}{Upper triangular Cholesky factor of \code{invAri}.}
#'   \item{Lri}{Lower triangular Cholesky factor of \code{invAri}.}
#'   \item{Lambda}{Numeric matrix. Centering matrix of dimension n x n.}
#'   \item{resAri}{Numeric. Residual sum of squares for the current
#'     model excluding the proposed variable.}
#'   \item{logratiodetT}{Numeric. Log ratio of determinants of the
#'     constraint matrices.}
#' }
#'
#' @seealso \code{\link{bsrmmgibbs}}
#'
#' @export
#'
#' @examples
#' ## bsrmmbf is typically called internally by bsrmmgibbs
#' ## See bsrmmgibbs for a complete working example
bsrmmbf <- function(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri,
                    n, tau, nu, omega, flag, keep) {

  keep  <- keep
  lri1  <- chol(invAri)
  Lri   <- t(lri1)
  Lambda <- diag(1, nrow = n, ncol = n) -
    matrix(1, nrow = n, ncol = n) / (n + 1)

  sqrtdetinvAri <- sum(log(diag(Lri)))
  resAri <- t(Y) %*% Lambda %*% Y -
    t(Y) %*% Lambda %*% Xri %*% invAri %*%
    t(Xri) %*% Lambda %*% Y

  if (flag) {
    ## ── Exclusion: proposed variable is in the model ──────────────────────
    sqrtdetinvAi <- keep$sqrtdetinvAi
    resAi        <- keep$resAi

    keep$resAri        <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri

    logratiodetT <- sum(log(diag(chol(t(Tri) %*% Tri)))) -
      sum(log(diag(chol(t(TIi) %*% TIi))))

    BF1 <- -log(tau) +
      (sqrtdetinvAri - sqrtdetinvAi) +
      logratiodetT +
      (n + nu) / 2 * log((nu * omega + resAri) /
                           (nu * omega + resAi))
    BF1 <- exp(BF1)

    return(list(
      BF1           = BF1,
      keep          = keep,
      lri1          = lri1,
      Lri           = Lri,
      Lambda        = Lambda,
      sqrtdetinvAi  = sqrtdetinvAi,
      resAri        = resAri,
      resAi         = resAi,
      logratiodetT  = logratiodetT
    ))

  } else {
    ## ── Inclusion: proposed variable is outside the model ─────────────────
    Srii <- t(Xri) %*% Lambda %*% Xi +
      tau^(-2) * (t(Tri) %*% Ti)
    sii  <- t(Xi) %*% Lambda %*% Xi +
      tau^(-2) * (t(Ti) %*% Ti)

    # Rank-1 update of inverse precision matrix
    v1 <- sqrt(1 / (sii * (1 - t(Srii) %*% invAri %*% Srii / sii)))
    v  <- v1[1, 1] * invAri %*% Srii

    A11 <- invAri + v %*% t(v)
    A12 <- -A11 %*% (Srii / sii[1, 1])
    A21 <- -(1 / sii) %*% t(Srii) %*% A11
    A22 <- 1 / sii[1, 1] +
      (1 / sii[1, 1]) %*% t(Srii) %*% A11 %*% Srii / sii[1, 1]

    invAi <- rbind(cbind(A11, A12), cbind(A21, A22))
    resAi <- t(Y) %*% Lambda %*% Y -
      t(Y) %*% Lambda %*% XIi %*% invAi %*%
      t(XIi) %*% Lambda %*% Y

    # Rank-1 Cholesky update
    tilLri <- t(chol(Lri %*% t(Lri) + v %*% t(v)))
    Lrii   <- forwardsolve(tilLri, A12)
    lii    <- sqrt(A22 - t(Lrii) %*% Lrii)
    Li     <- rbind(cbind(tilLri, matrix(0, nrow(tilLri), 1)),
                    cbind(t(Lrii), lii))

    sqrtdetinvAi <- sum(log(diag(Li)))

    keep$resAri        <- resAri
    keep$sqrtdetinvAri <- sqrtdetinvAri
    keep$resAi         <- resAi
    keep$sqrtdetinvAi  <- sqrtdetinvAi

    logratiodetT <- -sum(log(diag(chol(t(Tri) %*% Tri)))) +
      sum(log(diag(chol(t(TIi) %*% TIi))))

    BF <- -log(tau) +
      (sqrtdetinvAri - sqrtdetinvAi) +
      logratiodetT +
      (n + nu) / 2 * log((nu * omega + resAri) /
                           (nu * omega + resAi))
    BF <- exp(BF)

    return(list(
      BF           = BF,
      keep         = keep,
      invAi        = invAi,
      lri1         = lri1,
      Lri          = Lri,
      Lambda       = Lambda,
      sqrtdetinvAri = sqrtdetinvAri,
      resAri       = resAri,
      Srii         = Srii,
      v1           = v1,
      v            = v,
      A11          = A11,
      A12          = A12,
      A21          = A21,
      A22          = A22,
      resAi        = resAi
    ))
  }
}
