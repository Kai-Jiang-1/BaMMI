#' Gibbs Sampler for Bayesian Sparse Regression with Missing Value Imputation
#'
#' Implements a Gibbs sampler for Bayesian sparse regression for
#' microbiome-metabolite data integration (BSRMM). The sampler performs
#' simultaneous variable selection and missing value imputation using a
#' Bayesian Adaptive Zone Estimation (BAZE) approach for variable selection
#' and truncated normal sampling for missing value imputation.
#'
#' @importFrom stats rnorm rbeta rbinom pnorm
#' @importFrom truncnorm rtruncnorm
#' @importFrom MCMCpack rinvgamma
#' @importFrom mvnfast rmvt
#'
#' @param nburnin Integer. Number of burn-in iterations.
#' @param niter Integer. Number of posterior sampling iterations after
#'   burn-in.
#' @param p Integer. Number of predictors (microbiome features).
#' @param nop Integer. Initial number of predictors included in the model.
#'   If \code{nop > p}, it is reset to \code{floor(p/2)}.
#' @param Y Numeric matrix of dimension n x 1. Outcome vector with
#'   missing values set to 0.
#' @param X Numeric matrix of dimension n x p. Standardized predictor
#'   matrix of microbiome features.
#' @param N Numeric matrix of dimension (p+1) x p. Constraint matrix
#'   constructed as \code{rbind(diag(p), rep(c, p))} where \code{c} is
#'   a large constant to maintain compositionality.
#' @param a Numeric vector of length p. Log prior inclusion probabilities
#'   for each predictor. Typically initialized as \code{rep(a0, p)} where
#'   \code{a0} is a negative constant (e.g., \code{-12}).
#' @param Q Numeric matrix of dimension p x p. Prior correlation matrix
#'   encoding prior knowledge about predictor correlations. Use
#'   \code{diag(0, p)} for independent predictors.
#' @param n Integer. Number of observations.
#' @param display Logical. If \code{TRUE}, prints progress every 5000
#'   iterations. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma}{Numeric matrix of dimension (nburnin + niter + 1) x p.
#'     Binary variable selection indicators at each iteration. Use
#'     \code{colMeans(gamma[(nburnin+2):(nburnin+niter+1), ])} to compute
#'     posterior probabilities of inclusion (PPI).}
#'   \item{nselect}{Numeric vector of length (nburnin + niter + 1).
#'     Number of selected variables at each iteration.}
#'   \item{theta}{Numeric vector of length (nburnin + niter + 1).
#'     Posterior samples of the probability of MAR.}
#'   \item{phi}{Numeric matrix. Posterior probability of each missing
#'     value being below the LOD at each iteration.}
#'   \item{z}{Numeric matrix. Binary indicators of whether each missing
#'     value is below the LOD (1) or above (0) at each iteration.}
#'   \item{missing_index}{Numeric matrix. Missing value classification
#'     at each iteration: 0 = MAR (above LOD), 1 = MNAR (below LOD),
#'     2 = observed.}
#'   \item{Y_mvi}{Numeric matrix of dimension n x (nburnin + niter + 1).
#'     Imputed outcome values at each iteration.}
#'   \item{Ri}{Numeric vector of length n. Missingness indicator:
#'     0 = missing, 1 = observed.}
#'   \item{sigma2}{Numeric vector of length (nburnin + niter + 1).
#'     Posterior samples of the residual variance.}
#'   \item{beta_0}{Numeric vector of length (nburnin + niter + 1).
#'     Posterior samples of the intercept.}
#'   \item{tembeta_save}{Numeric matrix of dimension p x
#'     (nburnin + niter + 1). Posterior samples of regression
#'     coefficients for all predictors.}
#' }
#'
#' @seealso \code{\link{bsrmmbf}}, \code{\link{gen_missing_value}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Generate simulation data
#' set.seed(123)
#' mydata <- gen_missing_value(
#'   n = 100, p = 200, snr = 1,
#'   proportion = 0.30,
#'   type = "mnar0.33"
#' )
#'
#' ## Prepare inputs
#' predictor <- mydata$X
#' outcome   <- mydata$Y_miss
#'
#' ## Standardize predictors
#' stand     <- list()
#' stand$mux <- colMeans(predictor)
#' stand$Sx  <- apply(predictor, 2, sd)
#' predictor <- apply(predictor, 2, function(x) (x - mean(x)) / sd(x))
#'
#'
#' ## Train/test split
#' set.seed(33)
#' train_ind       <- sort(sample(nrow(predictor),
#'                          size    = round(0.7 * nrow(predictor)),
#'                          replace = FALSE))
#' predictor_train <- predictor[train_ind, ]
#' outcome_train   <- outcome[train_ind, ]
#'
#' X <- as.matrix(predictor_train)
#' Y <- as.matrix(outcome_train)
#'
#' ## Set up model parameters
#' n   <- nrow(X)
#' p   <- ncol(X)
#' c   <- 100
#' N   <- rbind(diag(x = 1, nrow = p, ncol = p),
#'              matrix(data = c, nrow = 1, ncol = p))
#' Q   <- mydata$Q
#' a0  <- -12
#' a   <- rep(a0, p)
#' nop <- floor(n / 2)
#'
#' nburnin <- 5000
#' niter   <- 10000
#'
#' ## Run Gibbs sampler
#' bsrmm <- bsrmmgibbs(
#'   nburnin = nburnin,
#'   niter   = niter,
#'   p       = p,
#'   nop     = nop,
#'   Y       = Y,
#'   X       = X,
#'   N       = N,
#'   a       = a,
#'   Q       = Q,
#'   n       = n,
#'   display = FALSE
#' )
#'
#' ## Posterior probabilities of inclusion
#' PPI      <- colMeans(bsrmm$gamma[(nburnin+2):(nburnin+niter+1), ])
#' selected <- which(PPI > 0.5)
#' cat("Selected features (PPI > 0.5):", selected, "\n")
#' cat("True features:", which(mydata$coeffs[, 1] == 1), "\n")
#'
#' ## Confusion table
#' tab <- table(
#'   Predicted = PPI > 0.5,
#'   Truth     = mydata$coeffs[, 1]
#' )
#' print(tab)}
bsrmmgibbs <- function(nburnin, niter, p, nop, Y, X, N, a, Q, n,
                       display = FALSE) {

  ## Fixed hyperparameters ----
  tau           <- 1
  nu            <- 0
  omega         <- 0
  theta_a_shape <- 1
  theta_b_shape <- 1

  ## Initialize ----
  if (nop > p) nop <- floor(p / 2)

  index           <- sort(sample(1:p, nop, replace = FALSE))
  nop             <- length(index)
  gamma           <- matrix(0, nrow = nburnin + niter + 1, ncol = p)
  gamma[1, index] <- 1
  sigma2          <- numeric(nburnin + niter + 1)
  theta           <- numeric(nburnin + niter + 1)
  theta[1]        <- rbeta(1, shape1 = theta_a_shape, shape2 = theta_b_shape)

  ## Initial precision matrix and cached quantities ----
  Xri    <- X[, index]
  Tri    <- N[, index]
  Lambda <- diag(1, nrow = n, ncol = n) -
    matrix(1, nrow = n, ncol = n) / (n + 1)
  Ari    <- t(Xri) %*% Lambda %*% Xri +
    tau^(-2) * (t(Tri) %*% Tri)
  invAri <- chol2inv(chol(Ari))
  Lri    <- t(chol(invAri))

  keep              <- list()
  keep$resAi        <- t(Y) %*% Lambda %*% Y -
    t(Y) %*% Lambda %*% X[, index] %*% invAri %*%
    t(X[, index]) %*% Lambda %*% Y
  keep$sqrtdetinvAi <- sum(log(diag(Lri)))

  ## Initial beta, sigma2, beta_0 ----
  tembeta <- as.matrix(t(rmvt(
    n     = 1,
    mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
    sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
    df    = n + nu
  )))

  tembeta_save           <- matrix(NA, nrow = ncol(X), ncol = nburnin + niter + 1)
  tembeta_save[index, 1] <- tembeta

  sigma2[1] <- rinvgamma(
    n     = 1,
    shape = 0.5 * (n + nu),
    scale = 0.5 * (as.numeric(keep$resAi) + nu * omega)
  )

  beta_0    <- numeric(nburnin + niter + 1)
  beta_0[1] <- rnorm(
    1,
    mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
    sd   = sqrt(sigma2[1] / (n + 1))
  )

  ## Missingness setup ----
  Ri         <- rep(0, nrow(Y))
  Ri[Y != 0] <- 1
  xi         <- sort(unique(Y[Y != 0]))[1]

  ## Initial missing value imputation ----
  phi <- matrix(NA, nrow = n, ncol = nburnin + niter + 1)
  phi[which(Ri == 0), 1] <- pnorm(
    xi,
    mean       = beta_0[1] + X[which(Ri == 0), index] %*% tembeta,
    sd         = rep(sqrt(sigma2[1]), sum(Ri == 0)),
    lower.tail = TRUE
  )

  phi_values  <- phi[which(Ri == 0), 1]
  prob_values <- 1 / (1 + (theta[1] * (1 - phi_values) / phi_values))

  z <- matrix(NA, nrow = n, ncol = nburnin + niter + 1)
  z[which(Ri == 0), 1] <- rbinom(length(prob_values), size = 1, prob_values)

  missing_index <- matrix(NA, nrow = n, ncol = nburnin + niter + 1)
  missing_index[, 1] <- sapply(1:n, function(j) {
    if (Ri[j] == 1)         return(2)  # observed
    else if (z[j, 1] == 1) return(1)  # MNAR: below LOD
    else                    return(0)  # MAR: above LOD
  })

  y_mvi <- matrix(NA, nrow = n, ncol = nburnin + niter + 1)
  y_mvi[, 1] <- sapply(1:n, function(s) {
    if (missing_index[s, 1] == 1) {
      rtruncnorm(1, a = -Inf, b = xi,
                 mean = beta_0[1] + X[s, index] %*% tembeta,
                 sd   = sqrt(sigma2[1]))
    } else if (missing_index[s, 1] == 0) {
      rtruncnorm(1, a = xi, b = Inf,
                 mean = beta_0[1] + X[s, index] %*% tembeta,
                 sd   = sqrt(sigma2[1]))
    } else {
      Y[s, 1]
    }
  })

  Y        <- as.matrix(y_mvi[, 1])
  theta_a  <- sum(missing_index[, 1] == 0) + 1
  theta_b  <- sum(missing_index[, 1] == 2) + 1
  theta[1] <- rbeta(1, theta_a, theta_b)

  nselect    <- numeric(nburnin + niter + 1)
  nselect[1] <- nop

  if (display) {
    k <- 1
    cat("Gibbs Sampling starting - progress printed every 5000 iterations.\n")
  }

  ## Main Gibbs loop ----
  for (i in 1:(nburnin + niter)) {

    proposeindx    <- sample(1:p, 1, replace = TRUE)
    gamma[i + 1, ] <- gamma[i, ]
    flag           <- any(index == proposeindx) & (length(index) > 1)

    if (flag) {
      ## Exclusion step ----
      indxtemp   <- index[index != proposeindx]
      Xri        <- X[, indxtemp]
      Xi         <- X[, proposeindx]
      XIi        <- cbind(Xri, Xi)
      Tri        <- N[, indxtemp]
      Ti         <- N[, proposeindx]
      TIi        <- cbind(Tri, Ti)
      invAritemp <- invAri
      tn         <- length(index)
      seq_idx    <- 1:tn

      if (proposeindx == max(index)) {
        idx <- list(seq_idx, seq_idx)
      } else if (proposeindx == min(index)) {
        idx <- list(c(2:tn, 1), c(2:tn, 1))
      } else {
        ti  <- seq_idx[index == proposeindx]
        idx <- list(c(1:(ti-1), (ti+1):tn, ti),
                    c(1:(ti-1), (ti+1):tn, ti))
      }

      invAitemp <- invAri[idx[[1]], idx[[2]]]
      invAri1   <- invAitemp[1:(tn-1), 1:(tn-1)]
      invAri2   <- invAitemp[1:(tn-1), tn, drop = FALSE]
      invAri3   <- invAitemp[tn, tn]
      invAri    <- invAri1 - invAri2 %*% t(invAri2) / invAri3

      result1  <- bsrmmbf(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri,
                          n, tau, nu, omega, flag, keep)
      BF1      <- result1$BF1
      keep     <- result1$keep

      pgammai1 <- exp(a[proposeindx] +
                        Q[proposeindx, indxtemp] %*% gamma[i, indxtemp])
      pcond    <- 1 / (1 + 1 / (BF1 * pgammai1))
      newgamma <- rbinom(1, 1, pcond)
      gamma[i + 1, proposeindx] <- newgamma

      if ((newgamma == 0) && (tn > 1)) {
        index             <- indxtemp
        keep$resAi        <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        tembeta <- as.matrix(t(rmvt(
          n     = 1,
          mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
          sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
          df    = n + nu
        )))
        tembeta_save[index, i + 1] <- tembeta
        sigma2[i + 1] <- rinvgamma(1, shape = 0.5*(n+nu),
                                   scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        beta_0[i + 1] <- rnorm(1,
                               mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
                               sd   = sqrt(sigma2[i+1] / (n + 1)))
        nselect[i + 1] <- tn - 1

      } else {
        invAri         <- invAritemp
        nselect[i + 1] <- nselect[i]
        tembeta <- as.matrix(t(rmvt(
          n     = 1,
          mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
          sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
          df    = n + nu
        )))
        tembeta_save[index, i + 1] <- tembeta
        sigma2[i + 1] <- rinvgamma(1, shape = 0.5*(n+nu),
                                   scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        beta_0[i + 1] <- rnorm(1,
                               mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
                               sd   = sqrt(sigma2[i+1] / (n + 1)))
      }

    } else if (any(index != proposeindx)) {
      ## Inclusion step ----
      indxtemp    <- index
      Xri         <- X[, indxtemp]
      Xi          <- X[, proposeindx]
      XIi         <- cbind(Xri, Xi)
      Tri         <- N[, indxtemp]
      Ti          <- N[, proposeindx]
      TIi         <- cbind(Tri, Ti)
      BayesResult <- bsrmmbf(Y, Xri, Xi, XIi, Tri, Ti, TIi, invAri,
                             n, tau, nu, omega, flag, keep)
      BF    <- BayesResult$BF
      keep  <- BayesResult$keep
      invAi <- BayesResult$invAi

      pgammai1 <- exp(a[proposeindx] +
                        Q[proposeindx, indxtemp] %*% gamma[i, indxtemp])
      pcond    <- 1 / (1 + 1 / (BF * pgammai1))
      newgamma <- rbinom(1, 1, pcond)
      gamma[i + 1, proposeindx] <- newgamma

      if (newgamma == 1) {
        index   <- sort(c(indxtemp, proposeindx))
        invAri  <- invAi
        tn      <- length(index)
        seq_idx <- 1:tn

        if (proposeindx > max(indxtemp)) {
          idx <- list(seq_idx, seq_idx)
        } else if (proposeindx < min(indxtemp)) {
          idx <- list(c(tn, 1:(tn-1)), c(tn, 1:(tn-1)))
        } else {
          ti  <- seq_idx[index == proposeindx]
          idx <- list(c(1:(ti-1), tn, ti:(tn-1)),
                      c(1:(ti-1), tn, ti:(tn-1)))
        }
        invAri <- invAri[idx[[1]], idx[[2]]]

        tembeta <- as.matrix(t(rmvt(
          n     = 1,
          mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
          sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
          df    = n + nu
        )))
        tembeta_save[index, i + 1] <- tembeta
        sigma2[i + 1] <- rinvgamma(1, shape = 0.5*(n+nu),
                                   scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        beta_0[i + 1] <- rnorm(1,
                               mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
                               sd   = sqrt(sigma2[i+1] / (n + 1)))
        nselect[i + 1] <- tn

      } else {
        keep$resAi        <- keep$resAri
        keep$sqrtdetinvAi <- keep$sqrtdetinvAri
        tembeta <- as.matrix(t(rmvt(
          n     = 1,
          mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
          sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
          df    = n + nu
        )))
        tembeta_save[index, i + 1] <- tembeta
        sigma2[i + 1] <- rinvgamma(1, shape = 0.5*(n+nu),
                                   scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
        beta_0[i + 1] <- rnorm(1,
                               mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
                               sd   = sqrt(sigma2[i+1] / (n + 1)))
        nselect[i + 1] <- nselect[i]
      }

    } else {
      ## No change ----
      keep$resAi        <- keep$resAri
      keep$sqrtdetinvAi <- keep$sqrtdetinvAri
      tembeta <- as.matrix(t(rmvt(
        n     = 1,
        mu    = invAri %*% t(X[, index]) %*% Lambda %*% Y,
        sigma = (as.numeric(keep$resAi) + nu * omega) / (n + nu) * invAri,
        df    = n + nu
      )))
      tembeta_save[index, i + 1] <- tembeta
      sigma2[i + 1] <- rinvgamma(1, shape = 0.5*(n+nu),
                                 scale = 0.5*(as.numeric(keep$resAi)+nu*omega))
      beta_0[i + 1] <- rnorm(1,
                             mean = sum(Y - X[, index] %*% tembeta) / (n + 1),
                             sd   = sqrt(sigma2[i+1] / (n + 1)))
      nselect[i + 1] <- nselect[i]
    }

    ## Missing value imputation ----
    phi[which(Ri == 0), i + 1] <- pnorm(
      xi,
      mean       = beta_0[i+1] + X[which(Ri == 0), index] %*% tembeta,
      sd         = rep(sqrt(sigma2[i+1]), sum(Ri == 0)),
      lower.tail = TRUE
    )

    phi_values  <- phi[which(Ri == 0), i + 1]
    prob_values <- 1 / (1 + (theta[i] * (1 - phi_values) / phi_values))
    z[which(Ri == 0), i + 1] <- rbinom(length(prob_values), 1, prob_values)

    missing_index[, i + 1] <- sapply(1:n, function(j) {
      if (Ri[j] == 1)             return(2)
      else if (z[j, i + 1] == 1) return(1)
      else                        return(0)
    })

    y_mvi[, i + 1] <- sapply(1:n, function(s) {
      if (missing_index[s, i + 1] == 1) {
        rtruncnorm(1, a = -Inf, b = xi,
                   mean = beta_0[i+1] + X[s, index] %*% tembeta,
                   sd   = sqrt(sigma2[i+1]))
      } else if (missing_index[s, i + 1] == 0) {
        rtruncnorm(1, a = xi, b = Inf,
                   mean = beta_0[i+1] + X[s, index] %*% tembeta,
                   sd   = sqrt(sigma2[i+1]))
      } else {
        y_mvi[s, i]
      }
    })

    Y            <- as.matrix(y_mvi[, i + 1])
    theta_a      <- sum(missing_index[, i + 1] == 0) + theta_a_shape
    theta_b      <- sum(missing_index[, i + 1] == 2) + theta_b_shape
    theta[i + 1] <- rbeta(1, theta_a, theta_b)

    if (display) {
      if (k %% 5000 == 0) cat("Iteration:", k, "\n")
      k <- k + 1
    }
  }

  if (display) cat("Gibbs Sampling complete.\n")

  ## Return ----
  return(list(
    gamma         = gamma,
    nselect       = nselect,
    theta         = theta,
    phi           = phi,
    z             = z,
    missing_index = missing_index,
    Y_mvi         = y_mvi,
    Ri            = Ri,
    sigma2        = sigma2,
    beta_0        = beta_0,
    tembeta_save  = tembeta_save
  ))
}
