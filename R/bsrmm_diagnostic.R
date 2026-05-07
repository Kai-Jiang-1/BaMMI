#' MCMC Diagnostics for Bayesian Sparse Regression
#'
#' Produces diagnostic plots and statistics for multiple chains from
#' \code{\link{bsrmmgibbs}}, including trace plots, density plots,
#' and Gelman-Rubin Rhat convergence diagnostics for sigma2, theta,
#' and selected beta coefficients.
#'
#' @importFrom stats density
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
#' @importFrom grDevices dev.new graphics.off adjustcolor
#' @importFrom graphics abline axis barplot layout legend lines mtext par plot.new points text
#'
#' @param bsrmm_list List. A list of at least 2 outputs from
#'   \code{\link{bsrmmgibbs}}, each representing one chain.
#' @param nburnin Integer. Number of burn-in iterations used in
#'   \code{\link{bsrmmgibbs}}.
#' @param niter Integer. Number of posterior sampling iterations used
#'   in \code{\link{bsrmmgibbs}}.
#' @param threshold Numeric. PPI threshold for feature selection used
#'   to determine which beta coefficients to include in Rhat.
#'   Default is 0.5.
#'
#' @return A list containing:
#' \describe{
#'   \item{rhat_sigma2}{Rhat value for sigma2.}
#'   \item{rhat_theta}{Rhat value for theta.}
#'   \item{rhat_beta}{Named numeric vector of Rhat values for each
#'     selected beta coefficient.}
#'   \item{ESS_sigma2}{Effective sample size for sigma2.}
#'   \item{ESS_theta}{Effective sample size for theta.}
#'   \item{selected}{Integer vector of selected feature indices.}
#' }
#'
#' @seealso \code{\link{bsrmmgibbs}}, \code{\link{bsrmm_forest}},
#'   \code{\link{bsrmm_ppi}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' mydata <- gen_missing_value(
#'   n = 100, p = 200, snr = 1,
#'   proportion = 0.30,
#'   type = "mnar0.33"
#' )
#' predictor <- mydata$X
#' outcome   <- mydata$Y_miss
#' stand     <- list()
#' stand$mux <- colMeans(predictor)
#' stand$Sx  <- apply(predictor, 2, sd)
#' predictor <- apply(predictor, 2, function(x) (x - mean(x)) / sd(x))
#' set.seed(33)
#' train_ind <- sort(sample(nrow(predictor),
#'                    size = round(0.7 * nrow(predictor)),
#'                    replace = FALSE))
#' X <- as.matrix(predictor[train_ind, ])
#' Y <- as.matrix(outcome[train_ind, ])
#' n <- nrow(X); p <- ncol(X); c <- 100
#' N <- rbind(diag(p), matrix(c, 1, p))
#' Q <- mydata$Q; a <- rep(-12, p); nop <- floor(n / 2)
#' nburnin <- 1000; niter <- 2000
#'
#' ## Run multiple chains
#' chain1 <- bsrmmgibbs(nburnin=nburnin, niter=niter, p=p, nop=nop,
#'   Y=Y, X=X, N=N, a=a, Q=Q, n=n, display=FALSE)
#' chain2 <- bsrmmgibbs(nburnin=nburnin, niter=niter, p=p, nop=nop,
#'   Y=Y, X=X, N=N, a=a, Q=Q, n=n, display=FALSE)
#'
#' diag_result <- bsrmm_diagnostic(
#'   bsrmm_list = list(chain1, chain2),
#'   nburnin    = nburnin,
#'   niter      = niter
#' )
#' }
bsrmm_diagnostic <- function(bsrmm_list, nburnin, niter,
                             threshold = 0.5) {

  ## Input validation ----
  if (!is.list(bsrmm_list)) {
    stop("bsrmm_list must be a list of bsrmmgibbs outputs.")
  }
  if (length(bsrmm_list) < 2) {
    stop("bsrmm_list must contain at least 2 chains for diagnostics.")
  }

  n_chains <- length(bsrmm_list)
  post_idx <- (nburnin + 2):(nburnin + niter + 1)

  ## Chain colors ----
  chain_cols <- c("#E53935", "#1E88E5", "#43A047",
                  "#FB8C00", "#8E24AA", "#00ACC1",
                  "#F4511E", "#3949AB")
  chain_cols <- rep(chain_cols, length.out = n_chains)

  ## Identify selected features across all chains ----
  PPI_all  <- lapply(bsrmm_list, function(x)
    colMeans(x$gamma[post_idx, ]))
  PPI_mean <- Reduce("+", PPI_all) / n_chains
  selected <- which(PPI_mean > threshold)

  cat("=== MCMC Diagnostics ===\n")
  cat(sprintf("Number of chains  : %d\n", n_chains))
  cat(sprintf("Burn-in           : %d\n", nburnin))
  cat(sprintf("Iterations        : %d\n", niter))
  cat(sprintf("Selected features : %s\n\n",
              ifelse(length(selected) == 0, "None",
                     paste(selected, collapse = ", "))))

  ## Save and restore graphics state on exit ----
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  # Reset graphics state
  graphics.off()
  dev.new()

  n_rhat  <- 2 + length(selected)
  layout(matrix(c(1, 2, 3,
                  4, 5, 6,
                  7, 7, 7),
                nrow = 3, byrow = TRUE))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 4, 0))

  ## Layout ----
  # Row 1: trace plots     (sigma2, theta)
  # Row 2: density plots   (sigma2, theta)
  # Row 3: chain legend    (shared, full width)
  # Row 4: Rhat plot       (full width)
  # Row 5: Rhat legend     (full width)
  layout(matrix(c(1, 2,
                  3, 4,
                  5, 5,
                  6, 6,
                  7, 7),
                nrow = 5, byrow = TRUE),
         heights = c(4, 4, 0.8, 4, 0.8))

  par(oma = c(0, 0, 4, 0))

  ## --- Trace plots ---

  # sigma2
  sigma2_chains <- lapply(bsrmm_list, function(x) x$sigma2[post_idx])
  y_rng <- range(unlist(sigma2_chains))
  par(mar = c(4, 4, 3, 1))
  plot(sigma2_chains[[1]], type = "l",
       col = chain_cols[1], lwd = 1,
       ylim = y_rng,
       xlab = "Iteration",
       ylab = expression(sigma^2),
       main = expression("Trace: " * sigma^2))
  for (k in 2:n_chains)
    lines(sigma2_chains[[k]], col = chain_cols[k], lwd = 1)

  # theta
  theta_chains <- lapply(bsrmm_list, function(x) x$theta[post_idx])
  y_rng <- range(unlist(theta_chains))
  par(mar = c(4, 4, 3, 1))
  plot(theta_chains[[1]], type = "l",
       col = chain_cols[1], lwd = 1,
       ylim = y_rng,
       xlab = "Iteration",
       ylab = expression(theta),
       main = expression("Trace: " * theta))
  for (k in 2:n_chains)
    lines(theta_chains[[k]], col = chain_cols[k], lwd = 1)

  ## --- Density plots ---

  # sigma2
  dens_list <- lapply(sigma2_chains, density)
  x_rng     <- range(sapply(dens_list, function(d) range(d$x)))
  y_rng     <- range(sapply(dens_list, function(d) range(d$y)))
  par(mar = c(4, 4, 3, 1))
  plot(dens_list[[1]], col = chain_cols[1], lwd = 2,
       xlim = x_rng, ylim = y_rng,
       main = expression("Density: " * sigma^2),
       xlab = expression(sigma^2))
  for (k in 2:n_chains)
    lines(dens_list[[k]], col = chain_cols[k], lwd = 2)

  # theta
  dens_list <- lapply(theta_chains, density)
  x_rng     <- range(sapply(dens_list, function(d) range(d$x)))
  y_rng     <- range(sapply(dens_list, function(d) range(d$y)))
  par(mar = c(4, 4, 3, 1))
  plot(dens_list[[1]], col = chain_cols[1], lwd = 2,
       xlim = x_rng, ylim = y_rng,
       main = expression("Density: " * theta),
       xlab = expression(theta))
  for (k in 2:n_chains)
    lines(dens_list[[k]], col = chain_cols[k], lwd = 2)

  ## --- Chain legend panel (row 3) ---
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = paste0("Chain ", 1:n_chains),
         col    = chain_cols[1:n_chains],
         lwd    = 2,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.85)

  ## --- Rhat plot (row 4) ---

  # sigma2
  rhat_sigma2 <- coda::gelman.diag(
    coda::mcmc.list(lapply(sigma2_chains, coda::mcmc))
  )$psrf[1, 1]

  # theta
  rhat_theta <- coda::gelman.diag(
    coda::mcmc.list(lapply(theta_chains, coda::mcmc))
  )$psrf[1, 1]

  # selected betas
  rhat_beta <- numeric(0)
  if (length(selected) > 0) {
    rhat_beta <- sapply(selected, function(j) {
      beta_chains <- lapply(bsrmm_list, function(x) {
        vals <- x$tembeta_save[j, post_idx]
        vals[is.na(vals)] <- 0
        coda::mcmc(vals)
      })
      coda::gelman.diag(
        coda::mcmc.list(beta_chains)
      )$psrf[1, 1]
    })
    names(rhat_beta) <- paste0("beta_", selected)
  }

  all_rhat  <- c(sigma2 = rhat_sigma2,
                 theta  = rhat_theta,
                 rhat_beta)
  rhat_cols <- ifelse(all_rhat < 1.1, "#43A047",
                      ifelse(all_rhat < 1.2, "#FB8C00", "#E53935"))

  par(mar = c(6, 5, 3, 2))
  barplot(all_rhat,
          col       = rhat_cols,
          border    = NA,
          ylim      = c(0, max(all_rhat, 1.3) * 1.1),
          ylab      = "Rhat",
          main      = "Gelman-Rubin Rhat Diagnostic",
          las       = 2,
          cex.names = 0.8)
  abline(h = 1.1, col = "#FB8C00", lty = 2, lwd = 2)
  abline(h = 1.2, col = "#E53935", lty = 2, lwd = 2)
  abline(h = 1.0, col = "gray50",  lty = 1, lwd = 1)

  ## --- Rhat legend panel (row 5) ---
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = c("Rhat < 1.1 (converged)",
                    "1.1 <= Rhat < 1.2 (warning)",
                    "Rhat >= 1.2 (not converged)"),
         fill   = c("#43A047", "#FB8C00", "#E53935"),
         border = NA,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.8)

  ## Overall title ----
  mtext("MCMC Diagnostics", outer = TRUE, cex = 1.3, font = 2)

  ## ESS ----
  ESS_sigma2 <- coda::effectiveSize(
    coda::mcmc.list(lapply(sigma2_chains, coda::mcmc)))
  ESS_theta  <- coda::effectiveSize(
    coda::mcmc.list(lapply(theta_chains, coda::mcmc)))

  ## Print summary ----
  cat("=== Rhat Summary ===\n")
  cat(sprintf("Rhat sigma2 : %.3f %s\n", rhat_sigma2,
              ifelse(rhat_sigma2 < 1.1, "[OK]",
                     ifelse(rhat_sigma2 < 1.2, "[WARNING]",
                            "[NOT CONVERGED]"))))
  cat(sprintf("Rhat theta  : %.3f %s\n", rhat_theta,
              ifelse(rhat_theta < 1.1, "[OK]",
                     ifelse(rhat_theta < 1.2, "[WARNING]",
                            "[NOT CONVERGED]"))))
  if (length(rhat_beta) > 0) {
    cat("Rhat beta:\n")
    for (nm in names(rhat_beta)) {
      cat(sprintf("  %-12s: %.3f %s\n", nm, rhat_beta[nm],
                  ifelse(rhat_beta[nm] < 1.1, "[OK]",
                         ifelse(rhat_beta[nm] < 1.2, "[WARNING]",
                                "[NOT CONVERGED]"))))
    }
  }
  cat("\n=== ESS Summary ===\n")
  cat(sprintf("ESS sigma2 : %.1f\n", ESS_sigma2))
  cat(sprintf("ESS theta  : %.1f\n", ESS_theta))

  ## Return ----
  return(invisible(list(
    rhat_sigma2 = rhat_sigma2,
    rhat_theta  = rhat_theta,
    rhat_beta   = rhat_beta,
    ESS_sigma2  = ESS_sigma2,
    ESS_theta   = ESS_theta,
    selected    = selected
  )))
}
