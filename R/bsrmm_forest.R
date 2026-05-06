#' Forest Plot for Selected Features
#'
#' Produces a forest plot of posterior mean estimates and 95 percent
#' credible intervals for features selected by \code{\link{bsrmmgibbs}}.
#'
#' @param bsrmm List. Output from \code{\link{bsrmmgibbs}}.
#' @param nburnin Integer. Number of burn-in iterations.
#' @param niter Integer. Number of posterior sampling iterations.
#' @param stand List. Standardization parameters from predictor
#'   preparation, containing \code{Sx} (column standard deviations).
#' @param threshold Numeric. PPI threshold for feature selection.
#'   Default is 0.5.
#'
#' @return A list containing:
#' \describe{
#'   \item{selected}{Integer vector of selected feature indices.}
#'   \item{beta_mean}{Numeric vector of posterior mean estimates.}
#'   \item{lower_cri}{Numeric vector of lower 95 percent credible
#'     interval bounds.}
#'   \item{upper_cri}{Numeric vector of upper 95 percent credible
#'     interval bounds.}
#' }
#'
#' @seealso \code{\link{bsrmmgibbs}}, \code{\link{bsrmm_diagnostic}},
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
#' bsrmm <- bsrmmgibbs(nburnin=nburnin, niter=niter, p=p, nop=nop,
#'   Y=Y, X=X, N=N, a=a, Q=Q, n=n, display=FALSE)
#' forest_result <- bsrmm_forest(bsrmm, nburnin, niter, stand)
#' }
bsrmm_forest <- function(bsrmm, nburnin, niter, stand,
                         threshold = 0.5) {

  ## Save and restore graphics state on exit ----
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)

  ## Extract posterior samples ----
  post_idx <- (nburnin + 2):(nburnin + niter + 1)
  PPI      <- colMeans(bsrmm$gamma[post_idx, ])
  selected <- which(PPI > threshold)

  if (length(selected) == 0) {
    message("No features selected at threshold = ", threshold)
    return(invisible(NULL))
  }

  ## Compute beta estimates ----
  beta_esti <- bsrmm$tembeta_save[selected, post_idx, drop = FALSE]
  beta_esti[is.na(beta_esti)] <- 0

  # Rescale by standard deviation
  beta_esti <- sweep(beta_esti, 1,
                     STATS = stand$Sx[selected], FUN = "/")

  beta_mean  <- rowMeans(beta_esti)
  lower_cri  <- apply(beta_esti, 1, function(x) quantile(x, 0.025))
  upper_cri  <- apply(beta_esti, 1, function(x) quantile(x, 0.975))
  feat_names <- paste0("Feature ", selected)

  ## Colors based on direction ----
  cols <- ifelse(beta_mean > 0, "#E53935", "#1E88E5")

  ## Layout: main plot + legend row ----
  layout(matrix(c(1, 2), nrow = 2), heights = c(8, 1))

  ## Plot setup ----
  n_feat  <- length(selected)
  y_pos   <- seq(n_feat, 1)
  x_range <- range(c(lower_cri, upper_cri, 0))
  x_pad   <- diff(x_range) * 0.15
  x_lim   <- c(x_range[1] - x_pad, x_range[2] + x_pad)

  par(mar = c(3, 8, 4, 2))

  ## Empty plot ----
  plot(NA, xlim = x_lim,
       ylim = c(0.5, n_feat + 0.5),
       xlab = "Posterior Mean (95% CrI)",
       ylab = "",
       main = "Forest Plot: Selected Features",
       yaxt = "n", bty = "l")

  ## Reference line at zero ----
  abline(v = 0, lty = 2, col = "gray50", lwd = 1.5)

  ## Grid lines ----
  for (y in y_pos) {
    abline(h = y, col = "gray90", lty = 1)
  }

  ## CrI lines ----
  for (k in seq_along(selected)) {
    lines(c(lower_cri[k], upper_cri[k]),
          c(y_pos[k], y_pos[k]),
          col = cols[k], lwd = 2)
    # Whisker ends
    lines(c(lower_cri[k], lower_cri[k]),
          c(y_pos[k] - 0.15, y_pos[k] + 0.15),
          col = cols[k], lwd = 2)
    lines(c(upper_cri[k], upper_cri[k]),
          c(y_pos[k] - 0.15, y_pos[k] + 0.15),
          col = cols[k], lwd = 2)
  }

  ## Point estimates ----
  points(beta_mean, y_pos,
         pch = 18, cex = 1.5, col = cols)

  ## Feature labels ----
  axis(2, at = y_pos, labels = feat_names,
       las = 2, cex.axis = 0.85)

  ## Legend panel ----
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = c("Positive", "Negative"),
         col    = c("#E53935", "#1E88E5"),
         pch    = 18,
         lwd    = 2,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.9)

  ## Return ----
  return(invisible(list(
    selected  = selected,
    beta_mean = beta_mean,
    lower_cri = lower_cri,
    upper_cri = upper_cri
  )))
}
