#' Posterior Probability of Inclusion Plot
#'
#' Plots posterior probabilities of inclusion (PPI) from a BSRMM fit.
#' By default only the selected features are shown, but any subset of
#' features (or all of them) can be plotted via \code{features}. Up to
#' \code{max_bar} features are drawn as a labelled bar plot; larger sets
#' are drawn as a needle plot against the feature index.
#'
#' @importFrom graphics abline axis barplot layout legend par plot
#' @importFrom graphics plot.new points strwidth text
#'
#' @param bsrmm List. Output from \code{\link{bsrmmgibbs}}.
#' @param nburnin Integer. Number of burn-in iterations.
#' @param niter Integer. Number of posterior sampling iterations.
#' @param threshold Numeric. PPI threshold for feature selection.
#'   Default is 0.5.
#' @param features Which features to plot: \code{"selected"} (default,
#'   features with PPI > \code{threshold}), \code{"all"}, or an integer
#'   vector of feature indices, e.g. \code{1:100}.
#' @param sort Logical. If \code{TRUE}, features are ordered by decreasing
#'   PPI; if \code{FALSE}, by feature index. Default \code{NULL} sorts only
#'   when \code{features = "selected"}.
#' @param true_features Optional integer vector of true feature indices
#'   (e.g. from simulated data). These are marked with a triangle on the
#'   plot and used to report true/false positives. Default \code{NULL}.
#' @param max_bar Integer. Maximum number of features drawn as a labelled
#'   bar plot; above this a needle plot is used. Default is 50.
#'
#' @return Invisibly, a list containing:
#' \describe{
#'   \item{PPI}{Numeric vector of posterior probabilities of inclusion
#'     for all features.}
#'   \item{selected}{Integer vector of selected feature indices.}
#'   \item{plotted}{Integer vector of plotted feature indices, in plot
#'     order.}
#' }
#'
#' @seealso \code{\link{bsrmmgibbs}}, \code{\link{bsrmm_diagnostic}},
#'   \code{\link{bsrmm_forest}}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' mydata <- gen_missing_value(
#'   n = 200, p = 200, snr = 1,
#'   proportion = 0.30,
#'   mnar_ratio = 1/3
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
#'
#' ## Selected features only (default)
#' ppi_result <- bsrmm_ppi(bsrmm, nburnin, niter)
#'
#' ## All features, marking the true ones
#' true_idx <- which(mydata$coeffs[, 1] == 1)
#' bsrmm_ppi(bsrmm, nburnin, niter, features = "all",
#'           true_features = true_idx)
#'
#' ## First 100 features only
#' bsrmm_ppi(bsrmm, nburnin, niter, features = 1:100,
#'           true_features = true_idx)
#' }
bsrmm_ppi <- function(bsrmm, nburnin, niter,
                      threshold     = 0.5,
                      features      = "selected",
                      sort          = NULL,
                      true_features = NULL,
                      max_bar       = 50) {

  ## Compute PPI ----
  post_idx <- (nburnin + 2):(nburnin + niter + 1)
  PPI      <- colMeans(bsrmm$gamma[post_idx, ])
  p        <- length(PPI)
  selected <- which(PPI > threshold)

  ## Resolve which features to plot ----
  check_idx <- function(x, arg) {
    if (!is.numeric(x) || length(x) == 0 || anyNA(x) ||
        any(x != round(x)) || any(x < 1 | x > p)) {
      stop(sprintf("%s must be integer indices between 1 and %d.", arg, p))
    }
    unique(as.integer(x))
  }

  if (is.character(features)) {
    features <- match.arg(features, c("selected", "all"))
    mode     <- features
    idx      <- if (mode == "all") seq_len(p) else selected
  } else {
    mode <- "custom"
    idx  <- check_idx(features, "features")
  }
  if (!is.null(true_features)) {
    true_features <- check_idx(true_features, "true_features")
  }
  if (is.null(sort)) sort <- (mode == "selected")

  if (length(idx) == 0) {
    message("No features selected at threshold = ", threshold)
    return(invisible(list(PPI = PPI, selected = selected,
                          plotted = integer(0))))
  }

  idx     <- if (sort) idx[order(PPI[idx], decreasing = TRUE)] else sort(idx)
  ppi     <- PPI[idx]
  is_sel  <- ppi > threshold
  is_true <- idx %in% true_features
  m       <- length(idx)

  ## Colours ----
  col_sel   <- "#E53935"
  col_unsel <- "#BDBDBD"
  col_true  <- "#1E88E5"
  cols      <- ifelse(is_sel, col_sel, col_unsel)

  main <- if (mode == "selected") {
    sprintf("PPI for Selected Features (threshold = %s)", threshold)
  } else {
    sprintf("PPI for %d of %d Features (threshold = %s)", m, p, threshold)
  }

  ## Save and restore graphics state on exit ----
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)

  ## Layout: main plot + legend row ----
  layout(matrix(c(1, 2), nrow = 2), heights = c(9, 1))

  if (m <= max_bar) {
    ## Bar plot with labels ----
    par(mar = c(6, 5, 4, 2))

    bp <- barplot(
      ppi,
      col       = cols,
      border    = NA,
      ylab      = "Posterior Probability of Inclusion (PPI)",
      main      = main,
      ylim      = c(0, 1.2),
      names.arg = rep("", m),
      xaxt      = "n"
    )

    axis(1,
         at       = bp,
         labels   = paste0("Feature ", idx),
         las      = 2,
         tcl      = -0.3,
         cex.axis = 0.8,
         col.axis = "#424242")

    text(x      = bp,
         y      = ppi + 0.04,
         labels = sprintf("%.2f", ppi),
         cex    = 0.75,
         col    = ifelse(is_sel, "#B71C1C", "#616161"),
         font   = 2)

    x_true <- bp[is_true]

  } else {
    ## Needle plot for many features ----
    par(mar = c(5, 5, 4, 2))

    x    <- if (sort) seq_len(m) else idx
    xlab <- if (sort) "Features (ranked by PPI)" else "Feature index"

    plot(x, ppi,
         type = "h",
         col  = cols,
         lwd  = if (m <= 200) 2 else 1,
         ylim = c(0, 1.2),
         xlab = xlab,
         ylab = "Posterior Probability of Inclusion (PPI)",
         main = main,
         yaxt = "n",
         las  = 1)
    axis(2, at = seq(0, 1, by = 0.2), las = 1)

    x_true <- x[is_true]
  }

  abline(h = threshold, lty = 2, col = "#757575")
  if (length(x_true) > 0) {
    points(x_true, rep(1.1, length(x_true)),
           pch = 17, col = col_true, cex = 0.8)
  }

  ## Legend panel (only categories present in the plot) ----
  leg <- data.frame(
    text = c(paste0("Selected (PPI > ", threshold, ")"),
             "Not selected", "True feature"),
    col  = c(col_sel, col_unsel, col_true),
    pch  = c(15, 15, 17),
    show = c(any(is_sel), any(!is_sel), any(is_true)),
    stringsAsFactors = FALSE
  )
  leg <- leg[leg$show, ]

  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend     = leg$text,
         col        = leg$col,
         pch        = leg$pch,
         pt.cex     = 1.5,
         horiz      = TRUE,
         bty        = "n",
         cex        = 0.9,
         text.width = strwidth(leg$text, cex = 0.9) + 0.03)

  ## Summary ----
  cat("=== PPI Summary ===\n")
  cat(sprintf("Total features    : %d\n", p))
  cat(sprintf("Plotted features  : %d\n", m))
  cat(sprintf("Selected features : %d\n", length(selected)))
  cat(sprintf("Selected indices  : %s\n",
              paste(selected, collapse = ", ")))
  if (!is.null(true_features)) {
    tp <- sum(selected %in% true_features)
    cat(sprintf("True features     : %d\n", length(true_features)))
    cat(sprintf("True positives    : %d\n", tp))
    cat(sprintf("False positives   : %d\n", length(selected) - tp))
    cat(sprintf("False negatives   : %d\n", length(true_features) - tp))
  }

  ## Return ----
  return(invisible(list(
    PPI      = PPI,
    selected = selected,
    plotted  = idx
  )))
}
