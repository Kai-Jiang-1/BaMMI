#' Posterior Probability of Inclusion Plot
#'
#' Produces a bar plot of posterior probabilities of inclusion (PPI)
#' for selected features only, with vertical feature index labels.
#'
#' @param bsrmm List. Output from \code{\link{bsrmmgibbs}}.
#' @param nburnin Integer. Number of burn-in iterations.
#' @param niter Integer. Number of posterior sampling iterations.
#' @param threshold Numeric. PPI threshold for feature selection.
#'   Default is 0.5.
#'
#' @return A list containing:
#' \describe{
#'   \item{PPI}{Numeric vector of posterior probabilities of inclusion
#'     for all features.}
#'   \item{selected}{Integer vector of selected feature indices.}
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
#' ppi_result <- bsrmm_ppi(bsrmm, nburnin, niter)
#' }
bsrmm_ppi <- function(bsrmm, nburnin, niter,
                      threshold = 0.5) {

  ## Save and restore graphics state on exit ----
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)

  ## Compute PPI ----
  post_idx <- (nburnin + 2):(nburnin + niter + 1)
  PPI      <- colMeans(bsrmm$gamma[post_idx, ])
  p        <- length(PPI)
  selected <- which(PPI > threshold)

  ## Check if any features selected ----
  if (length(selected) == 0) {
    message("No features selected at threshold = ", threshold)
    return(invisible(list(PPI = PPI, selected = selected)))
  }

  ## Extract PPI for selected features only ----
  PPI_selected <- PPI[selected]
  feat_labels  <- paste0("Feature ", selected)

  ## Sort by PPI descending ----
  ord          <- order(PPI_selected, decreasing = TRUE)
  PPI_selected <- PPI_selected[ord]
  feat_labels  <- feat_labels[ord]
  selected_ord <- selected[ord]

  ## Layout: main plot + legend row ----
  layout(matrix(c(1, 2), nrow = 2), heights = c(9, 1))

  ## Main plot ----
  par(mar = c(6, 5, 4, 2))

  bp <- barplot(
    PPI_selected,
    col       = "#E53935",
    border    = NA,
    ylab      = "Posterior Probability of Inclusion (PPI)",
    main      = sprintf("PPI for Selected Features (threshold = %s)",
                        threshold),
    ylim      = c(0, 1.15),
    names.arg = rep("", length(selected)),
    xaxt      = "n"
  )

  ## Vertical x-axis labels ----
  axis(1,
       at       = bp,
       labels   = feat_labels,
       las      = 2,
       tcl      = -0.3,
       cex.axis = 0.8,
       col.axis = "#B71C1C")

  ## PPI value labels on top of bars ----
  text(x      = bp,
       y      = PPI_selected + 0.04,
       labels = sprintf("%.2f", PPI_selected),
       cex    = 0.75,
       col    = "#B71C1C",
       font   = 2)

  ## Legend panel ----
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = paste0("Selected features (PPI > ", threshold, ")"),
         fill   = "#E53935",
         border = NA,
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.9)

  ## Summary ----
  cat("=== PPI Summary ===\n")
  cat(sprintf("Total features    : %d\n", p))
  cat(sprintf("Selected features : %d\n", length(selected)))
  cat(sprintf("Selected indices  : %s\n",
              paste(selected, collapse = ", ")))

  ## Return ----
  return(invisible(list(
    PPI      = PPI,
    selected = selected
  )))
}
