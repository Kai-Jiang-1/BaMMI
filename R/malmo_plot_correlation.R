#' Correlation Heatmaps for MALMO
#'
#' Produces two heatmaps: (1) outcome vs outcome correlation from
#' the observed data, and (2) outcome vs predictor correlation
#' from the observed data.
#'
#' @importFrom stats cor
#' @importFrom graphics par layout image axis text plot.new legend mtext
#' @importFrom grDevices colorRampPalette
#'
#' @param Y Numeric matrix of dimension n x q. Outcome matrix.
#'   Can use \code{mydata$Y_true} or \code{mydata$Y_miss}.
#' @param X Numeric matrix of dimension n x p. Predictor matrix.
#'   Can use \code{mydata$X}.
#'
#' @return A list containing:
#' \describe{
#'   \item{cor_Y}{Numeric matrix of dimension q x q. Outcome vs
#'     outcome correlation matrix.}
#'   \item{cor_YX}{Numeric matrix of dimension q x p. Outcome vs
#'     predictor correlation matrix.}
#' }
#'
#' @seealso \code{\link{malmo_gibbs}}, \code{\link{malmo_summarize}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' mydata <- gen_missing_value_mv(
#'   n = 50, p = 100, q = 3, snr = 1,
#'   proportion = 0.20, type = "mnar0.33"
#' )
#' malmo_plot_correlation(
#'   Y = mydata$Y_true,
#'   X = mydata$X
#' )
#' }
malmo_plot_correlation <- function(Y, X) {
  ## Setup ----
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)
  q <- ncol(Y)
  p <- ncol(X)
  ## Color palette ----
  col_div <- colorRampPalette(c("#1E88E5", "white", "#E53935"))(100)
  ## Layout: 2 heatmaps + 1 legend row ----
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE),
         heights = c(9, 1))
  par(oma = c(0, 0, 3, 0))
  ## Plot 1: Outcome vs Outcome correlation ----
  cor_Y <- cor(Y)
  par(mar = c(5, 6, 3, 2))
  image(1:q, 1:q,
        t(cor_Y[q:1, ]),
        col  = col_div,
        xaxt = "n",
        yaxt = "n",
        xlab = "Outcome",
        ylab = "Outcome",
        main = "Outcome vs Outcome",
        zlim = c(-1, 1))
  axis(1, at = 1:q,
       labels = paste0("Y", 1:q),
       las = 2, cex.axis = 0.9)
  axis(2, at = 1:q,
       labels = paste0("Y", q:1),
       las = 2, cex.axis = 0.9)
  for (i in 1:q) {
    for (j in 1:q) {
      text(i, q - j + 1,
           sprintf("%.2f", cor_Y[j, i]),
           cex  = 0.85,
           col  = ifelse(abs(cor_Y[j, i]) > 0.5,
                         "white", "black"),
           font = 2)
    }
  }
  ## Plot 2: Outcome vs Predictor correlation ----
  cor_YX    <- cor(Y, X)
  cor_range <- max(abs(cor_YX))
  if (cor_range == 0) cor_range <- 1
  par(mar = c(5, 6, 3, 2))
  image(1:p, 1:q,
        t(cor_YX),
        col  = col_div,
        xaxt = "n",
        yaxt = "n",
        xlab = "Predictor Index",
        ylab = "Outcome",
        main = "Outcome vs Predictor",
        zlim = c(-cor_range, cor_range))
  axis(2, at = 1:q,
       labels = paste0("Y", 1:q),
       las = 2, cex.axis = 0.9)
  x_ticks <- seq(1, p, by = max(1, round(p / 20)))
  axis(1, at = x_ticks,
       labels = x_ticks,
       las = 2, cex.axis = 0.7)
  ## Legend panel ----
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = c("Positive", "Zero/Neutral", "Negative"),
         fill   = c("#E53935", "white", "#1E88E5"),
         border = "gray50",
         horiz  = TRUE,
         bty    = "n",
         cex    = 0.9)
  ## Overall title ----
  mtext("MALMO Correlation Analysis",
        outer = TRUE, cex = 1.3, font = 2)
  ## Summary ----
  cat("=== Correlation Summary ===\n")
  cat(sprintf("Outcomes    : %d\n", q))
  cat(sprintf("Predictors  : %d\n", p))
  cat("\nOutcome correlation matrix:\n")
  print(round(cor_Y, 3))
  cat(sprintf("\nMax |cor(Y, X)| : %.3f\n", cor_range))
  cat(sprintf("Mean |cor(Y, X)|: %.3f\n", mean(abs(cor_YX))))
  ## Return ----
  return(invisible(list(
    cor_Y  = cor_Y,
    cor_YX = cor_YX
  )))
}
