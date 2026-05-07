#' MCMC Diagnostics for MALMO
#'
#' Produces diagnostic plots and statistics for a single chain from
#' \code{\link{malmo_gibbs}}, including trace plots, density plots,
#' autocorrelation plots, effective sample size (ESS) and Geweke
#' convergence test for sigma2 (all q outcomes) and theta.
#'
#' @importFrom stats acf density
#' @importFrom coda mcmc effectiveSize geweke.diag
#' @importFrom graphics par layout plot lines abline mtext legend
#'   polygon title
#' @importFrom grDevices adjustcolor pdf dev.off
#'
#' @param result List. Output from \code{\link{malmo_gibbs}}.
#' @param nburn Integer. Number of burn-in iterations.
#' @param niter Integer. Number of posterior sampling iterations.
#' @param save_plot Character or NULL. If a file path is provided
#'   (e.g. \code{"diagnostics.pdf"}), saves the plot to that file.
#'   Default is \code{NULL} (display in viewer).
#' @param width Numeric. Width of the saved plot in inches.
#'   Default is 18.
#' @param height Numeric or NULL. Height of the saved plot in inches.
#'   If \code{NULL}, computed automatically. Default is \code{NULL}.
#'
#' @return A list containing:
#' \describe{
#'   \item{ESS_sigma2}{Numeric vector of ESS for each sigma2.}
#'   \item{ESS_theta}{Numeric. ESS for theta.}
#'   \item{geweke_sigma2}{Numeric vector of Geweke z-scores for
#'     each sigma2.}
#'   \item{geweke_theta}{Numeric. Geweke z-score for theta.}
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
#' result <- malmo_gibbs(
#'   Y = mydata$Y_miss, X = mydata$X,
#'   nburn = 1000, niter = 2000
#' )
#' ## Display in viewer
#' diag_result <- malmo_diagnostic(result, nburn = 1000, niter = 2000)
#'
#' ## Save to PDF
#' diag_result <- malmo_diagnostic(result, nburn = 1000, niter = 2000,
#'                                 save_plot = "diagnostics.pdf")
#' }
malmo_diagnostic <- function(result, nburn, niter,
                             save_plot = NULL,
                             width     = 18,
                             height    = NULL) {

  ## Extract posterior samples ----
  post_idx <- (nburn + 1):(nburn + niter)
  q        <- ncol(result$Sigma_result[[1]])

  sigma2_chains <- lapply(1:q, function(h) {
    sapply(result$Sigma_result[post_idx], function(S) S[h, h])
  })
  theta_chain <- result$theta_result[post_idx]

  ## Colors ----
  sigma_cols <- c("#E53935", "#1E88E5", "#43A047", "#FB8C00",
                  "#8E24AA", "#00ACC1", "#F4511E", "#3949AB")
  sigma_cols <- rep(sigma_cols, length.out = q)
  theta_col  <- "#7B1FA2"

  ## ESS and Geweke ----
  ESS_sigma2    <- numeric(q)
  geweke_sigma2 <- numeric(q)

  for (h in 1:q) {
    mcmc_chain       <- coda::mcmc(sigma2_chains[[h]])
    ESS_sigma2[h]    <- coda::effectiveSize(mcmc_chain)
    geweke_sigma2[h] <- coda::geweke.diag(mcmc_chain)$z
  }
  ESS_theta    <- coda::effectiveSize(coda::mcmc(theta_chain))
  geweke_theta <- coda::geweke.diag(coda::mcmc(theta_chain))$z

  ## Large q warning ----
  large_q <- q > 6
  do_plot <- !large_q || !is.null(save_plot)

  if (large_q && is.null(save_plot)) {
    message(
      "\n[!] NOTE: You have q = ", q, " outcomes (> 6).\n",
      "    Displaying all diagnostic plots in the viewer may be\n",
      "    unreadable. Consider saving to a file instead.\n",
      "    Re-run with, for example:\n\n",
      "      malmo_diagnostic(result, nburn, niter,\n",
      "                       save_plot = 'diagnostics.pdf',\n",
      "                       width = 18, height = ", (q + 1) * 4, ")\n\n",
      "    Skipping plot - showing diagnostic statistics only."
    )
    do_plot <- FALSE
  }

  ## Plotting ----
  if (do_plot) {

    ## Open file device if save_plot provided ----
    if (!is.null(save_plot)) {
      auto_height <- if (is.null(height)) (q + 1) * 4 else height
      ext <- tolower(tools::file_ext(save_plot))
      if (ext == "pdf") {
        pdf(save_plot, width = width, height = auto_height)
      } else if (ext %in% c("png", "jpg", "jpeg")) {
        grDevices::png(save_plot,
                       width  = width * 100,
                       height = auto_height * 100,
                       res    = 100)
      } else {
        pdf(save_plot, width = width, height = auto_height)
      }
      on.exit(dev.off(), add = TRUE)
      message("Saving diagnostic plot to: ", save_plot)
    }

    ## Save and restore graphics state ----
    old_par <- par(no.readonly = TRUE)
    on.exit({
      layout(matrix(1))
      par(old_par)
    }, add = TRUE)

    ## Build parameter list ----
    all_params <- c(
      lapply(1:q, function(h) list(type = "sigma2", h = h)),
      list(list(type = "theta"))
    )
    n_params     <- length(all_params)
    max_per_page <- if (q <= 6) n_params else 4
    n_pages      <- ceiling(n_params / max_per_page)

    ## Loop over pages ----
    for (pg in seq_len(n_pages)) {

      idx_start <- (pg - 1) * max_per_page + 1
      idx_end   <- min(pg * max_per_page, n_params)
      params_pg <- all_params[idx_start:idx_end]
      n_rows_pg <- length(params_pg)

      layout(matrix(seq_len(n_rows_pg * 3),
                    nrow  = n_rows_pg,
                    ncol  = 3,
                    byrow = TRUE))

      row_mar <- max(2, min(4, round(12 / n_rows_pg)))
      par(mar = c(row_mar, 4, row_mar, 1),
          oma = c(1, 0, 5, 0))

      for (pm in params_pg) {

        if (pm$type == "sigma2") {
          h     <- pm$h
          chain <- sigma2_chains[[h]]
          col_h <- sigma_cols[h]

          ## Trace ----
          plot(chain, type = "l", col = col_h, lwd = 1,
               main = paste0("Trace: sigma2[", h, "]"),
               xlab = "Iteration", ylab = "")
          abline(h = mean(chain), col = "red", lty = 2)

          ## Density ----
          d <- density(chain)
          plot(d, col = col_h, lwd = 2,
               main = paste0("Density: sigma2[", h, "]"),
               xlab = "")
          polygon(d, col = adjustcolor(col_h, 0.3), border = col_h)
          abline(v = mean(chain), col = "red", lty = 2)

          ## ACF ----
          acf(chain, main = "", col = col_h, lwd = 2)
          title(main = paste0("ACF: sigma2[", h, "]"), line = 1)

        } else {

          ## Trace theta ----
          plot(theta_chain, type = "l", col = theta_col, lwd = 1,
               main = "Trace: theta",
               xlab = "Iteration", ylab = "")
          abline(h = mean(theta_chain), col = "red", lty = 2)

          ## Density theta ----
          d <- density(theta_chain)
          plot(d, col = theta_col, lwd = 2,
               main = "Density: theta", xlab = "")
          polygon(d, col = adjustcolor(theta_col, 0.3),
                  border = theta_col)
          abline(v = mean(theta_chain), col = "red", lty = 2)

          ## ACF theta ----
          acf(theta_chain, main = "", col = theta_col, lwd = 2)
          title(main = "ACF: theta", line = 1)
        }
      }

      mtext(
        paste0("MALMO MCMC Diagnostics (Page ", pg,
               " of ", n_pages, ")"),
        outer = TRUE, line = 3, cex = 1.3, font = 2
      )
    }
  }

  ## Print summary ----
  cat("=== MALMO MCMC Diagnostics ===\n")
  cat(sprintf("Iterations used : %d\n\n", niter))

  for (h in 1:q) {
    cat(sprintf("sigma2[%d] - ESS: %6.1f | Geweke z: %6.3f %s\n",
                h, ESS_sigma2[h], geweke_sigma2[h],
                ifelse(abs(geweke_sigma2[h]) < 2,
                       "[OK]", "[WARNING]")))
  }
  cat(sprintf("theta     - ESS: %6.1f | Geweke z: %6.3f %s\n",
              ESS_theta, geweke_theta,
              ifelse(abs(geweke_theta) < 2, "[OK]", "[WARNING]")))

  cat("\nConvergence guide:\n")
  cat("  ESS > 100    : acceptable\n")
  cat("  ESS > 1000   : good\n")
  cat("  |Geweke| < 2 : converged [OK]\n")
  cat("  |Geweke| > 2 : not converged [WARNING]\n")

  ## Return ----
  return(invisible(list(
    ESS_sigma2    = ESS_sigma2,
    ESS_theta     = ESS_theta,
    geweke_sigma2 = geweke_sigma2,
    geweke_theta  = geweke_theta
  )))
}
