#' Venn Diagram or UpSet Plot for Selected Features in MALMO
#'
#' For q <= 3 outcomes, produces a Venn diagram showing overlap of
#' selected features across outcomes. For q > 3 outcomes, produces
#' an UpSet plot using base R.
#'
#' @importFrom graphics par plot rect text segments barplot axis
#'   mtext title points plot.new legend layout polygon
#' @importFrom grDevices adjustcolor dev.size
#'
#' @param C_summary List. Output from \code{\link{malmo_summarize}}
#'   applied to \code{result$C_result}.
#'
#' @return A list containing:
#' \describe{
#'   \item{selected_per_outcome}{List of selected feature indices
#'     for each outcome.}
#'   \item{overlap_matrix}{Binary matrix showing which features
#'     are selected for which outcomes.}
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
#' C_summary <- malmo_summarize(result$C_result,
#'                              nburn = 1000, niter = 2000)
#' malmo_plot_venn(C_summary)
#' }
malmo_plot_venn <- function(C_summary) {

  ## Setup ----
  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(matrix(1))
    par(old_par)
  }, add = TRUE)

  C_est <- C_summary$mean
  q     <- ncol(C_est)
  p     <- nrow(C_est)

  ## Single accent color ----
  accent_col <- "#546E7A"

  ## Selected features per outcome ----
  selected_per_outcome <- lapply(1:q, function(j)
    which(C_est[, j] != 0))

  ## Overlap matrix ----
  all_selected <- sort(unique(unlist(selected_per_outcome)))

  if (length(all_selected) == 0) {
    message("No features selected - nothing to plot.")
    return(invisible(NULL))
  }

  overlap_mat <- matrix(0,
                        nrow = length(all_selected),
                        ncol = q)
  rownames(overlap_mat) <- paste0("F", all_selected)
  colnames(overlap_mat) <- paste0("Y", 1:q)

  for (j in 1:q)
    overlap_mat[paste0("F", selected_per_outcome[[j]]), j] <- 1

  ## Venn diagram for q <= 3 ----
  if (q <= 3) {

    dev_width  <- dev.size("in")[1]
    dev_height <- dev.size("in")[2]

    if (dev_width > dev_height) {
      margin <- (1 - dev_height / dev_width) / 2
      par(fig = c(margin, 1 - margin, 0, 1), new = FALSE)
    } else {
      margin <- (1 - dev_width / dev_height) / 2
      par(fig = c(0, 1, margin, 1 - margin), new = FALSE)
    }

    par(mar = c(3, 3, 4, 3))
    plot(0, 0, type = "n",
         xlim = c(-2, 2), ylim = c(-2, 2),
         xlab = "", ylab = "",
         main = "Selected Features: Venn Diagram",
         asp  = 1, axes = FALSE)

    if (q == 2) {
      cx <- c(-0.5, 0.5)
      cy <- c(0, 0)
    } else {
      cx <- c(-0.6, 0.6, 0)
      cy <- c(-0.4, -0.4, 0.6)
    }
    r <- 0.9

    theta_seq <- seq(0, 2 * pi, length.out = 200)
    for (j in 1:q) {
      polygon(cx[j] + r * cos(theta_seq),
              cy[j] + r * sin(theta_seq),
              col    = adjustcolor(accent_col, alpha.f = 0.3),
              border = accent_col,
              lwd    = 2)
      text(cx[j] + r * 1.15 * cos(pi/2 + (j-1) * 2*pi/q),
           cy[j] + r * 1.15 * sin(pi/2 + (j-1) * 2*pi/q),
           paste0("Y", j, "\n(n=", length(selected_per_outcome[[j]]), ")"),
           cex  = 0.9,
           font = 2,
           col  = accent_col)
    }

    if (q == 2) {
      s1   <- selected_per_outcome[[1]]
      s2   <- selected_per_outcome[[2]]
      both <- intersect(s1, s2)
      text(-0.8, 0, length(setdiff(s1, s2)), cex = 1.2, font = 2)
      text( 0.8, 0, length(setdiff(s2, s1)), cex = 1.2, font = 2)
      text( 0.0, 0, length(both),            cex = 1.2, font = 2)

    } else if (q == 3) {
      s1 <- selected_per_outcome[[1]]
      s2 <- selected_per_outcome[[2]]
      s3 <- selected_per_outcome[[3]]

      only1 <- length(setdiff(s1, union(s2, s3)))
      only2 <- length(setdiff(s2, union(s1, s3)))
      only3 <- length(setdiff(s3, union(s1, s2)))
      s12   <- length(setdiff(intersect(s1, s2), s3))
      s13   <- length(setdiff(intersect(s1, s3), s2))
      s23   <- length(setdiff(intersect(s2, s3), s1))
      s123  <- length(intersect(intersect(s1, s2), s3))

      text(-0.9, -0.2, only1, cex = 1.1, font = 2)
      text( 0.9, -0.2, only2, cex = 1.1, font = 2)
      text( 0.0,  0.9, only3, cex = 1.1, font = 2)
      text( 0.0, -0.5, s12,   cex = 1.1, font = 2)
      text(-0.5,  0.3, s13,   cex = 1.1, font = 2)
      text( 0.5,  0.3, s23,   cex = 1.1, font = 2)
      text( 0.0,  0.1, s123,  cex = 1.1, font = 2)
    }

    ## UpSet plot for q > 3 ----
  } else {

    combo_list  <- list()
    combo_sizes <- c()
    combo_names <- c()

    for (k in 1:(2^q - 1)) {
      bits    <- as.integer(intToBits(k))[1:q]
      members <- which(bits == 1)

      in_all <- Reduce(intersect, selected_per_outcome[members])
      not_in <- setdiff(1:q, members)
      if (length(not_in) > 0) {
        exclude <- Reduce(union, selected_per_outcome[not_in])
        in_all  <- setdiff(in_all, exclude)
      }
      if (length(in_all) > 0) {
        combo_list[[length(combo_list) + 1]] <- members
        combo_sizes <- c(combo_sizes, length(in_all))
        combo_names <- c(combo_names,
                         paste0("Y", members, collapse = "&"))
      }
    }

    if (length(combo_sizes) == 0) {
      message("No overlap to plot.")
      return(invisible(NULL))
    }

    ord         <- order(combo_sizes, decreasing = TRUE)
    combo_list  <- combo_list[ord]
    combo_sizes <- combo_sizes[ord]
    combo_names <- combo_names[ord]
    n_combos    <- length(combo_sizes)

    layout(matrix(c(1, 2, 3), nrow = 3),
           heights = c(2, 1, 1))

    par(mar = c(0, 6, 4, 2))
    bp <- barplot(combo_sizes,
                  col       = accent_col,
                  border    = NA,
                  ylab      = "Intersection Size",
                  main      = "Selected Features: UpSet Plot",
                  ylim      = c(0, max(combo_sizes) * 1.2),
                  names.arg = rep("", n_combos),
                  xaxt      = "n")

    text(bp, combo_sizes + max(combo_sizes) * 0.05,
         combo_sizes, cex = 0.8, font = 2)

    par(mar = c(4, 6, 0, 2))
    plot(NA,
         xlim = c(0.5, n_combos + 0.5),
         ylim = c(0.5, q + 0.5),
         xaxt = "n", yaxt = "n",
         xlab = "", ylab = "")

    axis(2, at = 1:q,
         labels   = paste0("Y", q:1),
         las      = 2,
         cex.axis = 0.85)

    for (j in 1:q)
      rect(0.5, j - 0.5, n_combos + 0.5, j + 0.5,
           col    = ifelse(j %% 2 == 0, "#F5F5F5", "white"),
           border = NA)

    for (k in 1:n_combos)
      for (j in 1:q)
        points(k, q - j + 1,
               pch = 16, cex = 1.8, col = "gray85")

    for (k in 1:n_combos) {
      members <- combo_list[[k]]
      y_pos   <- q - members + 1

      if (length(members) > 1)
        segments(k, min(y_pos), k, max(y_pos),
                 col = accent_col, lwd = 3)

      points(rep(k, length(members)), y_pos,
             pch = 16, cex = 1.8, col = accent_col)
    }

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
           legend = c("Selected combination", "Not in combination"),
           pch    = 16,
           col    = c(accent_col, "gray85"),
           horiz  = TRUE,
           bty    = "n",
           cex    = 0.85)
  }

  ## Summary ----
  cat("=== Selection Summary ===\n")
  for (j in 1:q) {
    cat(sprintf("Y%-2d: %d selected features\n",
                j, length(selected_per_outcome[[j]])))
  }
  all_union <- length(unique(unlist(selected_per_outcome)))
  all_inter <- length(Reduce(intersect, selected_per_outcome))
  cat(sprintf("Union        : %d features\n", all_union))
  cat(sprintf("Intersection : %d features\n", all_inter))

  ## Return ----
  return(invisible(list(
    selected_per_outcome = selected_per_outcome,
    overlap_matrix       = overlap_mat
  )))
}
