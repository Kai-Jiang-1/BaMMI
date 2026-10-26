## Internal helpers for simulating missing outcome values.
## Shared by gen_missing_value() and gen_missing_value_mv().

## Resolve the MNAR ratio from either `mnar_ratio` or the legacy `type`.
## Returns a list with the numeric ratio and a label for `data_type`.
.resolve_mnar_ratio <- function(type = NULL, mnar_ratio = NULL) {

  legacy <- c(mar = 0, mnar0.33 = 1/3, mnar0.67 = 2/3, mnar = 1)

  if (!is.null(mnar_ratio)) {
    if (!is.numeric(mnar_ratio) || length(mnar_ratio) != 1 ||
        is.na(mnar_ratio) || mnar_ratio < 0 || mnar_ratio > 1) {
      stop("mnar_ratio must be a single number between 0 and 1.")
    }
    if (!is.null(type)) {
      warning("Both `type` and `mnar_ratio` supplied; using `mnar_ratio`.")
    }
    label <- if (mnar_ratio == 0) {
      "mar"
    } else if (mnar_ratio == 1) {
      "mnar"
    } else {
      sprintf("mnar%.2f", mnar_ratio)
    }
    return(list(mnar_ratio = mnar_ratio, label = label))
  }

  if (is.null(type)) {
    stop("Supply `mnar_ratio` (a number between 0 and 1) or `type`.")
  }
  if (length(type) != 1 || !type %in% names(legacy)) {
    stop('type must be one of "mnar", "mar", "mnar0.33", "mnar0.67".')
  }

  list(mnar_ratio = unname(legacy[type]), label = type)
}

## Choose which entries of Y_true (vector or matrix) become missing.
## The lowest round(n_miss * mnar_ratio) values are MNAR (below the LOD);
## the rest of the missing entries are drawn at random (MAR) from the
## remaining values.
.gen_missing_idx <- function(Y_true, proportion, mnar_ratio) {

  n_total <- length(Y_true)
  n_miss  <- round(n_total * proportion)
  n_mnar  <- round(n_miss * mnar_ratio)
  n_mar   <- n_miss - n_mnar

  mnar_idx  <- order(Y_true)[seq_len(n_mnar)]
  remaining <- setdiff(seq_len(n_total), mnar_idx)
  mar_idx   <- remaining[sample.int(length(remaining), n_mar)]

  list(
    missing_idx = c(mnar_idx, mar_idx),
    mnar_idx    = mnar_idx,
    mar_idx     = mar_idx
  )
}
