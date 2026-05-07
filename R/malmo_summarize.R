#' Summarize Posterior Samples from MALMO
#'
#' Computes posterior mean and 95 percent credible intervals from
#' the Gibbs sampler output of \code{\link{malmo_gibbs}}. Posterior
#' mean estimates are set to zero when the 95 percent credible interval
#' crosses zero.
#'
#' @importFrom stats quantile
#'
#' @param posterior_list List. A list of posterior sample matrices
#'   from \code{\link{malmo_gibbs}} (e.g. \code{result$C_result}).
#' @param nburn Integer. Number of burn-in iterations to discard.
#' @param niter Integer. Number of posterior sampling iterations.
#' @param thin Integer. Thinning interval for posterior samples.
#'   Default is 1 (no thinning).
#'
#' @return A list containing:
#' \describe{
#'   \item{mean}{Matrix of posterior mean estimates. Entries are set
#'     to zero when the 95 percent credible interval crosses zero.}
#'   \item{lower}{Matrix of lower 95 percent credible interval bounds.}
#'   \item{upper}{Matrix of upper 95 percent credible interval bounds.}
#' }
#'
#' @seealso \code{\link{malmo_gibbs}}, \code{\link{gen_missing_value_mv}}
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
#'
#' ## Summarize C matrix
#' C_summary <- malmo_summarize(
#'   result$C_result,
#'   nburn = 1000,
#'   niter = 2000
#' )
#' C_est <- C_summary$mean
#'
#' ## Summarize imputed Y
#' Y_summary <- malmo_summarize(
#'   result$Y_mvi_result,
#'   nburn = 1000,
#'   niter = 2000
#' )
#' Y_mvi <- Y_summary$mean
#' }
malmo_summarize <- function(posterior_list, nburn, niter, thin = 1) {

  ## Thinning indices ----
  total_iter      <- nburn + niter
  posterior_index <- seq(nburn, total_iter, by = thin)

  ## Convert list to array ----
  posterior_array <- simplify2array(posterior_list[posterior_index])

  ## Compute summaries ----
  mean_est <- apply(posterior_array, c(1, 2), mean)
  lower_ci <- apply(posterior_array, c(1, 2), quantile, probs = 0.025)
  upper_ci <- apply(posterior_array, c(1, 2), quantile, probs = 0.975)

  ## Set to zero if 95% CI crosses zero ----
  cross_zero           <- which(lower_ci * upper_ci < 0, arr.ind = TRUE)
  mean_est[cross_zero] <- 0

  ## Return ----
  return(list(
    mean  = mean_est,
    lower = lower_ci,
    upper = upper_ci
  ))
}
