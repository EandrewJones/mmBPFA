#' calculate log-likelihood by integrating over the posterior distribution
#'
#' @keywords internal
#' 
#' @param k_lambda stored lambda matrix from mcmc iteration.
#' @param k_omega stored omega matrix from mcmc iteration.
#' @param alpha stored alpha vector from mcmc iteration.
#' @param dat data matrix.
#' 
#' @return log-likelihood
calculate_log_lik <- function(k_lambda, k_omega, alpha, dat) {
    # get dimensions
    c(n, p) %<-% get_dims_2d(dat)

    # get vector-wise indices for non-missing data
    where_data <- (is.na(as.vector(dat)) == FALSE)

    # calculate predictive (mean) matrix
    pred_mat <- t(k_lambda %*% t(k_omega)) -
        matrix(nrow = n, ncol = p, rep(alpha, n), byrow = TRUE)
    log_lik <- pnorm(
        0,
        mean = as.vector(pred_mat)[where_data],
        sd = 1,
        lower.tail = as.vector(dat)[where_data],
        log.p = TRUE
    )
    log_lik <- sum(log_lik)
    return(log_lik)
}