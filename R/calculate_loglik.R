#' estimate per-column thresholds implied by the extended rank likelihood
#'
#' The extended rank likelihood never estimates cutpoints, but any threshold in
#' the gap between adjacent levels of x is consistent with the observed ranks.
#' This takes the midpoint of each gap.
#'
#' @keywords internal
#'
#' @param x N x P latent Gaussian matrix.
#' @param dat data matrix.
#'
#' @return List with one numeric vector of thresholds per column, of length (number of observed levels - 1).
estimate_thresholds <- function(x, dat) {
    lapply(seq_len(ncol(dat)), function(j) {
        obs <- !is.na(dat[, j])
        levels <- sort(unique(dat[obs, j]))
        if (length(levels) < 2) return(numeric(0))
        vapply(seq_len(length(levels) - 1), function(v) {
            upper_of_lower <- max(x[obs & dat[, j] == levels[v], j])
            lower_of_upper <- min(x[obs & dat[, j] == levels[v + 1], j])
            (upper_of_lower + lower_of_upper) / 2
        }, numeric(1))
    })
}


#' log probability of each observed cell under an ordered probit with unit variance
#'
#' @keywords internal
#'
#' @param pred_mat N x P matrix of latent means.
#' @param dat data matrix.
#' @param thresholds list of per-column thresholds; levels are the sorted observed values of each column.
#'
#' @return N x P matrix of log probabilities, NA for missing cells.
cell_log_probs <- function(pred_mat, dat, thresholds) {
    out <- matrix(NA_real_, nrow(dat), ncol(dat))
    for (j in seq_len(ncol(dat))) {
        obs <- which(!is.na(dat[, j]))
        levels <- sort(unique(dat[obs, j]))
        cuts <- c(-Inf, thresholds[[j]][is.finite(thresholds[[j]])], Inf)
        v <- match(dat[obs, j], levels)
        mu <- pred_mat[obs, j]
        upper <- pnorm(cuts[v + 1] - mu)
        lower <- pnorm(cuts[v] - mu)
        out[obs, j] <- log(upper - lower)
    }
    out
}


#' calculate log-likelihood of the observed data at a posterior draw
#'
#' @keywords internal
#'
#' @param k_lambda stored lambda matrix from mcmc iteration.
#' @param k_omega stored omega matrix from mcmc iteration.
#' @param alpha stored alpha vector from mcmc iteration.
#' @param dat data matrix.
#' @param thresholds list of per-column thresholds, e.g. from \code{estimate_thresholds}.
#'
#' @return log-likelihood
calculate_log_lik <- function(k_lambda, k_omega, alpha, dat, thresholds) {
    # get dimensions
    c(n, p) %<-% get_dims_2d(dat)

    # calculate predictive (mean) matrix
    pred_mat <- t(k_lambda %*% t(k_omega)) -
        matrix(nrow = n, ncol = p, rep(alpha, n), byrow = TRUE)

    sum(cell_log_probs(pred_mat, dat, thresholds), na.rm = TRUE)
}
