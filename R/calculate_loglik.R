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
