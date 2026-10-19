#' sample X values according to Hoff (2007) and Murray et al. (2013) extended rank likelihood approach
#'
#' @keywords internal
#'
#' @param dat Data matrix.
#' @param mode Margin type. Must be one of "fixed", "mixed" or "multi."
#' @param d_mask Boolean mask of extended rank likelihood.
#' @param n_levels Numeric vector of number of unique values per margin.
#' @param margin_vals List of unique values per margin.
#' @param x N x P Gaussian copula matrix.
#' @param lambda P x K matrix of factor loadings.
#' @param alpha A numeric vector of item-level intercepts.
#' @param omega N x K matrix of factor scores.
#' 
#' @return List of x matrix and cutpoints (the lower and upper truncation bounds used for each cell).
#' 
sample_x <- function(
    dat,
    mode,
    d_mask,
    n_levels,
    margin_vals,
    x,
    lambda,
    alpha,
    omega
    ) {
    # get dimensions
    c(n, p) %<-% get_dims_2d(dat)

    # check if equal to Murray et al. (2013) specification
    alpha_mat <- matrix(nrow = n, ncol = p, alpha, byrow = TRUE)
    pred_mat <- t(lambda %*% t(omega)) - alpha_mat

    # extract D masks for lower and upper bounds
    c(d_mask_l, d_mask_u) %<-% d_mask

    # truncation bounds used for each cell (missing cells are unbounded)
    mat_a <- matrix(-Inf, nrow = n, ncol = p)
    mat_b <- matrix(Inf, nrow = n, ncol = p)
    observed <- !is.na(dat)

    # Within each column, levels are updated one at a time so every draw is
    # bounded by the current values of the neighbouring levels. Cells sharing a
    # level do not constrain each other, so each level is a valid Gibbs block.
    for (j in seq_len(p)) {
        vals <- if (mode == "fixed") margin_vals else margin_vals[[j]]
        for (v in seq_along(vals)) {
            cells <- which(observed[, j] & dat[, j] == vals[v])
            if (length(cells) == 0) next
            if (mode == "fixed") {
                below <- which(d_mask_l[[v]][[j]])
                above <- which(d_mask_u[[v]][[j]])
            } else {
                below <- which(d_mask_l[[j]][[v]])
                above <- which(d_mask_u[[j]][[v]])
            }
            a <- if (length(below)) max(x[below, j]) else -Inf
            b <- if (length(above)) min(x[above, j]) else Inf
            x[cells, j] <- truncnorm::rtruncnorm(
                length(cells),
                a = a,
                b = b,
                mean = pred_mat[cells, j]
            )
            mat_a[cells, j] <- a
            mat_b[cells, j] <- b
        }

        # draws from full normal for missing data
        missing <- which(!observed[, j])
        x[missing, j] <- rnorm(length(missing), mean = pred_mat[missing, j], sd = 1)
    }

    # store outputs
    output_list <- list()
    output_list[["x"]] <- x
    output_list[["cutpoints"]] <- list("a" = mat_a, "b" = mat_b)

    return(output_list)
}
