# sample X values according to Hoff (2007) and Murray et al. (2013) extended rank likelihood approach
# 
# @param dat Data matrix.
# @param mode Margin type. Must be one of "fixed", "mixed" or "multi."
# @param d_mask Boolean mask of extended rank likelihood.
# @param n_levels Numeric vector of number of unique values per margin.
# @param margin_vals List of unique values per margin.
# @param x N x P Gaussian copula matrix.
# @param lambda P x K matrix of factor loadings.
# @param alpha A numeric vector of item-level intercepts.
# @param omega N x K matrix of factor scores.
# 
# @return List of x matrix and cutpoints.
# 
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

    # convert inputs to tibbles for mapping
    c(x_tbl, dat_tbl, pred_mat_tbl) %<-% map(
        list(x, dat, pred_mat),
        tibble::as_tibble,
        .name_repair = "minimal"
    )

    # extract D masks for lower and upper bounds
    c(d_mask_l, d_mask_u) %<-% d_mask

    # Vary sampling approach by margin type
    # below implementation should be slightly faster for fixed margins
    # where entire x matrix can be draw in single step for each level
    if (mode == "fixed") {
        # get lower cutpoints
        margin_bounds_l <- map(
            d_mask_l,
            ~ map2_dbl(
                x_tbl, .x,
                function(x, y) max(x[y], na.rm = TRUE)
                )
            )
        # upper cutpoints
        margin_bounds_u <- map(
            d_mask_u,
            ~ map2_dbl(
                x_tbl, .x,
                function(x, y) min(x[y], na.rm = TRUE)
            )
        )
        # combine into single list
        boundary_lists <- map2(
            margin_bounds_l, margin_bounds_u,
            ~ list("a" = .x, "b" = .y)
        )

        # Draw samples for observed margins
        samples <- map2(
            boundary_lists, margin_vals,
            function(bounds, vals) {
                c(a, b) %<-% map(bounds, rep, each = n)
                draws <- matrix(
                    nrow = n,
                    ncol = p,
                    truncnorm::rtruncnorm(
                        n * p,
                        a = a,
                        b = b,
                        mean = as.vector(pred_mat)
                    )
                )

                # filter draws by mask
                draws <- ((dat == vals) & (is.na(dat) == FALSE)) * draws

                # filter cutpoints by mask
                c(mat_a, mat_b) %<-% map(
                    list(a, b),
                    function(x) {
                        x <- ((dat == vals) & (is.na(dat) == FALSE)) * x
                        x[which(is.nan(x))] <- 0
                        return(x)
                    }
                )

                # store output
                output_list <- list()
                output_list[["draws"]] <- draws
                output_list[["mat_a"]] <- mat_a
                output_list[["mat_b"]] <- mat_b
                return(output_list)
            }
        )
        draws <- map(samples, ~ .x[["draws"]])
        mat_a <- map(samples, ~ .x[["mat_a"]])
        mat_b <- map(samples, ~ .x[["mat_b"]])
        mat_a <- Reduce('+', mat_a)
        mat_b <- Reduce('+', mat_b)

        # draws for NA (if needed)
        if (any(is.na(dat))) {
            draws[[n_levels + 1]] <- matrix(
                nrow = n,
                ncol = p,
                rnorm(
                    n * p,
                    mean = as.vector(pred_mat),
                    sd = 1
                )
            )
            draws[[n_levels + 1]] <- (is.na(dat) == TRUE) * draws[[n_levels + 1]]
            mat_a[as.vector((is.na(dat) == TRUE))] <- -Inf
            mat_b[as.vector((is.na(dat) == TRUE))] <- Inf
        }

        # combine draws
        x_out <- Reduce('+', draws)

        # combine boundaries to get cutpoints list
        cutpoints_out <- list("a" = mat_a, "b" = mat_b)

    } else { # multimodal margins

        # get lower cutpoints
        margin_bounds_l <- map2(
            x_tbl,
            d_mask_l,
            function(x, y) {
                l_bounds <- map_dbl(
                    y,
                    function(y) max(x[y], na.rm = T)
                )
                unname(l_bounds)
            }
        )
        # get upper cutpoints
        margin_bounds_u <- map2(
            x_tbl,
            d_mask_u,
            function(x, y) {
                u_bounds <- map_dbl(
                    y,
                    function(y) min(x[y], na.rm = T)
                )
                unname(u_bounds)
            }
        )
        # combine into single list
        boundary_lists <- map2(
            margin_bounds_l, margin_bounds_u,
            ~ list('a' = .x, 'b' = .y)
        )

        # draw samples for observed data
        samples <- pmap(
            list(n_levels, boundary_lists, pred_mat_tbl, dat_tbl, margin_vals),
            function(n_levels, cutpoints, pred_mat, dat, margin_vals) {
                c(a, b) %<-% cutpoints
                draws <- mat_a <- mat_b <- list()
                for (c_k in 1:n_levels) {
                    draws[[c_k]] <- matrix(
                        nrow = n,
                        ncol = 1,
                        truncnorm::rtruncnorm(
                            n,
                            a = a[c_k],
                            b = b[c_k],
                            mean = pred_mat
                        )
                    )
                    draws[[c_k]] <- ((dat == margin_vals[c_k]) & is.na(dat) == FALSE) *
                        draws[[c_k]]
                    mat_a[[c_k]] <- matrix(
                        nrow = n,
                        ncol = 1,
                        ((dat == margin_vals[c_k]) & is.na(dat) == FALSE) *
                            a[c_k]
                    )
                    mat_a[[c_k]][which(is.nan(mat_a[[c_k]]))] <- 0
                    mat_b[[c_k]] <- matrix(
                        nrow = n,
                        ncol = 1,
                        ((dat == margin_vals[c_k]) & is.na(dat) == FALSE) *
                            b[c_k]
                    )
                    mat_b[[c_k]][which(is.nan(mat_b[[c_k]]))] <- 0
                }
                c(draws, mat_a, mat_b) %<-% map(
                    list(draws, mat_a, mat_b),
                    ~ Reduce("+", .x)
                )

                # store output
                output_list <- list()
                output_list[["draws"]] <- draws
                output_list[["mat_a"]] <- mat_a
                output_list[["mat_b"]] <- mat_b
                return(output_list)
            }
        )
        draws <- map(samples, ~ .x[["draws"]])
        mat_a <- map(samples, ~ .x[["mat_a"]])
        mat_b <- map(samples, ~ .x[["mat_b"]])

        # collapse col-wise draws into matrix
        c(x_out, mat_a, mat_b) %<-% map(
            list(draws, mat_a, mat_b),
            ~ do.call(cbind, .x)
        )
        # draws from full normal if needed
        if (any(is.na(dat))) {
            draws_na <- matrix(
                nrow = n,
                ncol = p,
                rnorm(
                    n * p,
                    mean = as.vector(pred_mat),
                    sd = 1
                )
            )
            draws_na <- (is.na(dat) == TRUE) * draws_na
            x_out <- x_out + draws_na
            mat_a[as.vector((is.na(dat) == TRUE))] <- -Inf
            mat_b[as.vector((is.na(dat) == TRUE))] <- Inf
        }

        # get cutpoints
        cutpoints_out <- list("a" = mat_a, "b" = mat_b)
    }

    # store outputs
    output_list <- list()
    output_list[["x"]] <- x_out
    output_list[["cutpoints"]] <- cutpoints_out

    return(output_list)
}
