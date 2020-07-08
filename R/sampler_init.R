# generate initial values for single chain
initialize_sampler <- function(
    dat,
    mode,
    n_levels,
    margin_vals,
    K_start = 100
    ) {
    # get dimensions
    c(n, p) %<-% get_dims_2d(dat)

    # set number of starting K
    K_start <- min(K_start, n - 1, p - 1)

    # sample factor loadings (p x K) from N(0, 1)
    lambda <- matrix(ncol = K_start, nrow = p, rnorm(K_start * p, 0, 1))
    # induce sparsity
    zeros  <- lambda
    zeros[abs(zeros) < .1]  <- 0
    zeros[abs(zeros) >= .1]  <- 1
    lambda  <- lambda * zeros
    # normalize
    lambda_rsums  <- sqrt(1 + apply(lambda^2, 1, sum))
    lambda  <- lambda / lambda_rsums

    # initialize alpha to 0 (p x 1)
    alpha  <- as.matrix(rep(0, p))

    # sample factor scores (n x K) from N(0, 1)
    omega <- matrix(nrow = n, ncol = K_start, rnorm(K_start * n, 0, 1))

    # sample x from truncated normal with fixed cutpoints based on n_levels
    alpha_mat <- matrix(nrow = n, ncol = p, alpha, byrow = TRUE)
    pred_mat <- t(lambda %*% t(omega)) - alpha_mat

    # sample for fixed margins
    if (mode == "fixed") {
        # calculate uniform cutpoints
        cutpoints <- seq(0, 1, 1 / n_levels)
        tnorm_cutpoints <- qnorm(cutpoints)

        # draws from truncated normal for observed data
        draws <- list()
		for (c_k in 1:n_levels) {
			draws[[c_k]] <- matrix(
                nrow = n,
                ncol = p,
                truncnorm::rtruncnorm(
                    n * p,
                    a = tnorm_cutpoints[c_k],
                    b = tnorm_cutpoints[c_k + 1],
                    mean = as.vector(pred_mat)
                )
            )
            draws[[c_k]] <- ((dat == margin_vals[c_k]) & is.na(dat) == FALSE) * draws[[c_k]]
        }
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
        }

        # collapse list of draws into matrix
        x <- Reduce("+", draws)

    } else {
        # sample draws for mixed margins
        # uniform cutpoints by margin
        cutpoints <- purrr::map(n_levels, function(x) seq(0, 1, 1 / x))
        tnorm_cutpoints <- lapply(cutpoints, qnorm)

        # convert data and means to tibbles for mapping
        dat_tbl <- tibble::as_tibble(dat, .name_repair = "minimal")
        pred_mat_tbl <- tibble::as_tibble(pred_mat, .name_repair = "minimal")

        # draws from truncated normals for observed data
        draws <- purrr::pmap(
            list(n_levels, tnorm_cutpoints, pred_mat_tbl, dat_tbl, margin_vals),
            function(n_levels, cutpoints, pred_mat, dat, margin_vals) {
                draws <- list()
                for (c_k in 1:n_levels) {
                    draws[[c_k]] <- matrix(
                        nrow = n,
                        ncol = 1,
                        truncnorm::rtruncnorm(
                            n,
                            a = cutpoints[c_k],
                            b = cutpoints[c_k + 1],
                            mean = pred_mat
                        )
                    )
                    draws[[c_k]] <- ((dat == margin_vals[c_k]) & is.na(dat) == FALSE) *
                        draws[[c_k]]
                }
                Reduce("+", draws)
            }
        )

        # collapse col-wise draws into matrix
        x <- do.call(cbind, draws)

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
            x  <- x + draws_na
        }
    }


    # gamma row vector (1 x K) set to 1
    gamma_k  <- rep(1, K_start)

    # d
    d <- 100

    # dc
    dc <- rep(0, K_start)

    # Indian Buffet alpha and beta hyperpriors (a, b)
    ibp_a <- ibp_b <- 1

    # store values
    output_list <- list()
    output_list[["x"]] <- x
    output_list[["lambda"]] <- lambda
    output_list[["zeros"]] <- zeros
    output_list[["omega"]] <- omega
    output_list[["alpha"]] <- alpha
    output_list[["gamma_k"]] <- gamma_k
    output_list[["d"]] <- d
    output_list[["dc"]] <- dc
    output_list[["ibp_a"]] <- ibp_a
    output_list[["ibp_b"]] <- ibp_b

    return(output_list)
}