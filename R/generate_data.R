#' Generate data for running simulations with mmBPFA
#' 
#' \code{generate_sim_data} produces simulated data that can be used for testing 
#' the mmBPFA sampler prior to production runs to ensure it will behave as expected.
#' 
#' @param mode Sets the type of margins for the simulated data. Must be one of "fixed", "multi" or "mixed." Default is "fixed."
#' @param n_levels Max number of unique values per margin. When mode is "fixed" n_levels is constant across margins. When mode is "multi" a margin may take on unique values up to n_levels. When mode is "mixed", n_levels determines max number of arbitrary distributions data may be drawn from. Default is 2.
#' @param n Number of observations. Default is 1000.
#' @param p Number of features. Default is 500.
#' @param K_true True number of latent dimensions. Default is 10.
#' @param sparsity Numeric value between 0 and 1 that determines amount of sparsity in the factor loadings. If sparsity argument is set, zeros are uniformly distributed across all dimensions. Default is NULL.
#' @param a A hyperparameter determining the distribution of sparsity across the factor loadings. Default is 10
#' @param b A hyperparameter determining the distribution of sparsity across the factor loadings. Default is 1.
#' @param d A hyperparameter determining the shape of the factor precisions (gamma_k). Default is 2.
#' @param e A hyperparameter determining the scale of the factor precisions (gamma_k). Default is 1.
#' @param prop_missing Numeric value between 0 and 1 determining the proportion of missingness in the data. Default is 0.
#' @param noise Numeric value between 0 and 1 determining how much noise should be induced into the data. Default is 0.05.
#' @param seed Random seed. Default is 123.
#' 
#' @return S3 object of class `mmBPFA.simulated.data` containing generated data and true underlying values. Attributes note the user-given arguments.
#' 
#' @export
#' 
generate_sim_data <- function(
    mode = "fixed",
    n_levels = 2,
    n = 1000,
    p = 500,
    K_true = 10,
    sparsity = NULL,
    a = 10,
    b = 1,
    d = 2,
    e = 1,
    prop_missing = 0,
    noise = 0.05,
    seed = 123
    ) {
    # Assertions
    avail_modes <- c("fixed", "multi", "mixed")
    if (!(mode %in% avail_modes)) stop("Mode parameter must be 'fixed', 'multi' or 'mixed'.")
    if (n_levels < 2) stop("n_levels parameter must be greater than 2.")
    if (mode == "multi" & n_levels == 2) {
        warning(
            "Mode = 'multi' with n_levels = 2 is the same as 'fixed' binary margins.
                            There cannot be less than two choices per column."
        )
    }
    if (mode == "mixed" & n_levels > 14) {
        warning("There are only 14 possible distributions to sample from for 'mixed'-margin type data. n_levels will be set to 14.")
        n_levels <- 14
    }
    if (!is.null(sparsity)) {
        if (sparsity < 0 | sparsity >= 1) stop("Sparsity parameter must be 0 <= x < 1")
        cat("Sparsity parameter included. Ignoring parameters a and b.")
    }
    if (noise < 0) stop("Noise parameter must be 0 or greater.")
    if(a <= 0 | b <= 0 | d <= 0 | e <= 0) stop("Parameters a, b, d, and e must be positive values.")

    # seed
    set.seed(seed)

    # get feature levels
    if (mode == "multi") feature_levels <- sample(2:n_levels, p, replace = T)

    # Should we use exact model or unif sampling for t_lambda?
    # generate factor precisions and lambda
    t_gamma_k <- rgamma(K_true, shape = d, scale = e)
    t_lambda <- matrix(
        nrow = p,
        ncol = K_true,
        rnorm(p * K_true, mean = 0, sd = 1 / t_gamma_k),
        byrow = T
    )

    # generate factor precisions
    if (!is.null(sparsity)) {
        t_pi_k <- sparsity
    } else {
        t_pi_k <- rbeta(K_true, a / K_true, b * (K_true - 1) / K_true)
    }

    # generate sparsity matrix
    t_zero <- matrix(
        nrow = p,
        ncol = K_true,
        rbinom(K_true * p, 1, t_pi_k),
        byrow = T
    )

    # induce sparsity
    t_lambda <- t_lambda * t_zero
    lambda_rsums <- sqrt(1 + apply(t_lambda^2, 1, sum))
    t_lambda <- t_lambda / lambda_rsums

    # generate lambda
    t_omega <- matrix(
        nrow = n,
        ncol = K_true,
        rnorm(K_true * n, 0, 1)
    )
    t_alpha <- matrix(
        nrow = p,
        ncol = 1,
        runif(p, -0.1, 0.1)
    )
    t_am <- matrix(
        nrow = n,
        ncol = p,
        rep(t_alpha, n),
        byrow = T
    )
    t_aug <- t(t_lambda %*% t(t_omega)) - t_am

    # copula transformation from normal margins to uniform margins
    uniform_margins <- pnorm(t_aug)

    # if binary use q_binom
    if (mode == "fixed") {
        if (n_levels == 2) {
            y <- qbinom(p = uniform_margins, size = 1, prob = 0.5)
            # center at 0
            y <- apply(y, 2, function(x) x - 0.5)
        } else {
            # transform back to multinomial margins based on cutpoints
            uniform_margins <- tibble::as_tibble(
                uniform_margins,
                .name_repair = "minimal"
            )

            # generate cuptoints based on number of feature levels
            cutpoints <- seq(0, 1, 1 / n_levels)[-1]
            y <- map(
                uniform_margins,
                vectorized_which.max,
                y = cutpoints) %>%
                dplyr::bind_cols() %>%
                as.matrix()
            colnames(y) <- NULL

            # center at 0
            y <- apply(y, 2, function(x) x - (n_levels - 1) / 2)
        }
    } else if (mode == "multi") {
        if (n_levels == 2) {
            y <- qbinom(p = uniform_margins, size = 1, prob = 0.5)
            # center at 0
            y <- apply(y, 2, function(x) x - 0.5)
        } else {
            # transform back to multinomial margins based on cutpoints
            uniform_margins <- tibble::as_tibble(
                uniform_margins,
                .name_repair = "minimal"
            )

            # generate a sequence of cutpoints based on number feature
            # levels per column
            cutpoints <- lapply(feature_levels, function(x) seq(0, 1, 1 / x)[-1])
            y <- purrr::map2(
                uniform_margins,
                cutpoints,
                vectorized_which.max) %>%
                purrr::map2(
                    .,
                    feature_levels,
                    ~ .x - ((.y - 1) / 2)
                ) %>%
                dplyr::bind_cols() %>%
                as.matrix()
            colnames(y) <- NULL
        }
    } else {
        # convert copula to tibble for mapping
        uniform_margins <- tibble::as_tibble(
                uniform_margins,
                .name_repair = "minimal"
            )
        dists <- c(
            "qnorm", "qlnorm", "qgamma", "qbeta", "qchisq",
            "qt", "qcauchy", "qexp", "qlogis", "qweibull",
            "qbinom", "qnbinom", "qpois", "qgeom"
        )
        # set hyperparameters
        beta_prior <- rbeta(1, 1, 1)
        gamma_prior <- rgamma(1, d, 1)
        df_prior <- rpois(1, gamma_prior + 10)
        hyperparams <- list(
            "qnorm" = list("mean" = 0, "sd" = 1),
            "qlnorm" = list("meanlog" = 0, "sdlog" = 1),
            "qgamma" = list("shape" = d, "scale" = 1),
            "qbeta" = list("shape1" = a, "shape2" = b),
            "qchisq" = list("df" = df_prior),
            "qt" = list("df" = df_prior),
            "qcauchy" = list("location" = 0, "scale" = 1),
            "qexp" = list("rate" = gamma_prior),
            "qlogis" = list("location" = 0, "scale" = 1),
            "qweibull" = list("shape" = 2, "scale" = 1),
            "qbinom" = list("size" = 1, "prob" = 0.5),
            "qnbinom" = list("size" = 1, "prob" = 0.2),
            "qpois" = list("lambda" = gamma_prior),
            "qgeom" = list("prob" = beta_prior)
        )

        # draw possible margin distributions
        dist_draw <- sample(
            dists,
            size = n_levels,
            replace = F,
            prob = rep(1 / length(dists), length(dists))
        )

        # draw margin distributions
        margin_dists <- sample(
            dist_draw,
            size = ncol(uniform_margins),
            replace = T,
            prob = rep(1 / n_levels, n_levels)
        )

        # function to draw samples
        margin_draws <- function(dist_type, copula) {
            # setup arguments
            f <- get(dist_type)
            args <- hyperparams[[dist_type]]
            margins <- list("p" = copula)

            # make draws
            margin_draws <- purrr::lift_dl(f)(c(margins, args))
            return(margin_draws)
        }

        # make draws
        y <- purrr::map2(
            margin_dists,
            uniform_margins,
            ~ margin_draws(.x, .y)) %>%
            dplyr::bind_cols() %>%
            as.matrix()
        colnames(y) <- NULL
    }

    # induce noise into observations
    if (noise > 0) {

        # sample random noise mask
        noise_idx <- matrix(
            nrow = n,
            ncol = p,
            sample(c(T, F), n * p, replace = T, prob = c(noise, 1 - noise))
        )

        # get random samplings from margins
        if (mode == "multi") {
            feature_level_centers <- sapply(feature_levels, function(x) {
                (x - 1) / 2
            })
            noise_samples <- sapply(feature_levels, function(x) {
                sample(
                    0:(x - 1),
                    size = n,
                    prob = seq(0, 1, 1 / x)[-1],
                    replace = T
                )
            })
            noise_samples <- t(t(noise_samples) - feature_level_centers)
        } else if (mode == "fixed") {
            feature_level_centers <- (n_levels - 1) / 2
            noise_samples <- matrix(
                nrow = n,
                ncol = p,
                sample(
                    0:(n_levels - 1),
                    size = n * p,
                    prob = seq(0, 1, 1 / n_levels)[-1],
                    replace = T
                )
            )
            noise_samples <- noise_samples - feature_level_centers
        } else {
            # subset margins with added noise
            n_draws_margin <- colSums(noise_idx)
            to_draw <- which(n_draws_margin > 0)
            n_draws <- n_draws_margin[to_draw]
            margin_dists_noise <- margin_dists[to_draw]

            # function to draw noise
            margin_noise_draws <- function(dist_type, n_draws) {
                # setup arguments
                f <- get(sub("^q", "r", dist_type))
                args <- hyperparams[[dist_type]]
                n <- list("n" = n_draws)

                # make draws
                margin_noise_draws <- purrr::lift_dl(f)(c(n, args))
                return(margin_noise_draws)
            }

            noise_samples <- purrr::map2(
                margin_dists_noise,
                n_draws,
                ~ margin_noise_draws(.x, .y)
            ) %>%
            do.call(c, .)
        }

        # add noise where mask == T
        if (mode == "mixed") {
            y[noise_idx]  <- noise_samples
        } else {
            y[noise_idx] <- noise_samples[noise_idx]
        }
    }

    # induce missingness
    if (prop_missing > 0) {
        missing_mask <- matrix(
            nrow = n,
            ncol = p,
            sample(
                c(1, NA),
                size = n * p,
                replace = T,
                prob = c(1 - prop_missing, prop_missing)
            )
        )
        y <- missing_mask * y
    }

    # Store outputs
    true_vals <- list()
    true_vals[["x"]] <- t_aug
    true_vals[["lambda"]] <- t_lambda
    true_vals[["zeros"]] <- t_zero
    true_vals[["omega"]] <- t_omega
    true_vals[["alpha"]] <- t_alpha
    true_vals[["pi_k"]] <- t_pi_k
    true_vals[["gamma_k"]] <- t_gamma_k

    output_list <- list()
    output_list[["data"]] <- y
    output_list[["true_values"]] <- true_vals

    # add attributes
    class(output_list) <- "mmBPFA.simulated.data"
    attr(output_list, "mode") <- mode
    attr(output_list, "n_levels") <- n_levels
    if (mode == "mixed") attr(output_list, "margin_dists") <- margin_dists
    attr(output_list, "n") <- n
    attr(output_list, "p") <- p
    attr(output_list, "K_true") <- K_true
    attr(output_list, "sparsity") <- sparsity
    attr(output_list, "a") <- a
    attr(output_list, "b") <- b
    attr(output_list, "d") <- d
    attr(output_list, "e") <- e
    attr(output_list, "prop_missing") <- prop_missing
    attr(output_list, "noise") <- noise
    attr(output_list, "seed") <- seed

    return(output_list)
}


#' zeallot destructure function
#' @keywords internal
destructure.mmBPFA.simulated.data <- function(x) {
    list(x$data, x$true_values)
}