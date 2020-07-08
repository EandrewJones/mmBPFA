#' Process results from mmBPFA sampling run
#' 
#' \code{process_results} ingests raw results from \code{mmBPFA_sampler} and returns processed estimates of posterior quantities of interest.
#' 
#' @param mcmc_object An S3 object of class `bpfa.results` from mmBPFA sampling run.
#' @param SE A boolean indicating whether to calculate standard error credibility intervals around posterior estimates.
#' @param p a vector indicating the lower and upper probability bounds for the credibility intervals. Default is c(0.05, 0.95).
#' @param to_tibble A boolean indicating whether to convert output to tibbles for use in tidy workflow. Default is TRUE.
#' 
#' @return An S3 object of class `mcmc.output.processed` containing the processed posterior margins.
#' 
#' @export
#' 
process_results <- function(
    mcmc_object,
    SE = T,
    p = c(0.05, 0.95),
    to_tibble = T
    ) {
    # Assertions
    if (class(mcmc_object) != "bpfa.results") stop("mcmc_object must be the output of bpfa sampling run.")

    # extract results
    n_dims <- mcmc_object[[1]]
    log_liks <- mcmc_object[[2]]
    results <- mcmc_object[[3]]

    # get mode
    mode <- attributes(mcmc_object)$mode

    # Process n_dims
    median_K <- median(n_dims)
    stable_K <- min(n_dims)

    out_msg <- paste(
        "The sampler discovered", median_K, "dimensions (median),", stable_K, "of which are stable dimensions. Processing results for first", stable_K, "stable dimensions.\n"
    )
    cat(out_msg)

    # process log_liks
    if (to_tibble) {
        log_liks_out <- tibble::tibble(
            "mcmc_iter" = 1:length(log_liks),
            "log_liks" = log_liks
        )
    } else {
        log_liks_out <- log_liks
    }

    # Process results
    c(zeros, lambda, omega, alpha, gamma_k, cutpoints) %<-% results

    # 1) zeros
    zero_csums_array <- t(abind(map(zeros, ~ .x[, 1:stable_K] %>% colSums()), along = 2))
    zero_csums_mean <- colMeans(zero_csums_array)
    dim_order <- rev(order(zero_csums_mean))

    # 2) lambda (matrix)
    lambda_out <- extract_omega_lambda_mean_se(
        x = lambda,
        dim_order = dim_order,
        stable_K = stable_K,
        p = p,
        se = SE
    )

    # 3) omega (matrix)
    omega_out <- extract_omega_lambda_mean_se(
        x = omega,
        dim_order = dim_order,
        stable_K = stable_K,
        p = p,
        se = SE
    )

    # 4) process alpha
    alpha_mat <- do.call(cbind, alpha)
    alpha_mean <- rowMeans(alpha_mat)
    alpha_out <- round(alpha_mean, 7)
    if (SE) {
        se_mat <- apply(alpha_mat, 1, quantile, p = p)
        se_long <- t(se_mat)
        alpha_out <- cbind(1:length(alpha_mean), alpha_mean, se_long)
        names(alpha_out) <- c("row", "mean", paste0("ci_", p[1]), paste0("ci_", p[2]))
        alpha_out <- round(alpha_out, 7)
        alpha_out <- tibble::as_tibble(alpha_out, .name_repair = "minimal")
    }
    if (to_tibble & !SE) {
        alpha_out <- tibble::tibble(
            "item" = 1:length(alpha_out),
            "alpha" = as.vector(alpha_out)
        )
    }

    # 4) process gamma_k
    gamma_k_mat <- do.call(cbind, map(gamma_k, ~ .x[1:stable_K]))
    gamma_k_mat <- gamma_k_mat[dim_order, ]
    gamma_k_mean <- rowMeans(gamma_k_mat)
    gamma_k_out <- round(gamma_k_mean, 7)
    if (SE) {
        se_mat <- apply(gamma_k_mat, 1, quantile, p = p)
        se_long <- t(se_mat)
        gamma_k_out <- cbind(1:length(gamma_k_mean), gamma_k_mean, se_long)
        names(gamma_k_out) <- c("row", "mean", paste0("ci_", p[1]), paste0("ci_", p[2]))
        gamma_k_out <- round(gamma_k_out, 7)
        gamma_k_out <- tibble::as_tibble(gamma_k_out, .name_repair = "minimal")
    }
    if (to_tibble & !SE) {
        gamma_k_out <- tibble::tibble(
            "dimension" = 1:length(gamma_k_out),
            "gamma_k" = gamma_k_out
        )
    }

    # 5) process cutpoints
    c(a_mean, b_mean) %<-% map(
        c("a", "b"),
        function(x) {
            cutpoint_array <- abind(map(cutpoints, ~ .x[[x]]), along = 3)
            cutpoint_mean <- apply(cutpoint_array, 1:2, mean, simplify = FALSE)
            cutpoint_mean <- apply(cutpoint_mean, 2, unique)
            if (mode == "multi") {
                cutpoint_mean
            } else {
                as_tibble(cutpoint_mean, .name_repair = "minimal")
            }
        }
    )
    a_sort <- map(a_mean, ~ sort(.x)[-1])
    b_sort <- map(b_mean, ~ sort(.x)[-length(.x)])
    cutpoints <- map2(a_sort, b_sort, ~ c(colMeans(rbind(.x, .y)), Inf))

    # store outputs
    output_list <- list()
    output_list[["Stable Dimensions"]] <- stable_K
    output_list[["Median Number Dimensions"]] <- median_K
    output_list[["Factor Loadings"]] <- lambda_out
    output_list[["Factor Scores"]] <- omega_out
    output_list[["Factor Precisions"]] <- gamma_k_out
    output_list[["Item-level Intercepts"]] <- alpha_out
    output_list[["log-Likelihoods"]] <- log_liks_out
    output_list[["Cutpoints"]] <- cutpoints

    class(output_list) <- "mcmc.output.processed"
    attr(output_list, "Include standard errors") <- SE
    attr(output_list, "Conf. Interval") <- p
    attr(output_list, "Message") <- out_msg

    return(output_list)
}


# zeallot destructure function
destructure.mcmc.output.processed <- function(x) {
    list(
        x[[1]], x[[2]], x[[3]], x[[4]], x[[5]], x[[6]], x[[7]], x[[8]]
    )
}


#' Function to process benchmarks
#' 
#' \code{process_benchmarks} ingests the raw results from \code{mmBPFA_sampler} and 
#' returns the raw benchmarks and average time per sampling step. Only works if the sampler was run with
#' \code{benchmarks = TRUE}.
#'
#' @param mcmc_object A S3 object of class `bpfa.results` containing benchmarks.
#' 
#' @return A list containing the raw benchmarks and average timing per sample step.
#' 
#' @export
# TODO #6 update functionality for multi-chain outputs?
process_benchmarks <- function(mcmc_object) {
    # Assertions

    if (!attributes(mcmc_object)$benchmark) {
        stop("mcmc_object must include timing Benchmarks for each sampling step. Please re-run the bpfa sampling with benchmark = T")
    }

    # get benchmarks
    benchmarks <- mcmc_object[["Benchmarks"]]

    # reshape raw benchmarks
    raw_times_out <- benchmarks %>%
        bind_rows(.id = "phase") %>%
        unnest(time) %>%
        mutate(
            time = map(time, ~ .x %>%
            unlist() %>%
            t() %>%
            as_tibble(.name_repair = "minimal")),
            phase = factor(phase, levels = c("pre_sparse", "warmup", "mcmc"))
        ) %>%
        unnest(time) %>%
        mutate_at(
            .vars = vars(tic.elapsed, toc.elapsed),
            .funs = as.numeric
        ) %>%
        mutate(
            time = toc.elapsed - tic.elapsed,
        ) %>%
        group_by(
            phase
        ) %>%
        mutate(
            phase_iter = rep(1:(n() / 6), each = 6)
        ) %>%
        ungroup() %>%
        select(
            phase, phase_iter, sample_step = msg, time
        ) %>%
        pivot_wider(
            names_from = sample_step,
            values_from = time
        ) %>%
        rowid_to_column(
            "iter"
        )

    # calculate average times
    avg_times_out <- raw_times_out %>%
        group_by(phase) %>%
        summarize_at(
            .vars = vars(3:8),
            .funs = mean
        )

    # store outputs
    output_list <- list()
    output_list[["Raw Benchmarks"]] <- raw_times_out
    output_list[["Sampling Step Averages"]] <- avg_times_out

    return(output_list)
}


#' Function to examine posterior log-Likelihoods
#' 
#' \code{check_logliks} ingests the raw vector of posterior log-likelihoods and returns quantiles of interests, the mean, and
#' a plot of the distribution if \code{plot = TRUE}. A well-behaving posterior distribution should be unimodal and peak at the mean.
#'
#' @param log_liks A numeric vector of posterior draws from mmBPFA sampler. Must not be processed.
#' @param probs A numeric vector of quantiles to display for log-likelihoods distribution. Default is c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95).
#' @param plot A boolean indicating whether to plot the log-likelihood distribution. Default is TRUE.
#' 
#' @return If plot is TRUE, returns none, else returns mean log-likelihood and quantiles.
#' 
#' @export
#' 
check_logliks <- function(
    log_liks,
    probs = c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95),
    plot = TRUE
    ) {

    # Assertions
    if(!is.vector(log_liks)) stop("log_liks must be a vector of draws!")

    # get measures of central tendency
    mean_log_liks <- round(mean(log_liks), 3)
    quantile_log_liks <- quantile(log_liks, p = probs)

    cat("Mean log-likelihood: \n", mean_log_liks, "\n")
    cat("Log-likelihood Quantiles: \n")
    print(quantile_log_liks)

    # plot
    if (plot) {
        plot(
            density(log_liks),
            main = "Log-Likelihood Distribution",
            xlab = "Draws"
        )
        abline(v = median(log_liks), lty = 2)
    } else {
        output_list <- list(
            "mean_log_liks" = mean_log_liks,
            "log_lik_quantiles" = quantile_log_liks
        )
        return(output_list)
    }
}

# TODO Add function to check acf


#' Function to check reconstruction accuracy
#' 
#' \code{check_accuracy} ingests the cleaned results from \code{process_results} without standard errors and 
#' returns the accuracy in reconstructing the observed data. If \code{true_vals} from \code{generate_sim_data} 
#' is supplied, the function will also return the RMSE and MAE of posterior predictions for the latent values. 
#' A confusion matrix, precision, recall, and Kappa score are also available via \code{confusion = TRUE} for data
#' with fixed marginal distributions. 
#' 
#' @param dat Data matrix.
#' @param results processed results from \code{process_results}. Must not include standard errors.
#' @param true_vals Optional list of true latent values from \code{generate_sim_data}.
#' @param mode Margin type. Must be one of "fixed", "multi" or "mixed".
#' @param confusion Boolean indicating whether to return a confusion matrix. Only works for "fixed" margin data.
#' 
#' @return List of observed accuracy and accuracy metrics (rmse and mae). If confusion is TRUE, also returns confusion matrix results.
#' 
#' @export
#' 
check_accuracy <- function(
    dat,
    results,
    true_vals = NULL,
    mode,
    confusion = T
    ) {
    # Assertions
    class_check <- class(results) == "mcmc.output.processed"
    if (!class_check) stop("Results must be of class `mcmc.output.process`.")
    compat_check <- attributes(results)$`Include standard errors`
    if (compat_check) stop("Processed results must not include standard errors. Please re-run process_results() with SE set to FALSE.")

    # get results
    stable_K <- results[["Stable Dimensions"]]
    lambda <- results[["Factor Loadings"]]
    lambda <- as.matrix(lambda[, 2:ncol(lambda)])
    omega <- results[["Factor Scores"]]
    omega <- as.matrix(omega[, 2:ncol(omega)])
    alpha <- results[["Item-level Intercepts"]]
    alpha <- as.matrix(alpha[,2])
    alpha_mat <- matrix(
        nrow = nrow(omega),
        ncol = nrow(lambda),
        rep(alpha, nrow(omega)),
        byrow = TRUE
    )

    # calculate X
    x <- t(lambda %*% t(omega)) - alpha_mat
    x_tbl <- as_tibble(x, .name_repair = "minimal")

    # predict y
    cutpoints <- results[["Cutpoints"]]
    c(..., n_levels, margin_vals) %<-% calculate_d_mask(dat = dat, mode = mode)
    if (mode == "fixed") {
        preds <- pmap_dfr(
            list(x_tbl, cutpoints),
            function(x, y) {
                map_dbl(x, ~ margin_vals[which.max(.x < y)])
            }
        )
    } else {
        preds <- pmap_dfr(
            list(x_tbl, margin_vals, cutpoints),
            function(x, y, z) {
                map_dbl(x, ~ y[which.max(.x < z)])
            }
    )
    }
    preds <- as.matrix(preds)

    # check reconstruction accuracy
    if (mode == "mixed") {
        accuracy <- sqrt(mean((preds - dat)^2, na.rm = T))
    } else {
        accuracy <- sum(dat == preds, na.rm = T) / length(dat[!is.na(dat)])
    }

    # store output
    output_list <- list()
    output_list[["Observed Accuracy"]] <- accuracy

    # check accuracy against true vals
    if (!is.null(true_vals)) {
        # function to calculate rmse and mae
        rmse_mae <- function(actual, predicted) {
            rmse <- sqrt(mean((actual - predicted)^2))
            mae <- mean(abs(actual - predicted))
            return(list("rmse" = rmse, "mae" = mae))
        }

        # extract true parameters
        c(
            true_x, true_lambda, true_zeros,
            true_omega, true_alpha, true_pi_k, true_gamma_k
        ) %<-% true_vals
        true_K <- ncol(true_lambda)

        # get correct dimensions from sampled parameters
        K <- min(true_K, stable_K)
        lambda <- lambda[, 1:K]
        true_lambda <- true_lambda[, 1:K]
        omega <- omega[, 1:K]
        true_omega <- true_omega[, 1:K]
        gamma_k <- results[["Factor Precisions"]][[2]][1:K]
        true_gamma_k <- true_gamma_k[1:K]

        # calculate metrics
        rmse_mae_results <- tibble:::tibble(
            "parameters" = c("x", "lambda", "omega", "alpha", "gamma_k"),
            "actual" = list(true_x, true_lambda, true_omega, true_alpha, true_gamma_k),
            "predicted" = list(x, lambda, omega, alpha, gamma_k)
        )
        rmse_mae_results <- rmse_mae_results %>%
            mutate(results = map2(actual, predicted, ~ rmse_mae(.x, .y))) %>%
            unnest_wider(results) %>%
            select(parameters, rmse, mae)
        output_list[["Accuracy Metrics"]] <- rmse_mae_results
    }

    # produce confusion matrix
    if (confusion) {
        if (mode != "fixed") stop("Cannot produce confusion matrix for multi/mixed-modal data.")

        # confusion matrix
        values <- sort(unique(as.vector(dat)))
        confusion_mat <- matrix(0, nrow = length(values), ncol = length(values))
        rownames(confusion_mat) <- colnames(confusion_mat) <- values
        for (i in 1:length(values)) {
            for (j in 1:length(values)) {
                confusion_mat[i, j] <- sum(which(preds == values[i]) %in% which(dat == values[j]))
            }
        }
        names(dimnames(confusion_mat)) <- c("Predicted", "Actual")

        # precision and recall
        diags  <- diag(confusion_mat)
        precision <- diags / rowSums(confusion_mat)
        recall <- diags / colSums(confusion_mat)

        pc_df <- tibble::as_tibble(
            cbind(values, precision, recall),
            .name_repair = "minimal"
        )

        # Cohen's Kappa
        agree <- sum(diags) / length(preds)
        p_actual <- colSums(confusion_mat) / length(preds)
        p_pred <- rowSums(confusion_mat) / length(preds)
        chance <- sum(p_actual * p_pred)
        kappa <- (agree - chance) / (1 - chance)

        # store outputs
        output_list[["Confusion Matrix"]] <- confusion_mat
        output_list[["Class-wise Precision/Recall"]] <- pc_df
        output_list[["Kappa Score"]] <- kappa
    }

    return(output_list)
}

# internal function to estimate factor means and se"s
extract_omega_lambda_mean_se <- function(
    x,
    dim_order,
    stable_K,
    p,
    se = TRUE,
    to_tibble = TRUE
    ) {
    # reshape samples
    x_array <- abind(map(x, ~ .x[, 1:stable_K]), along = 3)
    x_array <- x_array[, dim_order, ]
    x_mean <- apply(x_array, 1:2, mean, simplify = F)

    # no se output
    x_out <- tibble::as_tibble(x_mean, .name_repair = "minimal")
    names(x_out) <- paste0("dim_", 1:stable_K)

    # Add standard errors
    if (se) {
        se_array <- apply(x_array, 1:2, quantile, p = p, simplify = F)
        se_perm <- aperm(se_array, c(2, 3, 1))
        se_long <- apply(se_perm, 3, matrix, nrow = prod(dim(se_perm)[1:2]), ncol = 1)
        x_mean_long <- matrix(x_mean, nrow = prod(dim(x_mean)), ncol = 1)
        x_out <- tibble::as_tibble(
            cbind(x_mean_long, se_long),
            .name_repair = "minimal"
        )
        names(x_out) <- c("mean", paste0("ci_", p[1]), paste0("ci_", p[2]))
        x_out$dimension <- rep(1:stable_K, each = nrow(x_mean))
        x_out <- select(x_out, dimension, everything())
    }

    # add row (feature or observation) numbers and round
    if (se) {
        x_out$row <- rep(1:nrow(x_mean), stable_K)
    } else {
        x_out$row <- 1:nrow(x_out)
    }

    x_out <- select(x_out, row, everything())
    x_out <- round(x_out, 7)

    if (!to_tibble) {
        x_out <- as.matrix(x_out)
    }

    # output
    return(x_out)
}