#' Process results from an mmBPFA sampling run
#'
#' \code{process_results} pools the kept draws of all chains and returns
#' posterior summaries of the counted dimensions: columns that are switched on
#' and have at least three active items (design §4.4). Draws are relabelled
#' against a pivot draw with the median number of dimensions, and each matched
#' dimension is sign-aligned to the pivot before averaging, so label switching
#' and sign flips (within or across chains) do not cancel.
#'
#' The model's correlation uses every column; the loadings and scores reported
#' here use only the counted dimensions. "Uncounted Correlation Share" reports
#' how much of the copula correlation the uncounted columns carry.
#'
#' @param mcmc_object Output of \code{mmBPFA_sampler} (class \code{bpfa.results}).
#' @param SE Whether to add credible intervals.
#' @param p Lower and upper probabilities of the credible intervals.
#' @param to_tibble Whether to return tibbles.
#'
#' @return Object of class \code{mcmc.output.processed} with the stable
#'   dimensions, the K distribution (and R-hat for K with several chains), loadings
#'   on the sampler and copula scales, factor scores, slab variances, item
#'   intercepts, log-likelihoods, level cutpoints, the uncounted correlation share
#'   and the priors.
#'
#' @export
process_results <- function(mcmc_object, SE = TRUE, p = c(0.05, 0.95), to_tibble = TRUE) {
    if (!inherits(mcmc_object, "bpfa.results")) stop("mcmc_object must be the output of mmBPFA_sampler().")
    H <- attr(mcmc_object, "priors")$H
    n_chains <- length(mcmc_object)

    # pool draws across chains
    pool <- function(field) unlist(lapply(mcmc_object, `[[`, field), recursive = FALSE)
    lambda <- pool("lambda"); z <- pool("z"); cc <- pool("c"); omega <- pool("omega")
    alpha <- pool("alpha"); theta <- pool("theta")
    chain_of <- rep(seq_len(n_chains), sapply(mcmc_object, function(ch) length(ch$K)))

    # counted dimensions of each draw: on (c_h > h) with at least three items
    counted <- purrr::map2(cc, z, ~ which(.x > seq_len(H) & colSums(.y) >= 3))
    K <- lengths(counted)
    median_K <- stats::median(K)

    # relabel against a pivot draw with the median K, then align signs
    pivot <- which.min(abs(K - median_K))
    perms <- align_dimensions(purrr::map2(z, counted, ~ .x[, .y, drop = FALSE]), pivot = pivot)
    piv_lambda <- lambda[[pivot]][, counted[[pivot]], drop = FALSE]
    aligned <- function(mat, t, scale = NULL) {
        m <- mat[, counted[[t]], drop = FALSE][, perms[[t]], drop = FALSE]
        if (!is.null(scale)) m <- m * rep(scale, each = nrow(m))
        m
    }
    signs <- lapply(seq_along(lambda), function(t) {
        l <- aligned(lambda[[t]], t)
        s <- sign(colSums(l * piv_lambda))
        s[is.na(s) | s == 0] <- 1
        s
    })
    lambda_al <- lapply(seq_along(lambda), function(t) aligned(lambda[[t]], t, signs[[t]]))
    omega_al <- lapply(seq_along(omega), function(t) aligned(omega[[t]], t, signs[[t]]))
    copula_al <- lapply(seq_along(lambda), function(t)
        aligned(lambda[[t]] / sqrt(1 + rowSums(lambda[[t]]^2)), t, signs[[t]]))
    z_al <- lapply(seq_along(z), function(t) aligned(z[[t]], t))
    theta_al <- lapply(seq_along(theta), function(t) theta[[t]][counted[[t]]][perms[[t]]])

    stable_K <- find_stable_dims(perms)
    zero_csums <- t(sapply(z_al, function(m) colSums(m[, stable_K, drop = FALSE])))
    if (length(stable_K) == 1) zero_csums <- matrix(zero_csums, ncol = 1)
    dim_order <- order(colMeans(zero_csums, na.rm = TRUE), decreasing = TRUE)

    out_msg <- paste("The sampler found", median_K, "dimensions (median);", length(stable_K),
                     "are present in at least half of the draws and are reported.\n")
    cat(out_msg)

    # convergence of K across chains
    rhat_K <- if (n_chains > 1) rhat(K, chain_of) else NA_real_
    if (!is.na(rhat_K) && rhat_K > 1.1) {
        warning(sprintf(paste0("Chains disagree on the number of dimensions (R-hat for K = %.2f). ",
                               "Run longer chains or more chains from different K_start before relying on K."), rhat_K))
    }

    summarise <- function(draws) extract_omega_lambda_mean_se(draws, dim_order, stable_K, p, se = SE)
    alpha_mat <- do.call(cbind, alpha)
    alpha_out <- tibble::tibble(item = seq_len(nrow(alpha_mat)), alpha = rowMeans(alpha_mat))
    if (SE) {
        q <- t(apply(alpha_mat, 1, stats::quantile, p = p))
        alpha_out[[paste0("ci_", p[1])]] <- q[, 1]
        alpha_out[[paste0("ci_", p[2])]] <- q[, 2]
    }
    theta_mat <- do.call(cbind, lapply(theta_al, function(th) th[stable_K][dim_order]))
    theta_out <- tibble::tibble(dimension = seq_along(stable_K),
                                theta = rowMeans(matrix(theta_mat, nrow = length(stable_K)), na.rm = TRUE))

    log_liks <- unlist(lapply(mcmc_object, `[[`, "log_lik"))
    log_liks_out <- tibble::tibble(mcmc_iter = seq_along(log_liks), log_liks = log_liks, chain = chain_of)

    output_list <- list()
    output_list[["Stable Dimensions"]] <- stable_K
    output_list[["Median Number Dimensions"]] <- median_K
    output_list[["K Distribution"]] <- table(K)
    output_list[["R-hat K"]] <- rhat_K
    output_list[["Factor Loadings"]] <- summarise(lambda_al)
    output_list[["Copula Loadings"]] <- summarise(copula_al)
    output_list[["Factor Scores"]] <- summarise(omega_al)
    output_list[["Slab Variances"]] <- theta_out
    output_list[["Item-level Intercepts"]] <- alpha_out
    output_list[["log-Likelihoods"]] <- log_liks_out
    output_list[["Cutpoints"]] <- pooled_cutpoints(mcmc_object)
    output_list[["Uncounted Correlation Share"]] <- uncounted_share(lambda, counted)
    output_list[["Priors"]] <- attr(mcmc_object, "priors")

    if (!to_tibble) {
        for (nm in c("Factor Loadings", "Copula Loadings", "Factor Scores")) {
            output_list[[nm]] <- as.matrix(output_list[[nm]])
        }
    }
    class(output_list) <- "mcmc.output.processed"
    attr(output_list, "Include standard errors") <- SE
    attr(output_list, "Conf. Interval") <- p
    attr(output_list, "Message") <- out_msg
    output_list
}


#' zeallot destructure function
#' @keywords internal
destructure.mcmc.output.processed <- function(x) {
    unname(unclass(x))
}


#' Gelman-Rubin potential scale reduction for one scalar
#'
#' @keywords internal
#'
#' @param x draws pooled over chains.
#' @param chain chain index of each draw.
#'
#' @return R-hat; 1 when every draw is identical, Inf when chains are constant but differ.
rhat <- function(x, chain) {
    groups <- split(x, chain)
    n <- min(lengths(groups))
    groups <- lapply(groups, function(g) g[seq_len(n)])
    means <- sapply(groups, mean)
    W <- mean(sapply(groups, stats::var))
    B <- n * stats::var(means)
    if (W == 0) return(if (B == 0) 1 else Inf)
    sqrt(((n - 1) / n * W + B / n) / W)
}


#' Posterior-mean level cutpoints per item, pooled over chains
#'
#' Each chain stores the posterior mean of the thresholds between adjacent
#' observed levels, flat and column by column. Returns one vector per item with
#' its thresholds followed by Inf (the format \code{check_accuracy} and \code{gmp} use).
#'
#' @keywords internal
pooled_cutpoints <- function(mcmc_object) {
    n_th <- pmax(lengths(attr(mcmc_object, "levels")) - 1L, 0L)
    th <- Reduce(`+`, lapply(mcmc_object, `[[`, "thresholds")) / length(mcmc_object)
    split_th <- split(th, rep(seq_along(n_th), n_th))
    lapply(seq_along(n_th), function(j) c(if (n_th[j] > 0) split_th[[as.character(j)]] else numeric(0), Inf))
}


#' Share of the copula correlation carried by uncounted columns
#'
#' Mean absolute difference between the copula correlation of all columns and
#' that of the counted columns, averaged over up to 50 evenly spaced draws and,
#' for more than 1000 items, over 1000 randomly chosen items.
#'
#' @keywords internal
uncounted_share <- function(lambda, counted) {
    cor_of <- function(L) { C <- tcrossprod(L) + diag(nrow(L)); D <- 1 / sqrt(diag(C)); C * outer(D, D) }
    draws <- unique(round(seq(1, length(lambda), length.out = min(50, length(lambda)))))
    items <- seq_len(nrow(lambda[[1]]))
    if (length(items) > 1000) items <- sort(sample(items, 1000))
    off <- upper.tri(diag(length(items)))
    mean(sapply(draws, function(t) {
        L <- lambda[[t]][items, , drop = FALSE]
        mean(abs(cor_of(L) - cor_of(L[, counted[[t]], drop = FALSE]))[off])
    }))
}


#' Function to examine posterior log-likelihoods
#'
#' @param results processed results from \code{process_results}.
#' @param probs quantiles of the log-likelihood distribution to display.
#' @param plot whether to plot the log-likelihood density.
#'
#' @export
check_logliks <- function(results, probs = c(0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95), plot = TRUE) {
    if (!inherits(results, "mcmc.output.processed")) stop("Results must be object of class 'mcmc.output.processed'")
    log_liks <- results[["log-Likelihoods"]]$log_liks
    mean_log_liks <- round(mean(log_liks), 3)
    quantile_log_liks <- stats::quantile(log_liks, p = probs)
    cat("Mean log-likelihood: \n", mean_log_liks, "\n")
    cat("Log-likelihood Quantiles: \n")
    print(quantile_log_liks)
    if (plot) {
        plot(stats::density(log_liks), main = "Log-Likelihood Distribution", xlab = "Draws")
        abline(v = stats::median(log_liks), lty = 2)
    } else {
        list("mean_log_liks" = mean_log_liks, "log_lik_quantiles" = quantile_log_liks)
    }
}


#' Function to check reconstruction accuracy
#'
#' \code{check_accuracy} predicts each observed cell from the posterior-mean
#' loadings, scores, intercepts and level cutpoints, and reports exact-match
#' accuracy, RMSE, MAE and the geometric mean probability of the observed data.
#' If \code{true_vals} from \code{generate_sim_data} is supplied, it also reports
#' RMSE and MAE of the latent quantities. A confusion matrix, precision, recall
#' and Kappa are available when every item has the same set of observed values.
#'
#' @param dat Data matrix.
#' @param results processed results from \code{process_results}, without standard errors.
#' @param true_vals Optional list of true latent values from \code{generate_sim_data}.
#' @param confusion Whether to return a confusion matrix.
#' @param return_preds Whether to also return the predicted data matrix.
#'
#' @return List of accuracy metrics.
#'
#' @export
check_accuracy <- function(dat, results, true_vals = NULL, confusion = TRUE, return_preds = FALSE) {
    if (!inherits(results, "mcmc.output.processed")) stop("Results must be of class `mcmc.output.processed`.")
    if (isTRUE(attributes(results)$`Include standard errors`)) {
        stop("Processed results must not include standard errors. Please re-run process_results() with SE set to FALSE.")
    }
    stable_K <- results[["Stable Dimensions"]]
    lambda <- as.matrix(results[["Factor Loadings"]][, -1])
    omega <- as.matrix(results[["Factor Scores"]][, -1])
    alpha <- results[["Item-level Intercepts"]][[2]]
    x <- omega %*% t(lambda) - matrix(alpha, nrow(omega), length(alpha), byrow = TRUE)

    cutpoints <- results[["Cutpoints"]]
    levels <- lapply(seq_len(ncol(dat)), function(j) sort(unique(dat[!is.na(dat[, j]), j])))
    preds <- sapply(seq_len(ncol(dat)), function(j) {
        if (length(levels[[j]]) == 0) return(rep(NA_real_, nrow(dat)))
        levels[[j]][findInterval(x[, j], utils::head(cutpoints[[j]], -1)) + 1L]
    })

    output_list <- list()
    acc <- rmse_mae(actual = dat, predicted = preds)
    output_list[["RMSE"]] <- acc$rmse
    output_list[["Accuracy"]] <- sum(dat == preds, na.rm = TRUE) / sum(!is.na(dat))
    output_list[["MAE"]] <- acc$mae
    output_list[["GMP"]] <- gmp(dat = dat, pred = x, cutpoints = cutpoints)

    if (!is.null(true_vals)) {
        c(true_x, true_lambda, true_zeros, true_omega, true_alpha, true_pi_k, true_gamma_k) %<-% true_vals
        K <- min(ncol(true_lambda), length(stable_K))
        theta <- results[["Slab Variances"]][[2]][1:K]
        rmse_mae_results <- tibble::tibble(
            "parameters" = c("x", "lambda", "omega", "alpha", "gamma_k"),
            "actual" = list(true_x, true_lambda[, 1:K], true_omega[, 1:K], true_alpha, true_gamma_k[1:K]),
            "predicted" = list(x, lambda[, 1:K], omega[, 1:K], alpha, 1 / theta)
        )
        rmse_mae_results <- rmse_mae_results %>%
            dplyr::mutate(results = purrr::map2(.data$actual, .data$predicted, ~ rmse_mae(.x, .y))) %>%
            tidyr::unnest_wider(results) %>%
            dplyr::select("parameters", "rmse", "mae")
        output_list[["Accuracy Metrics"]] <- rmse_mae_results
    }

    same_levels <- length(unique(lapply(levels, as.character))) == 1 && length(levels[[1]]) <= 20
    if (confusion && !same_levels) {
        warning("A confusion matrix needs every item to have the same observed values. Returning results without it.")
    } else if (confusion) {
        values <- levels[[1]]
        confusion_mat <- table(factor(preds, levels = values), factor(dat, levels = values),
                               dnn = c("Predicted", "Actual"))
        diags <- diag(confusion_mat)
        total <- sum(confusion_mat)
        p_actual <- colSums(confusion_mat) / total; p_pred <- rowSums(confusion_mat) / total
        chance <- sum(p_actual * p_pred)
        output_list[["Confusion Matrix"]] <- confusion_mat
        output_list[["Class-wise Precision/Recall"]] <- tibble::tibble(
            values = values, precision = diags / rowSums(confusion_mat), recall = diags / colSums(confusion_mat))
        output_list[["Kappa Score"]] <- (sum(diags) / total - chance) / (1 - chance)
    }
    if (return_preds) output_list[["Predictions"]] <- preds
    output_list
}


#' Posterior means (and credible intervals) of aligned loadings or scores
#'
#' @keywords internal
extract_omega_lambda_mean_se <- function(x, dim_order, stable_K, p, se = TRUE, to_tibble = TRUE) {
    x_array <- abind::abind(purrr::map(x, ~ .x[, stable_K, drop = FALSE]), along = 3)
    x_array <- x_array[, dim_order, , drop = FALSE]
    x_mean <- apply(x_array, 1:2, mean, na.rm = TRUE)
    if (!is.matrix(x_mean)) x_mean <- matrix(x_mean, ncol = length(stable_K))
    if (!se) {
        x_out <- tibble::as_tibble(x_mean, .name_repair = "minimal")
        names(x_out) <- paste0("dim_", seq_along(stable_K))
        x_out <- tibble::add_column(x_out, row = seq_len(nrow(x_mean)), .before = 1)
        return(x_out)
    }
    lo <- apply(x_array, 1:2, stats::quantile, p = p[1], na.rm = TRUE)
    hi <- apply(x_array, 1:2, stats::quantile, p = p[2], na.rm = TRUE)
    x_out <- tibble::tibble(
        row = rep(seq_len(nrow(x_mean)), length(stable_K)),
        dimension = rep(seq_along(stable_K), each = nrow(x_mean)),
        mean = as.vector(x_mean), lo = as.vector(lo), hi = as.vector(hi)
    )
    names(x_out)[4:5] <- paste0("ci_", p)
    x_out
}


#' match each draw's dimensions to those of a pivot draw
#'
#' Dimensions are paired greedily by the Euclidean distance between their
#' sparsity patterns, closest pairs first.
#'
#' @keywords internal
#'
#' @param zeros list of P x K binary matrices, one per draw.
#' @param pivot index of the reference draw.
#'
#' @return list of integer vectors, one per draw, giving for each pivot dimension
#' the matching column of that draw (NA when the draw has no dimension left to match).
align_dimensions <- function(zeros, pivot) {
  ref <- zeros[[pivot]]
  purrr::map(zeros, function(z) {
    perm <- rep(NA_integer_, ncol(ref))
    if (ncol(z) == 0 || ncol(ref) == 0) return(perm)
    dist <- rdist::cdist(t(ref), t(z), "euclidean")
    for (step in seq_len(min(ncol(ref), ncol(z)))) {
      best <- which(dist == min(dist), arr.ind = TRUE)[1, ]
      perm[best[1]] <- best[2]
      dist[best[1], ] <- Inf
      dist[, best[2]] <- Inf
    }
    perm
  })
}


#' function to identify stable dimensions from sampling run
#'
#' @keywords internal
#'
#' @param perms list of pivot-to-draw matchings from \code{align_dimensions}.
#' @param threshold minimum share of draws in which a pivot dimension must be present.
#'
#' @return indices of the stable pivot dimensions.
find_stable_dims <- function(perms, threshold = 0.5) {
  if (length(perms[[1]]) == 0) return(integer(0))
  matched <- do.call(rbind, purrr::map(perms, ~ !is.na(.x)))
  which(colMeans(matched) >= threshold)
}


#' function to calculate rmse and mae
#'
#' @keywords internal
rmse_mae <- function(actual, predicted) {
  rmse <- sqrt(mean((actual - predicted)^2, na.rm = T))
  mae <- mean(abs(actual - predicted), na.rm = T)
  return(list("rmse" = rmse, "mae" = mae))
}

#' function to calculate geometric mean probability
#'
#' @keywords internal
gmp <- function(dat, pred, cutpoints) {
  # geometric mean probability of the observed cells: exp(mean log probability)
  log_probs <- cell_log_probs(pred, dat, cutpoints)
  exp(mean(log_probs, na.rm = TRUE))
}
