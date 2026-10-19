#' mmBPFA Gibbs sampler
#'
#' \code{mmBPFA_sampler} fits the Gaussian-copula factor model with the
#' extended rank likelihood (Hoff 2007), so any mix of binary, ordinal, count
#' and continuous items can be analysed without modelling their margins. The
#' number of dimensions is learned with a cumulative shrinkage prior over
#' columns (Legramanti, Durante & Dunson 2020), and each item's inclusion on
#' each dimension is learned within columns (Schiavon, Canale & Dunson 2022).
#' Only observed cells are used: missing cells are treated as missing at random
#' and are never imputed, so the cost of a sweep scales with the number of
#' observed cells.
#'
#' @param dat Numeric data matrix (N x P); NA marks a missing cell.
#' @param priors A preset name or a list of prior settings; see \code{\link{resolve_priors}}.
#' @param H Truncation: maximum number of columns. Default min(floor((P - 1) / 2), 50).
#' @param K_start Number of columns switched on at the start. Capped at H.
#' @param warmup Number of warmup iterations.
#' @param iter Number of sampling iterations after warmup.
#' @param thin Keep every \code{thin}-th sampling iteration.
#' @param chains Number of chains.
#' @param n_cores Number of chains run in parallel.
#' @param threads Threads used inside each chain (OpenMP). Results do not depend on it.
#' @param seed Random seed; chain c uses \code{seed + c - 1}.
#' @param max_memory_gb Stop before sampling if the kept draws would need more memory than this.
#'
#' @return Object of class \code{bpfa.results}: a list with one element per chain,
#'   each holding the kept draws (\code{lambda}, \code{z}, \code{c}, \code{omega},
#'   \code{alpha}, \code{theta}, \code{a_cusp}), the number of counted dimensions
#'   \code{K} and the log-likelihood \code{log_lik} per draw, and posterior-mean level
#'   \code{thresholds}. Attributes record the priors (see \code{\link{get_priors}})
#'   and the run settings.
#'
#' @export
mmBPFA_sampler <- function(
    dat,
    priors = "uniform",
    H = NULL,
    K_start = 10,
    warmup = 1000,
    iter = 5000,
    thin = 5,
    chains = 1,
    n_cores = 1,
    threads = 1,
    seed = 101473,
    max_memory_gb = 4
    ) {
    if (warmup < 0 || iter < 1 || thin < 1 || chains < 1) {
        stop("warmup must be >= 0, and iter, thin and chains must be positive")
    }
    data <- prepare_data(dat)
    prior <- resolve_priors(priors, data$p, H)
    if (K_start > prior$H) {
        message(sprintf("K_start (%d) exceeds H (%d); using %d", as.integer(K_start), prior$H, prior$H))
        K_start <- prior$H
    }

    # kept draws store lambda and z (P x H), omega (N x H), alpha (P) and column vectors
    H <- prior$H
    gb <- 8 * chains * floor(iter / thin) * (2 * data$p * H + data$n * H + data$p + 3 * H) / 1e9
    if (gb > max_memory_gb) {
        stop(sprintf(paste0("Stored draws would need about %.2f GB (limit %.2f GB); ",
                            "increase `thin`, reduce `iter`, or raise `max_memory_gb`."), gb, max_memory_gb))
    }

    run <- function(chain) run_one_chain(data, prior, K_start, warmup, iter, thin, seed + chain - 1, threads)
    chain_out <- if (chains == 1 || n_cores == 1) {
        lapply(seq_len(chains), run)
    } else if (.Platform$OS.type == "windows") {
        run_chains_psock(chains, n_cores, run)
    } else {
        parallel::mclapply(seq_len(chains), run, mc.cores = n_cores)
    }

    structure(
        chain_out,
        class = "bpfa.results",
        priors = prior,
        n = data$n, p = data$p, levels = level_values(data),
        warmup = warmup, iter = iter, thin = thin, chains = chains,
        K_start = K_start, seed = seed
    )
}


#' Run one chain from a seeded starting point
#'
#' @keywords internal
run_one_chain <- function(data, prior, K_start, warmup, iter, thin, seed, threads = 1) {
    set.seed(seed)
    H <- prior$H; p <- data$p; n <- data$n
    on <- seq_len(H) <= K_start
    z <- matrix(0L, p, H)
    z[, on] <- rbinom(p * sum(on), 1, 0.5)
    state <- list(
        lambda = z * matrix(rnorm(p * H), p, H),
        z = z,
        rho = rep(0.5, H),
        c = as.integer(ifelse(on, H, 1L)),
        v = c(rep(0.3, H - 1), 1),
        a_cusp = 1,
        theta = ifelse(on, 1, prior$spike_var)
    )
    omega <- matrix(rnorm(n * H), n, H)
    .run_chain(data, initial_x(data), state, omega, rep(0, p), prior, warmup, iter, thin, seed, threads)
}


#' Distinct observed values of each item, in increasing order
#'
#' @keywords internal
level_values <- function(data) {
    lapply(seq_len(data$p), function(j) {
        cells <- seq_len(data$col_ptr[j + 1] - data$col_ptr[j]) + data$col_ptr[j]
        unique(data$value[cells])
    })
}


#' Run chains on a socket cluster
#'
#' Used for parallel chains on Windows, where forking is unavailable.
#'
#' @keywords internal
#'
#' @param chains Number of chains.
#' @param n_cores Number of worker processes.
#' @param chain_fun Function of the chain index returning one chain's output.
#'
#' @return List with one element per chain.
run_chains_psock <- function(chains, n_cores, chain_fun) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(mmBPFA))
    parallel::parLapply(cl, seq_len(chains), chain_fun)
}
