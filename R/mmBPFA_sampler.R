#' mmBPFA Gibbs sampler function
#' 
#' \code{mmBPFA_sampler} is the primary functional interface for the mmBPFA package. 
#' The sampler estimates a latent posterior distribution from the data using a Gibbs sampling procedure. 
#' Each conditional sampling step is written in C++ for improved sampling speed.
#' 
#' @param dat Data matrix.
#' @param mode Sets the type of margins for the simulated data. Must be one of "fixed", "multi" or "mixed." Default is "fixed."
#' @param share_noise Boolean indicating whether to sample the d hyperparameter for the factor precision draws. Setting to TRUE allows for more accurate dimension-varying precisions but takes longer to mix. Default is FALSE.
#' @param K_start Number of dimensions K to initialize sampler with. Default is 100.
#' @param pre_sparse Number of iterations sampler should take prior to turning on the beta-process prior. A minimum of 100 is suggested as a default to improve mixing and avoid the posterior dimensionality from collapsing to zero.
#' @param warmup Number of burnin iterations. Default 1000.
#' @param iter Number of posterior draws. 1000 draws usually suffices if runtime is an issue. Default is 5000.
#' @param thin Indicates how often posterior draws should be stored. Simulations suggest autocorrelation is present up to 5 lags, but tapers off at around 3. Defaults to 5 to minimize autocorrelation.
#' @param chains Number of chains to run. Implementation currently only allows for a single chain. Default is 1.
#' @param parallel Boolean indicating whether to use parallel processing for multiple chain runs.
#' @param n_cores Number of cores to use for parallel chain runs.
#' @param seed Random seed.
#' @param benchmark Boolean indicating whether to return benchmark timings of each sampling step. Primarily used for development purposes if user wants to further optimize the sampler for improved runtime. Default is FALSE. 
#' 
#' @return S3 object of class `bpfa.results` containing log-likelihoods, number of sampled dimensions at each step, and the posterior draws.
#' 
#' @export
#' 
mmBPFA_sampler <- function(
    dat,
    mode,
    share_noise = FALSE,
    K_start = 100,
    pre_sparse = 100,
    warmup = 1000,
    iter = 5000,
    thin = 5,
    chains = 1,
    parallel = TRUE,
    n_cores = 1,
    seed = 101473,
    benchmark = FALSE
    ) {
    set.seed(seed)

    # assertions
    #if (is.null(dat) == TRUE | is.matrix(data) == FALSE) stop("Data must be a numeric matrix or dataframe")
    if (prod(apply(dat, 2, class) == "numeric" | apply(dat, 2, class) == "integer") != 1) {
          stop("All elements in the data must be numeric or integers")
    }
    if (warmup < 1 | iter < 1 | pre_sparse < 1 | chains < 1) {
        stop("Warmup, iter, and pre_sparse and chains must all be positive integers")
    }

    # get dimensions
    c(n, p) %<-% get_dims_2d(x = dat, warn = TRUE)
    if (p > n) stop("Data matrix is rank deficient (P > N). Please re-check data.")

    # Check K_start against dimensions
    K_start <- check_k_start(K_start = K_start, p = p, n = n)

    # calculate D mask, margin cardinality, and unique margin values
    c(d_mask, n_levels, margin_vals) %<-% calculate_d_mask(dat = dat, mode = mode)

    # generate initial values
    c(x, lambda, zeros, omega, alpha, gamma_k, d, dc, ibp_a, ibp_b) %<-% initialize_sampler(
        dat = dat,
        mode = mode,
        n_levels = n_levels,
        margin_vals = margin_vals,
        K_start = K_start
    )

    # Single chain sampler
    if (chains == 1) {
        chain_out <- chain_call(
            pre_sparse = pre_sparse,
            warmup = warmup,
            iter = iter,
            thin = thin,
            dat = dat,
            mode = mode,
            share_noise = share_noise,
            d_mask = d_mask,
            n_levels = n_levels,
            margin_vals = margin_vals,
            x = x,
            lambda = lambda,
            zeros = zeros,
            omega = omega,
            alpha = alpha,
            gamma_k = gamma_k,
            d = d,
            dc = dc,
            ibp_a = ibp_a,
            ibp_b = ibp_b,
            benchmark = benchmark
            )
    } else if (!parallel & chains > 1) {
        chain_out <- list()
        for (i in 1:chains) {
            chain_out[[i]] <- chain_call(
                pre_sparse = pre_sparse,
                warmup = warmup,
                iter = iter,
                thin = thin,
                dat = dat,
                mode = mode,
                share_noise = share_noise,
                d_mask = d_mask,
                n_levels = n_levels,
                margin_vals = margin_vals,
                x = x,
                lambda = lambda,
                zeros = zeros,
                omega = omega,
                alpha = alpha,
                gamma_k = gamma_k,
                d = d,
                dc = dc,
                ibp_a = ibp_a,
                ibp_b = ibp_b,
                benchmark = benchmark
            )
        }
    } else {
        # Determine parallel implementation based on OS
        sys_name <- unname(Sys.info()[1])
        
        # Use snow functionality for windows
        if (grepl("windows", sys_name)) {
            # make cluster
            cl <- parallel::makeCluster(n_cores)
            parallel::clusterExport(cl, ls(.GlobalEnv))
            junk <- parallel::clusterEvalQ(cl, expr = {library(mmBFPA)})
            
            # set seeds
            paraellel::clusterSetRNGStream(cl, 123)
            
            # run sampler
            chain_out <- parallel::parLapply(
                cl,
                X = seq_len(chains),
                fun = chain_call,
                pre_sparse = pre_sparse,
                warmup = warmup,
                iter = iter,
                thin = thin,
                dat = dat,
                mode = mode,
                share_noise = share_noise,
                d_mask = d_mask,
                n_levels = n_levels,
                margin_vals = margin_vals,
                x = x,
                lambda = lambda,
                zeros = zeros,
                omega = omega,
                alpha = alpha,
                gamma_k = gamma_k,
                d = d,
                dc = dc,
                ibp_a = ibp_a,
                ibp_b = ibp_b,
                benchmark = benchmark
            )
            stopCluster(cl)
        } else { 
            # set RNG
            RNGkind("L'Ecuyer-CMRG")
            
            # otherwise FORK
            chain_out <- parallel::mclapply(
                1:chains,
                FUN = function(X) {
                    chain_call(
                        pre_sparse = pre_sparse,
                        warmup = warmup,
                        iter = iter,
                        thin = thin,
                        dat = dat,
                        mode = mode,
                        share_noise = share_noise,
                        d_mask = d_mask,
                        n_levels = n_levels,
                        margin_vals = margin_vals,
                        x = x,
                        lambda = lambda,
                        zeros = zeros,
                        omega = omega,
                        alpha = alpha,
                        gamma_k = gamma_k,
                        d = d,
                        dc = dc,
                        ibp_a = ibp_a,
                        ibp_b = ibp_b,
                        benchmark = benchmark
                    )
                },
                mc.preschedule = FALSE,
                mc.set.seed = TRUE,
                mc.cores = n_cores
            )
        }
    }

    class(chain_out) <- "bpfa.results"
    attr(chain_out, "mode") <- mode
    attr(chain_out, "pre_sparse") <- pre_sparse
    attr(chain_out, "warmup") <- warmup
    attr(chain_out, "iter") <- iter
    attr(chain_out, "chains") <- chains
    attr(chain_out, "thin") <- thin
    attr(chain_out, "parallel") <- parallel
    attr(chain_out, "n_cores") <- n_cores
    attr(chain_out, "benchmark") <- benchmark
    return(chain_out)
}


#' Meta-function that bundles the sampler into a chain call
#' 
#' Performs single chain run of the mmBPFA algorithm.
#' 
#' @keywords internal
#' 
chain_call <- function(
    pre_sparse,
    warmup,
    iter,
    thin,
    dat,
    mode,
    share_noise,
    d_mask,
    n_levels,
    margin_vals,
    x,
    lambda,
    zeros,
    omega,
    alpha,
    gamma_k,
    d,
    dc,
    ibp_a,
    ibp_b,
    benchmark
    ) {
    # ================ #
    # Pre-sparse steps #
    # ================ #
    # initialize progress bar
    print(paste("Beginning pre-Beta process phase. Will run for", pre_sparse, "steps."))
    pb <- utils::txtProgressBar(min = 0, max = pre_sparse, style = 3)

    # initialize benchmark container
    if (benchmark) pre_time_log <- tibble::tibble("time" = list(NA), .rows = pre_sparse)

    for (step in 1:pre_sparse) {
        # get new samples
        c(x, zeros, lambda, omega, alpha, gamma_k, d, dc, ibp_a, ibp_b, time_log) %<-% sampler_phase(
            phase = "pre-sparse",
            dat = dat,
            mode = mode,
            share_noise = share_noise,
            d_mask = d_mask,
            n_levels = n_levels,
            margin_vals = margin_vals,
            x = x,
            lambda = lambda,
            alpha = alpha,
            omega = omega,
            zeros = zeros,
            gamma_k = gamma_k,
            d = d,
            dc = dc,
            ibp_a = ibp_a,
            ibp_b = ibp_b,
            benchmark = benchmark
        )
        # update progress bar
        utils::setTxtProgressBar(pb, step)

        # update log
        if (benchmark) pre_time_log$time[[step]] <- time_log
    }
    close(pb)

    # ============= #
    # warm-up steps #
    # ============= #
    # initialize new progress bar
    print(paste("Beginning warmup phase. Will run for", warmup, "steps."))
    pb <- utils::txtProgressBar(min = 0, max = warmup, style = 3)

    # initialize benchmark container
    if (benchmark) warmup_time_log <- tibble::tibble("time" = list(NA), .rows = warmup)

    for (step in 1:warmup) {
        # get new samples
        c(x, zeros, lambda, omega, alpha, gamma_k, d, dc, ibp_a, ibp_b, time_log) %<-% sampler_phase(
            phase = "warmup",
            dat = dat,
            mode = mode,
            share_noise = share_noise,
            d_mask = d_mask,
            n_levels = n_levels,
            margin_vals = margin_vals,
            x = x,
            lambda = lambda,
            alpha = alpha,
            omega = omega,
            zeros = zeros,
            gamma_k = gamma_k,
            d = d,
            dc = dc,
            ibp_a = ibp_a,
            ibp_b = ibp_b,
            benchmark = benchmark
            )
        # update progress bar
        utils::setTxtProgressBar(pb, step)

        # update log
        if (benchmark) warmup_time_log$time[[step]] <- time_log
    }
    close(pb)

    # ========== #
    # mcmc steps #
    # ========== #
    # initialize results containers
    mcmc_store <- tibble::tibble(
        "zeros" = list(NA),
        "lambda" = list(NA),
        "omega" = list(NA),
        "alpha" = list(NA),
        "gamma_k" = list(NA),
        "cutpoints" = list(NA),
        .rows = iter / thin
    )
    num_dims <- vector("numeric", length = iter / thin)
    log_liks <- vector("numeric", length = iter / thin)

    # initialize new progress bar
    print(paste("Beginning main sampling phase. Will run for", iter, "steps."))
    pb <- utils::txtProgressBar(min = 0, max = iter, style = 3)
    print_k <- floor(iter / 70)

    # initialize benchmark container
    if (benchmark)  mcmc_time_log <- tibble::tibble("time" = list(NA), .rows = iter)

    for (step in 1:iter) {
        # get new samples
        c(x, zeros, lambda, omega, alpha, gamma_k, d, dc, ibp_a,
        ibp_b, k_zeros, k_lambda, k_omega, k_gamma_k, k_dc,
        log_lik, n_active_dims, cutpoints, time_log
        ) %<-% sampler_phase(
            phase = "main",
            dat = dat,
            mode = mode,
            share_noise = share_noise,
            d_mask = d_mask,
            n_levels = n_levels,
            margin_vals = margin_vals,
            x = x,
            lambda = lambda,
            alpha = alpha,
            omega = omega,
            zeros = zeros,
            gamma_k = gamma_k,
            d = d,
            dc = dc,
            ibp_a = ibp_a,
            ibp_b = ibp_b,
            benchmark = benchmark
            )

        # store main phase results
        if (step %% thin == 0) {
            mcmc_store$zeros[[step / thin]] <- k_zeros
            mcmc_store$lambda[[step / thin]] <- k_lambda
            mcmc_store$omega[[step / thin]] <- k_omega
            mcmc_store$alpha[[step / thin]] <- alpha
            mcmc_store$gamma_k[[step / thin]] <- k_gamma_k
            mcmc_store$cutpoints[[step / thin]] <- cutpoints

            num_dims[step / thin] <- n_active_dims
            log_liks[step / thin] <- log_lik
        }

        # update progress bar
        if (step %% print_k == 0) cat("\rNumber of active factors: ", n_active_dims)
        utils::setTxtProgressBar(pb, step)

        # update log
        if (benchmark) mcmc_time_log$time[[step]] <- time_log
    }
    close(pb)

    # benchmark outputs
    if (benchmark) {
        benchmarks <- list(
            "pre_sparse" = pre_time_log,
            "warmup" = warmup_time_log,
            "mcmc" = mcmc_time_log
        )
    }

    # Chain outputs
    output_list <- list()
    output_list[["Number of Dimensions"]] <- num_dims
    output_list[["Log-Likelihoods"]] <- log_liks
    output_list[["MCMC Results"]] <- mcmc_store
    if (benchmark) output_list[["Benchmarks"]] <- benchmarks

    return(output_list)
}


#' All-in-one function for the sampling procedure 
#' 
#' @keywords internal
#' 
sampler_phase  <- function(
    phase,
    dat,
    mode,
    share_noise,
    d_mask,
    n_levels,
    margin_vals,
    x,
    lambda,
    alpha,
    omega,
    zeros,
    gamma_k,
    d,
    dc,
    ibp_a,
    ibp_b,
    benchmark
    ) {
    # initialize sparsity parameters for each phase
    if (phase == "pre-sparse") {
        sparse <- infinite <- FALSE
    } else {
        sparse <- infinite <- TRUE
    }

    # 1. sample x
    if (benchmark) tictoc::tic("sample x")
    c(x, cutpoints) %<-% sample_x(
        dat = dat,
        mode = mode,
        d_mask = d_mask,
        n_levels = n_levels,
        margin_vals = margin_vals,
        x = x,
        lambda = lambda,
        alpha = alpha,
        omega = omega
    )
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # 2. sample lambdas and zeros without beta-process prior
    if (benchmark) tictoc::tic("samples zeros and lambda")
    c(zeros, lambda, dc) %<-% sample_lambda_and_zeros(
        x = x,
        lambda = lambda,
        zeros = zeros,
        omega = omega,
        alpha = alpha,
        gamma_k = gamma_k,
        dc = dc,
        ibp_a = ibp_a,
        ibp_b = ibp_b,
        tau = 5,
        sparse = sparse,
        infinite = infinite
    )
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # Remove inactive features for warmup and main phases &
    # sample IBP parameters
    if (benchmark) tictoc::tic("remove inactive features and sample IBP")
    if (phase == "warmup" | phase == "main") {
        if (any(colSums(zeros) == 0)) {
            c(zeros, lambda, omega, gamma_k, dc) %<-% remove_inactive(
                zeros = zeros,
                lambda = lambda,
                omega = omega,
                gamma_k = gamma_k,
                dc = dc
                )
        }
        ibp_a <- sample_IBP_a(ibp_b = ibp_b, zeros = zeros)
        ibp_b <- sample_IBP_b(ibp_a = ibp_a, ibp_b = ibp_b, zeros = zeros)

    }
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # 3. sample omega
    if (benchmark) tictoc::tic("sample omega")
    omega <- sample_omega(
        lambda = lambda,
        x = x,
        alpha = alpha
    )
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # 4. sample alpha
    if (benchmark) tictoc::tic("sample alpha")
    alpha <- sample_alpha(
        x = x,
        lambda = lambda,
        omega = omega
    )
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # 5. sample gamma_k
    if (benchmark) tictoc::tic("sample gamma_k")
    gamma_k <- sample_gamma_k(
        zeros = zeros,
        lambda = lambda,
        d = d
    )
    if (benchmark) tictoc::toc(log = T, quiet = T)

    # 6. sample d
    if (share_noise) {
        if (benchmark) tictoc::tic("sample d")
        d <- sample_d(
            zeros = zeros,
            gamma_k = gamma_k,
            c = 1,
            c0 = 0,
            d0 = 1
        )
        if (benchmark) tictoc::toc(log = T, quiet = T)
    }

    # 7. remove weakly active features during main phase
    # calculate log-likelihood
    # calculate number of active dimensions
    if (phase == "main") {
        c(k_zeros, k_lambda, k_omega, k_gamma_k, k_dc) %<-% remove_weakly_active(
            zeros = zeros,
            lambda = lambda,
            omega = omega,
            gamma_k = gamma_k,
            dc = dc
        )
        log_lik <- calculate_log_lik(
            k_lambda = k_lambda,
            k_omega = k_omega,
            alpha = alpha,
            dat = dat
        )
        n_active_dims <- ncol(k_lambda)
    }

    # store values
    output_list <- list()
    output_list[["x"]] <- x
    output_list[["zeros"]] <- zeros
    output_list[["lambda"]] <- lambda
    output_list[["omega"]] <- omega
    output_list[["alpha"]] <- alpha
    output_list[["gamma_k"]] <- gamma_k
    output_list[["d"]] <- d
    output_list[["dc"]] <- dc
    output_list[["ibp_a"]] <- ibp_a
    output_list[["ibp_b"]] <- ibp_b

    if (phase == "main") {
        output_list[["k_zeros"]] <- k_zeros
        output_list[["k_lambda"]] <- k_lambda
        output_list[["k_omega"]] <- k_omega
        output_list[["k_gamma_k"]] <- k_gamma_k
        output_list[["k_dc"]] <- k_dc
        output_list[["log_lik"]] <- log_lik
        output_list[["n_active_dims"]] <- n_active_dims
        output_list[["cutpoints"]] <- cutpoints
    }

    # store timings
    if (benchmark) {
        time_log <- tictoc::tic.log(format = F)
        tictoc::tic.clearlog()
        output_list[["time_log"]] <- time_log
    } else {
        output_list[["time_log"]] <- NA
    }

    return(output_list)
}