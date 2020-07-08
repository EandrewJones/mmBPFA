# not in function
`%nin%`  <- Negate(`%in%`)


#' Indian Buffet Process log-beta draws
#' @keywords internal
ibp_lbeta <- function(m_k, p, b) {
    return(lbeta(m_k, p - m_k + b))
}


#' Indian Buffet point process parameter function
#' @keywords internal
Hpb <- function(b, p) {
    return(sum(b / (b + seq(1, p) - 1)))
}
Hpb  <- Vectorize(Hpb, "b")


#' function to get extract dimensions of 2d objects
#' @keywords internal
get_dims_2d <- function(x, warn = FALSE) {
    d1 <- dim(x)[1]
    d2 <- dim(x)[2]
    if (warn == TRUE) {
        if (d1 < 100 | d2 < 100) {
            warning("The number of observations and/or the number of covariates is lower than 100.  mmBPFA performs best when both are relatively large.  Be aware that your final results may have a dimensionality much lower than the truth.")
        }
    }
    return(list(d1, d2))
}


#' Selects appropriate K_start based on rank of data
#' @keywords internal
check_k_start <- function(K_start, p, n) {
    n_ks <- min(K_start, p, n)
    if (n_ks != K_start) {
        warning(paste0("K_start has been reduced from ", K_start, " to ", n_ks, " to keep K < p,n"))
        K_start <- n_ks
    }
    return(K_start)
}


#' vectorized which.max function for D mask
#' @keywords internal
vectorized_which.max <- Vectorize(
            function(x, y) {
                which.max(x < y) - 1
            },
            "x"
        )
