#' Resolve prior settings
#'
#' Turns a preset name, or a list of raw hyperparameters (optionally naming a
#' preset), into the complete set of priors the sampler uses. Raw values
#' override the preset.
#'
#' Presets:
#' \describe{
#'   \item{uniform}{slab InvGamma(1, 1/2): for an item on one active factor, the copula
#'     loading l / sqrt(1 + l^2) is Uniform(-1, 1).}
#'   \item{conservative}{slab InvGamma(1, 1/4): smaller loadings a priori.}
#' }
#' Both use spike variance 0.05 (the CUSP default), a Gamma(1, 1) prior on the
#' CUSP concentration, and Beta(1, 1) item inclusion within columns.
#'
#' @param priors a preset name, or a list of hyperparameters (\code{slab_shape},
#'   \code{slab_rate}, \code{spike_var}, \code{a_cusp_shape}, \code{a_cusp_rate},
#'   \code{rho_a}, \code{rho_b}) with an optional \code{preset}.
#' @param p number of items.
#' @param H truncation (maximum number of columns). Default min(floor((p - 1) / 2), 50).
#'
#' @return Named list of every prior setting, including \code{preset} and \code{H}.
#'
#' @export
resolve_priors <- function(priors = "uniform", p, H = NULL) {
    presets <- list(
        uniform = list(slab_shape = 1, slab_rate = 0.5),
        conservative = list(slab_shape = 1, slab_rate = 0.25)
    )
    shared <- list(spike_var = 0.05, a_cusp_shape = 1, a_cusp_rate = 1, rho_a = 1, rho_b = 1,
                   cut_sd = 1e300)   # SPIKE: cutpoint prior sd (explicit-cutpoint items only)

    raw <- if (is.character(priors)) list(preset = priors) else as.list(priors)
    preset <- if (is.null(raw$preset)) "uniform" else raw$preset
    if (!preset %in% names(presets)) stop("Unknown prior preset: ", preset)
    out <- c(presets[[preset]], shared)
    unknown <- setdiff(names(raw), c(names(out), "preset"))
    if (length(unknown)) stop("Unknown prior field: ", paste(unknown, collapse = ", "))
    out[names(raw)[names(raw) != "preset"]] <- raw[names(raw) != "preset"]

    if (is.null(H)) H <- min(floor((p - 1) / 2), 50)
    if (H < 1) stop("H must allow at least one column; got H = ", H, " for p = ", p)
    c(list(preset = preset), out, list(H = as.integer(H)))
}


#' Retrieve every prior setting used for a fit
#'
#' @param x sampler output (class \code{mosaicfa.results}) or processed results.
#'
#' @return Named list of prior settings.
#'
#' @export
get_priors <- function(x) {
    if (inherits(x, "mcmc.output.processed")) return(x[["Priors"]])
    attr(x, "priors")
}
