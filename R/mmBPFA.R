#' @keywords internal
#' @import rlang
#' @import abind
#' @import dplyr
#' @import graphics
#' @import mvtnorm
#' @import parallel
#' @import tictoc
#' @import tibble
#' @import tidyr
#' @import truncnorm
#' @import zeallot
#' @importFrom purrr accumulate discard every keep map map2 map2_chr
#'   map2_dbl map2_df map2_int map2_lgl map_at map_chr map_dbl map_df
#'   map_if map_int map_lgl pmap pmap_chr pmap_dbl pmap_df pmap_int
#'   pmap_lgl reduce some transpose
#' @importFrom Rcpp evalCpp
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @useDynLib mmBPFA, .registration = TRUE
"_PACKAGE"

globalVariables(".")