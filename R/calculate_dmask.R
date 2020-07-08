# calculates the D subset (rank orderings mask), and number of levels and unique values per margin
#
# @param dat data matrix.
# @param mode Margin type, must be one of "fixed", "multi" or "mixed".
# 
# @return List containing the D mask, number of unique values per margin, and unique margin values
# 
calculate_d_mask <- function(dat, mode) {

    # calculate D mask, margin cardinality, and unique margin values
    # for fixed margins
    if (mode == "fixed") {
        margin_vals <- sort(unique(as.vector(dat)))
        n_levels <- length(margin_vals)
        d_mask <- lapply(margin_vals, function(x) {
            list(
                "lower" = tibble::as_tibble(dat < x, .name_repair = "minimal"),
                "upper" = tibble::as_tibble(dat > x, .name_repair = "minimal")
            )
        })
        d_mask <- list(
            purrr::map(d_mask, ~ .x[["lower"]]),
            purrr::map(d_mask, ~ .x[["upper"]])
            )
    } else { # for mixed margins
        dat_tbl <- tibble:::as_tibble(dat, .name_repair = "minimal")
        margin_vals <- purrr::map(dat_tbl, ~ sort(unique(.x)))
        n_levels <- purrr::map_int(margin_vals, length)
        d_mask <- purrr::map2(
            margin_vals,
            dat_tbl,
            ~ purrr::map(
                .x,
                function(x) {
                    list(
                        "lower" = .y < x,
                        "upper" = .y > x
                    )
                }
            )
        )
        d_mask <- list(
            purrr::map_depth(d_mask, 2, ~ .x[["lower"]]),
            purrr::map_depth(d_mask, 2, ~ .x[["upper"]])
        )
    }

    # store outputs
    output_list <- list()
    output_list[["d_mask"]] <- d_mask
    output_list[["n_levels"]] <- n_levels
    output_list[["margin_vals"]] <- margin_vals

    return(output_list)
}