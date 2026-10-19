#' Compressed layout of the observed cells
#'
#' Stores only observed cells, column by column (CSC), sorted by value within
#' each column so the levels of the extended rank likelihood are contiguous.
#' A row-wise index (CSR) into the same cells supports row updates. Sampler
#' kernels work on this layout so their cost scales with the number of
#' observed cells rather than N x P.
#'
#' @keywords internal
#'
#' @param dat numeric data matrix; NA marks a missing cell.
#'
#' @return List with
#'   \code{n}, \code{p};
#'   \code{row_idx}, \code{col_idx}, \code{value}, \code{level}: per cell (0-based row and column, observed value, 0-based level within its column);
#'   \code{col_ptr}: length p + 1, cells of column j are \code{col_ptr[j] + 1} to \code{col_ptr[j + 1]};
#'   \code{n_levels}: distinct observed values per column;
#'   \code{level_ptr}, \code{lvl_col_ptr}: level v of column j spans cells
#'     \code{level_ptr[lvl_col_ptr[j] + v]} up to (not including) \code{level_ptr[lvl_col_ptr[j] + v + 1]};
#'   \code{row_ptr}, \code{cell_of_row}: the cells of row i, as 0-based cell indices.
prepare_data <- function(dat) {
    dat <- as.matrix(dat)
    if (!is.numeric(dat)) stop("All elements in the data must be numeric or integers")
    n <- nrow(dat); p <- ncol(dat)
    obs <- which(!is.na(dat))
    row <- as.integer((obs - 1L) %% n)
    col <- as.integer((obs - 1L) %/% n)
    val <- as.numeric(dat[obs])
    ord <- order(col, val)
    row <- row[ord]; col <- col[ord]; val <- val[ord]
    nnz <- length(val)

    new_level <- if (nnz) c(TRUE, col[-1] != col[-nnz] | val[-1] != val[-nnz]) else logical(0)
    n_levels <- tabulate(col[new_level] + 1L, p)
    lvl_col_ptr <- c(0L, cumsum(n_levels))
    level <- as.integer(cumsum(new_level) - 1L - lvl_col_ptr[col + 1L])

    list(
        n = n, p = p,
        row_idx = row, col_idx = col, value = val, level = level,
        col_ptr = c(0L, cumsum(tabulate(col + 1L, p))),
        n_levels = n_levels,
        level_ptr = c(which(new_level) - 1L, nnz),
        lvl_col_ptr = lvl_col_ptr,
        row_ptr = c(0L, cumsum(tabulate(row + 1L, n))),
        cell_of_row = order(row, col) - 1L
    )
}


#' Scatter per-cell values back into an N x P matrix
#'
#' @keywords internal
#'
#' @param data output of \code{prepare_data}.
#' @param values numeric vector with one entry per observed cell, in the layout's order.
#'
#' @return N x P matrix with NA in missing cells.
cells_to_matrix <- function(data, values) {
    out <- matrix(NA_real_, data$n, data$p)
    col <- rep.int(seq_len(data$p), diff(data$col_ptr))
    out[cbind(data$row_idx + 1L, col)] <- values
    out
}


#' Rank-consistent starting values for the latent data
#'
#' Each level of a column gets the normal score of its mid-rank, so values
#' increase strictly across levels. The x-step requires a starting point that
#' satisfies the ordering: from an inconsistent one, crossed bounds collapse
#' draws onto each other.
#'
#' @keywords internal
#'
#' @param data output of \code{prepare_data}.
#'
#' @return Numeric vector with one value per observed cell.
initial_x <- function(data) {
    nnz <- length(data$value)
    if (nnz == 0) return(numeric(0))
    col <- rep.int(seq_len(data$p), diff(data$col_ptr))
    n_j <- diff(data$col_ptr)[col]
    pos <- seq_len(nnz) - data$col_ptr[col]               # 1-based position within column
    level_id <- data$lvl_col_ptr[col] + data$level         # global level index
    mid <- ave(pos, level_id, FUN = mean)
    qnorm(mid / (n_j + 1))
}
