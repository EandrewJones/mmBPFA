// C++ view of the compressed observed-cell layout built by prepare_data() in R.
#ifndef MOSAICFA_DATA_STRUCTURE_H
#define MOSAICFA_DATA_STRUCTURE_H

#include <RcppArmadillo.h>

namespace mosaicfa {

struct Data {
    int n, p;
    Rcpp::IntegerVector row_idx, col_idx, col_ptr, n_levels, level_ptr, lvl_col_ptr, row_ptr, cell_of_row;

    explicit Data(const Rcpp::List& d)
        : n(Rcpp::as<int>(d["n"])), p(Rcpp::as<int>(d["p"])),
          row_idx(d["row_idx"]), col_idx(d["col_idx"]), col_ptr(d["col_ptr"]), n_levels(d["n_levels"]),
          level_ptr(d["level_ptr"]), lvl_col_ptr(d["lvl_col_ptr"]),
          row_ptr(d["row_ptr"]), cell_of_row(d["cell_of_row"]) {}

    // cells [level_begin(j, v), level_begin(j, v + 1)) hold level v of column j
    int level_begin(int j, int v) const { return level_ptr[lvl_col_ptr[j] + v]; }
};

} // namespace mosaicfa

#endif
