#include <Rcpp.h>
#include "threads.h"

//' Whether the package was built with OpenMP
//' @keywords internal
// [[Rcpp::export(name = ".openmp_available")]]
bool openmp_available() {
#ifdef _OPENMP
    return true;
#else
    return false;
#endif
}

//' Set the number of threads used inside a chain
//' @keywords internal
// [[Rcpp::export(name = ".set_threads")]]
void set_threads(int n) { mosaicfa::num_threads() = n < 1 ? 1 : n; }
