// Number of OpenMP threads used inside a chain. Every parallel unit (column,
// row, item) draws from its own RNG stream and reductions run in a fixed
// order, so results do not depend on the thread count.
#ifndef MOSAICFA_THREADS_H
#define MOSAICFA_THREADS_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mosaicfa {
inline int& num_threads() { static int n = 1; return n; }
}

#endif
