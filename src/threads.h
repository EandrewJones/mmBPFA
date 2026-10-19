// Number of OpenMP threads used inside a chain. Every parallel unit (column,
// row, item) draws from its own RNG stream and reductions run in a fixed
// order, so results do not depend on the thread count.
#ifndef MMBPFA_THREADS_H
#define MMBPFA_THREADS_H

#ifdef _OPENMP
#include <omp.h>
#endif

namespace mmbpfa {
inline int& num_threads() { static int n = 1; return n; }
}

#endif
