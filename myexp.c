#include "omp.h"

#include "myexp.h"

#pragma omp declare simd
double myexp (double x) {
    return x;
}
