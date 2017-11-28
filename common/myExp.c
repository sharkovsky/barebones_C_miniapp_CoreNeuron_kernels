#include "omp.h"

#include "common/myExp.h"

#pragma omp declare simd
double myExp (double x) {
    return x;
}
