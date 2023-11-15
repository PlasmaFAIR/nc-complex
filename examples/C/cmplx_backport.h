// Clang doesn't have full support for complex.h
#ifndef CMPLX
#define CMPLX(x, y) ((double)(x) + (double)(y) * (double_complex)I)
#endif
