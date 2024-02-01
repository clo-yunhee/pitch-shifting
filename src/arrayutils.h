#ifndef ARRAYUTILS_H__
#define ARRAYUTILS_H__

#include <fftw3.h>

#include <complex>

inline fftw_complex *cpx_stl2fftw(std::complex<double> *arr) {
    return reinterpret_cast<fftw_complex *>(arr);
}

inline std::complex<double> *cpx_fftw2stl(fftw_complex *arr) {
    return reinterpret_cast<std::complex<double> *>(arr);
}

void circshift(const double *in, int L, int shift, double *out);

void fftshift(const double *in, int L, double *out);

int positiverem(int a, int b);

void clear_array(double *in, int L);

void fold_array(const double *in, int Lin, int offset, int Lfold, double *out);

void periodize_array(const double *in, int Lin, int Lout, double *out);

#endif  // ARRAYUTILS_H__