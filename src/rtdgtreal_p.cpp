#include "rtdgtreal_p.h"

#include <fftw3.h>

#include "rtpghi.h"

#define cpx_stl2fftw(arr) reinterpret_cast<fftw_complex *>(arr)
#define cpx_fftw2stl(arr) reinterpret_cast<std::complex<double> *>(arr)

static void circshift(const double *in, int L, int shift, double *out) {
    int p;
    rtpghi_assert(L > 0, "L must be positive");

    // Fix shift
    p = (L - shift) % L;

    if (p < 0) p += L;

    if (in == out) {
        if (p) {  // Do nothing if no shift is needed
            int m, count, i, j;

            for (m = 0, count = 0; count != L; ++m) {
                double t = in[m];

                for (i = m, j = m + p; j != m;
                     i = j, j = j + p < L ? j + p : j + p - L, ++count)
                    out[i] = out[j];

                out[i] = t;
                ++count;
            }
        }
    } else {
        // Still ok if p == 0
        std::copy(in + p, in + p + (L - p), out);
        std::copy(in, in + p, out + L - p);
    }
}

static void fftshift(const double *in, int L, double *out) {
    circshift(in, L, (L / 2), out);
}

static int positiverem(int a, int b) {
    int c = a % b;
    return (c < 0 ? c + b : c);
}

static void clear_array(double *in, int L) {
    rtpghi_assert(L >= 0, "L must be nonnegative");
    if (L > 0) std::fill(in, in + L, 0);
}

static void fold_array(const double *in, int Lin, int offset, int Lfold,
                       double *out) {
    int startIdx;

    rtpghi_assert(Lin > 0, "Lin must be positive");
    rtpghi_assert(Lfold > 0, "Lfold must be positive");

    // Sanitize offset.
    startIdx = positiverem(offset, Lfold);

    // Clear output, we will use it as an accumulator
    if (in != out)
        clear_array(out, Lfold);
    else if (Lfold > Lin)
        clear_array(out + Lin, Lfold - Lin);

    if (!startIdx) {
        // Common code for no offset
        int startAt = (in == out) ? Lfold : 0;

        for (int ii = startAt; ii < Lin;)
            for (int kk = 0; ii < Lin && kk < Lfold; ++ii, ++kk)
                out[kk] += in[ii];
    } else {
        if (in == out) {
            // We cannot avoid the (slow) inplace circshift anyway.
            // Do it after the folding.
            fold_array(in, Lin, 0, Lfold, out);
            circshift(in, Lfold, startIdx, out);
        } else {
            // We avoid the inplace circshift by effectively
            // doing circshift of all blocks.
            for (int ii = 0; ii < Lin;) {
                int kk = startIdx;
                for (; kk < Lfold && ii < Lin; ++ii, ++kk) out[kk] += in[ii];

                for (kk = 0; kk < startIdx && ii < Lin; ++ii, ++kk)
                    out[kk] += in[ii];
            }
        }
    }
}

static void periodize_array(const double *in, int Lin, int Lout, double *out) {
    rtpghi_assert(Lin > 0, "Lin must be positive");
    rtpghi_assert(Lout > 0, "Lout must be positive");

    // Do nothing if there is no place where to put periodized samples
    if (Lout <= Lin) {
        if (in != out) {
            std::copy(in, in + Lout, out);
        }
    } else {
        int periods = Lout / Lin;
        int lastL = Lout - periods * Lin;
        int startPer = (in == out) ? 1 : 0;

        for (int ii = startPer; ii < periods; ++ii)
            std::copy(in, in + Lin, out + ii * Lin);

        std::copy(in, in + lastL, out + periods * Lin);
    }
}

rtdgtreal_priv::rtdgtreal_priv(const double *g, int gl, int M,
                               const rtdgt_phase_t          ptype,
                               const dgt_transformdirection tradir) {
    int M2;

    rtpghi_assert(gl > 0, "gl must be positive");
    rtpghi_assert(M > 0, "M must be positive");

    M2 = M / 2 + 1;
    _fftBufLen = gl > 2 * M2 ? gl : 2 * M2;

    _g = fftw_alloc_real(gl);
    _fftBuf = fftw_alloc_real(_fftBufLen);
    _fftBuf_cpx = cpx_fftw2stl(fftw_alloc_complex(M2));
    _gl = gl;
    _M = M;
    _ptype = ptype;

    fftshift(g, gl, _g);

    if (tradir == DGT_FORWARD) {
        _pfft = fftw_plan_dft_r2c_1d(M, _fftBuf, cpx_stl2fftw(_fftBuf_cpx),
                                     FFTW_MEASURE);
    } else if (tradir == DGT_INVERSE) {
        _pfft = fftw_plan_dft_c2r_1d(M, cpx_stl2fftw(_fftBuf_cpx), _fftBuf,
                                     FFTW_MEASURE);
    }
}

rtdgtreal_priv::~rtdgtreal_priv() {
    fftw_destroy_plan(_pfft);
    fftw_free(_g);
    fftw_free(_fftBuf);
    fftw_free(cpx_stl2fftw(_fftBuf_cpx));
}

void rtdgtreal_priv::execute_fwd(const double *f, int W,
                                 std::complex<double> *c) {
    int M2;

    rtpghi_assert(W > 0, "W must be positive");

    M2 = _M / 2 + 1;

    for (int w = 0; w < W; ++w) {
        const double         *fchan = f + w * _gl;
        std::complex<double> *cchan = c + w * M2;

        if (_g) {
            for (int ii = 0; ii < _gl; ++ii) {
                _fftBuf[ii] = fchan[ii] * _g[ii];
            }
        }

        if (_M > _gl) {
            std::fill(_fftBuf + _gl, _fftBuf + _gl + (_M - _gl), 0);
        }

        if (_gl > _M) {
            fold_array(_fftBuf, _gl, _M, 0, _fftBuf);
        }

        if (_ptype == RTDGTPHASE_ZERO) {
            circshift(_fftBuf, _M, -(_gl / 2), _fftBuf);
        }

        fftw_execute(_pfft);

        std::copy(_fftBuf_cpx, _fftBuf_cpx + M2, cchan);
    }
}

void rtdgtreal_priv::execute_inv(const std::complex<double> *c, int W,
                                 double *f) {
    int M2;

    rtpghi_assert(W > 0, "W must be positive");

    M2 = _M / 2 + 1;

    for (int w = 0; w < W; ++w) {
        const std::complex<double> *cchan = c + w * M2;
        double                     *fchan = f + w * _gl;

        std::copy(cchan, cchan + M2, _fftBuf_cpx);

        fftw_execute(_pfft);

        if (_ptype == RTDGTPHASE_ZERO) {
            circshift(_fftBuf, _M, _gl / 2, _fftBuf);
        }

        if (_gl > _M) {
            periodize_array(_fftBuf, _M, _gl, _fftBuf);
        }

        if (_g) {
            for (int ii = 0; ii < _gl; ++ii) {
                _fftBuf[ii] *= _g[ii];
            }
        }

        std::copy(_fftBuf, _fftBuf + _gl, fchan);
    }
}