#include "rtdgtreal_p.h"

#include <fftw3.h>

#include "arrayutils.h"
#include "rtpghi.h"

rtdgtreal_priv::rtdgtreal_priv(const double *g, int gl, int M,
                               const rtdgt_phase_t            ptype,
                               const dgt_transformdirection_t tradir) {
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