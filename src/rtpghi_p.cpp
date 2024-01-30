#include "rtpghi_p.h"

#include <algorithm>
#include <complex>
#include <cstdio>

#include "rtpghi.h"
#include "rtpghi_heap.h"

static double princarg(double in) {
    return (in - 2.0 * M_PI * std::round(in / (2.0 * M_PI)));
}

static void shiftcolsleft(double *cols, int height, int N,
                          const double *newcol) {
    // rotate left by 1*height is same as rotate left by 1 column
    if (newcol) {
        std::rotate(cols, cols + height, cols + (N - 1) * height);
        std::copy(newcol, newcol + height, cols + (N - 1) * height);
    } else {
        std::rotate(cols, cols + height, cols + N * height);
    }
}

static void rtpghi_abs(const std::complex<double> *in, int height,
                       double *out) {
    for (int ii = 0; ii < height; ++ii) {
        out[ii] = std::abs(in[ii]);
    }
}

static void rtpghi_phase(const std::complex<double> *in, int height,
                         double *out) {
    for (int ii = 0; ii < height; ++ii) {
        out[ii] = std::arg(in[ii]);
    }
}

/** Compute phase frequency gradient by differentiation in time */
static void rtpghi_fgrad(const double *phase, int M, double stretch,
                         double *fgrad);

/** Compute phase time gradient by differentiation in frequency */
static void rtpghi_tgrad(const double *phase, int aanaprev, int aananext, int M,
                         double stretch, double *tgrad);

static void rtpghi_magphase(const double *s, const double *phase, int L,
                            std::complex<double> *c);

rtpghi_priv::rtpghi_priv(int W, int a, int M, double tol) {
    int M2;

    rtpghi_assert(W > 0, "W must be positive");
    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(tol > 0 && tol < 1, "tol must be in range ]0,1[");

    M2 = M / 2 + 1;

    _p = new rtpghi_update_plan(M, W, tol);
    _s.resize(3 * M2 * W);
    _tgrad.resize(2 * M2 * W);
    _fgrad.resize(1 * M2 * W);
    _phase.resize(1 * M2 * W);
    _phasein.resize(3 * M2 * W);

    _M = M;
    _a = a;
    _W = W;
    _stretch = 1.0;
}

rtpghi_priv::~rtpghi_priv() { delete _p; }

double rtpghi_priv::get_stretch() const { return _stretch; }

void rtpghi_priv::set_tolerance(double tol) {
    rtpghi_assert(tol > 0 && tol < 1, "tol must be in range ]0,1[");
    _p->_tol = tol;
}

void rtpghi_priv::reset(const double **sinit) {
    int M2 = _M / 2 + 1;

    std::fill(_s.begin(), _s.end(), 0);
    std::fill(_tgrad.begin(), _tgrad.end(), 0);
    std::fill(_fgrad.begin(), _fgrad.end(), 0);
    std::fill(_phase.begin(), _phase.end(), 0);
    std::fill(_phasein.begin(), _phasein.end(), 0);

    if (sinit) {
        for (int w = 0; w < _W; ++w) {
            if (sinit[w]) {
                std::copy(sinit[w], sinit[w] + 2 * M2, _s.begin() + 2 * w * M2);
            }
        }
    }
}

void rtpghi_priv::execute(const std::complex<double> *cin, double stretch,
                          std::complex<double> *cout) {
    // n, n-1, n-2 frames
    // s is n-th
    int    M2, W, asyn, aanaprev, aananext;
    double stretchmid;

    M2 = _M / 2 + 1;
    W = _W;
    asyn = _a;
    aanaprev = std::round(asyn / _stretch);  // old stretch
    aananext = std::round(asyn / stretch);   // new stretch

    for (int w = 0; w < W; ++w) {
        double *sCol = _s.data() + 3 * w * M2;
        double *tgradCol = _tgrad.data() + 2 * w * M2;
        double *fgradCol = _fgrad.data() + 1 * w * M2;
        double *phaseCol = _phase.data() + 1 * w * M2;
        double *phaseinCol = _phasein.data() + 3 * w * M2;

        shiftcolsleft(sCol, M2, 3, nullptr);
        shiftcolsleft(tgradCol, M2, 2, nullptr);
        shiftcolsleft(phaseinCol, M2, 3, nullptr);

        rtpghi_abs(cin + w * M2, M2, sCol + 2 * M2);
        rtpghi_phase(cin + w * M2, M2, phaseinCol + 2 * M2);

        rtpghi_tgrad(phaseinCol, aanaprev, aananext, _M, _stretch,
                     tgradCol + M2);

        if (std::abs(stretch - 1.0) < 1e-4) {
            // Bypass if no stretching is done
            std::copy(phaseinCol + M2, phaseinCol + 2 * M2, phaseCol);
        } else {
            rtpghi_fgrad(phaseinCol + M2, _M, _stretch, fgradCol);
            _p->execute(sCol, tgradCol, fgradCol, phaseCol, phaseCol);
        }

        // Combine phase with amplitude
        rtpghi_magphase(sCol + M2, phaseCol, M2, cout + w * M2);
    }

    // Only update stretch for the next frame
    _stretch = stretch;
}

rtpghi_update_plan::rtpghi_update_plan(int M, int W, double tol)
    : _rand(_rd()), _dist(-2.0 * M_PI, 2.0 * M_PI) {
    int M2 = M / 2 + 1;

    _donemask.resize(M2);

    _tol = tol;
    _M = M;
    _h = std::make_unique<rtpghi_heap_t>(2 * M2, nullptr);
}

void rtpghi_update_plan::execute_with_mask(const double *s, const double *tgrad,
                                           const double *fgrad,
                                           const double *startphase,
                                           const int8_t *mask, double *phase) {
    int M2 = _M / 2 + 1;
    std::copy(mask, mask + M2, _donemask.begin());

    return execute_common(s, tgrad, fgrad, startphase, phase);
}

// s: M2 x 2
// tgrad: M2 x 2
// fgrad: M2 x 1
// startphase: M2 x 1
// phase: M2 x 1
// donemask: M2 x 1
// heap must be able to hold 2 * M2 values
void rtpghi_update_plan::execute(const double *s, const double *tgrad,
                                 const double *fgrad, const double *startphase,
                                 double *phase) {
    std::fill(_donemask.begin(), _donemask.end(), 0);

    return execute_common(s, tgrad, fgrad, startphase, phase);
}

void rtpghi_update_plan::execute_common(const double *s, const double *tgrad,
                                        const double *fgrad,
                                        const double *startphase,
                                        double       *phase) {
    int M2 = _M / 2 + 1;
    int quickbreak = M2;
    // We only need to compute M2 values, so perform quick exit
    // if we have them, but the heap is not yet empty.
    // (deleting from heap involves many operations)
    const double *slog2 = s + M2;

    // Find max and the absolute threshold
    double logabstol = s[0];
    for (int m = 1; m < 2 * M2; ++m) {
        if (s[m] > logabstol) {
            logabstol = s[m];
        }
    }

    logabstol *= _tol;

    _h->reset(s);

    for (int m = 0; m < M2; ++m) {
        if (_donemask[m] > 0) {
            // We already know this one
            _h->push(m + M2);
            quickbreak -= 1;
        } else {
            if (slog2[m] <= logabstol) {
                // This will get a randomly generated phase
                _donemask[m] = -1;
                quickbreak -= 1;
            } else {
                _h->push(m);
            }
        }
    }

    int w = -1;
    while ((quickbreak > 0) && (w = _h->pop()) >= 0) {
        if (w >= M2) {
            // Next frame
            int wprev = w - M2;

            if (wprev != M2 - 1 && !_donemask[wprev + 1]) {
                phase[wprev + 1] =
                    phase[wprev] + (fgrad[wprev] + fgrad[wprev + 1]) / 2.0;
                _donemask[wprev + 1] = 1;

                _h->push(w + 1);
                quickbreak -= 1;
            }

            if (wprev != 0 && !_donemask[wprev - 1]) {
                phase[wprev - 1] =
                    phase[wprev] - (fgrad[wprev] + fgrad[wprev - 1]) / 2.0;
                _donemask[wprev - 1] = 1;

                _h->push(w - 1);
                quickbreak -= 1;
            }
        } else {
            // Current frame
            if (!_donemask[w]) {
                int wnext = w + M2;
                phase[w] = startphase[w] + (tgrad[w] + tgrad[wnext]) / 2.0;
                _donemask[w] = 1;

                _h->push(wnext);
                quickbreak -= 1;
            }
        }
    }

    // Fill in values below tol
    for (int ii = 0; ii < M2; ++ii) {
        if (_donemask[ii] < 0) {
            phase[ii] = random_phase();
        }
    }
}

double rtpghi_update_plan::random_phase() { return _dist(_rand); }

static void rtpghi_fgrad(const double *phase, int M, double stretch,
                         double *fgrad) {
    int M2 = M / 2 + 1;

    for (int m = 1; m < M2 - 1; ++m) {
        fgrad[m] = (princarg(phase[m + 1] - phase[m]) +
                    princarg(phase[m] - phase[m - 1])) /
                   (2.0) * stretch;
    }

    fgrad[0] = phase[0] * stretch;
    fgrad[M2 - 1] = phase[M2 - 1] * stretch;
}

static void rtpghi_tgrad(const double *phase, int aanaprev, int aananext, int M,
                         double stretch, double *tgrad) {
    int M2 = M / 2 + 1;
    // a is asyn
    double asyn = aanaprev * stretch;

    double const1prev = 2.0 * M_PI * ((double)aanaprev) / M;
    double const1next = 2.0 * M_PI * ((double)aananext) / M;
    double const2 = 2.0 * M_PI * asyn / M;

    const double *pcol0 = phase;
    const double *pcol1 = phase + 1 * M2;
    const double *pcol2 = phase + 2 * M2;

    for (int m = 0; m < M2; ++m) {
        tgrad[m] = asyn * (princarg(pcol2[m] - pcol1[m] - const1next * m) /
                               (2.0 * aananext) -
                           princarg(pcol1[m] - pcol0[m] - const1prev * m) /
                               (2.0 * aanaprev)) +
                   const2 * m;
    }
}

static void rtpghi_magphase(const double *s, const double *phase, int L,
                            std::complex<double> *c) {
    for (int l = 0; l < L; ++l) {
        c[l] = std::polar(s[l], phase[l]);
    }
}

#ifndef NDEBUG
void __rtpghi_assert(const char *expr_str, bool expr, const char *file,
                     int line, const char *msg) {
    if (!expr) {
        fprintf(stderr,
                "Assert failed:\t%s\nExpected:\t%s\nSource:\t\t%s, line %d\n",
                msg, expr_str, file, line);
    }
}
#endif