#include "gabdual_painless.h"

#include "arrayutils.h"
#include "rtpghi.h"

#define GABDIAGAPPLY(gg)                                                      \
    do {                                                                      \
        for (ltfat_int ii = 0; ii < domod.quot + domod.rem; ii++)             \
            (gg)[ii] = g[ii] * d[ii % a];                                     \
                                                                              \
        for (ltfat_int ii = gl - 1, jj = a - 1; ii >= domod.quot + domod.rem; \
             ii--, jj--) {                                                    \
            if (jj < 0) jj = a - 1;                                           \
            (gg)[ii] = g[ii] * d[jj];                                         \
        }                                                                     \
    } while (false)

/*
#define CHECKGABPAINLESS do{ \
    CHECK(LTFATERR_BADSIZE, gl > 0, "gl must be positive");\
    CHECK(LTFATERR_NOTPOSARG, a > 0, "a must be positive");\
    CHECK(LTFATERR_NOTPOSARG, M > 0, "M must be positive");\
    CHECK(LTFATERR_NOTAFRAME, M > a && gl >= a, "Not a frame. Check if M > a &&
gl >= a (passed M=%td, a=%td)",M,a);\
    CHECK(LTFATERR_NOTPAINLESS, M >= gl, "Not painless. Check if M>=gl");\
}while(0)
*/

void gabframediag(const double *g, int gl, int a, int M, int dl, double *d) {
    std::div_t domod;
    int        amax;

    rtpghi_assert(gl > 0, "gl must be positive");
    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(dl > 0, "dl must be positive");

    amax = a > dl ? dl : a;
    // d is used as an accumulator
    std::fill(d, d + dl, 0);

    domod = std::div(gl, 2);

    // First half
    for (int aIdx = 0; aIdx < amax; ++aIdx) {
        for (int ii = aIdx; ii < domod.quot + domod.rem; ii += a) {
            double gabs = g[ii];
            d[aIdx] += gabs * gabs;
        }
    }

    // Second half from the back
    for (int aIdx = amax - 1; aIdx >= 0; --aIdx) {
        for (int ii = gl - (a - aIdx); ii >= domod.quot + domod.rem; ii -= a) {
            double gabs = g[ii];
            d[aIdx] += gabs * gabs;
        }
    }

    for (int aIdx = 0; aIdx < amax; ++aIdx) {
        d[aIdx] *= M;
    }

    // frame diagonal is a-periodic
    if (dl > a) {
        periodize_array(d, a, dl, d);
    }
}

void gabdual_painless(const double g[], int gl, int a, int M, double gd[]) {
    std::div_t domod;

    rtpghi_assert(gl > 0, "gl must be positive");
    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(M > a && gl >= a, "Not a frame. Check if M > a && gl >= a");
    rtpghi_assert(M >= gl, "Not painless. Check if M >= gl");

    // This temporary array cannot be substituted by gd.
    auto d = std::make_unique<double[]>(a);

    gabframediag(g, gl, a, M, a, d.get());

    domod = std::div(gl, 2);

    // Invert the diagonal
    for (int ii = 0; ii < a; ++ii) {
        d[ii] = 1.0 / d[ii];
    }

    for (int ii = 0; ii < domod.quot + domod.rem; ii++) {
        gd[ii] = g[ii] * d[ii % a];
    }

    for (int ii = gl - 1, jj = a - 1; ii >= domod.quot + domod.rem;
         ii--, jj--) {
        if (jj < 0) {
            jj = a - 1;
        }
        gd[ii] = g[ii] * d[jj];
    }
}
