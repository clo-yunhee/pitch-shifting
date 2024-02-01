#define _USE_MATH_DEFINES
#include "firwin.h"

#include <cmath>
#include <cstdlib>

#include "rtpghi.h"

#define FIRWIN_RESETCOUNTER()                                \
    do {                                                     \
        if (ii == domod.quot + domod.rem) posInt = startInt; \
    } while (false)

void firwin(firwin_t win, int gl, double *g) {
    double     step, startInt, posInt;
    std::div_t domod;

    rtpghi_assert(gl > 0, "gl must be positive");

    step = 1.0 / gl;
    // for gl even;
    startInt = -0.5;
    domod = std::div(gl, 2);

    if (domod.rem) startInt = -0.5 + step / 2.0;

    posInt = 0;

    switch (win) {
        case FIRWIN_HANN: {
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.5 + 0.5 * cos(2.0 * M_PI * posInt);
                posInt += step;
            }
            break;
        }

        case FIRWIN_SQRTHANN:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = sqrt(0.5 + 0.5 * cos(2.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_HAMMING:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.54 + 0.46 * cos(2.0 * M_PI * posInt);
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL01:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.53836 + 0.46164 * cos(2 * M_PI * posInt);
                posInt += step;
            }
            break;

        case FIRWIN_RECT:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = (fabs(posInt) < 0.5 ? 1.0 : 0.0);
                posInt += step;
            }
            break;

        case FIRWIN_TRIANGULAR:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 1.0 - 2.0 * fabs(posInt);
                posInt += step;
            }
            break;

        case FIRWIN_SQRTTRIA:
            firwin(FIRWIN_TRIA, gl, g);
            for (int ii = 0; ii < gl; ++ii) g[ii] = sqrt(g[ii]);
            break;

        case FIRWIN_BLACKMAN:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.42 + 0.5 * cos(2 * M_PI * posInt +
                                         0.08 * cos(4.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_BLACKMAN2: {
            double denomfac = 1.0 / 18608.0;
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                double tmp = 7938.0 + 9240.0 * cos(2.0 * M_PI * posInt) +
                             1430.0 * cos(4.0 * M_PI * posInt);
                g[ii] = tmp * denomfac;
                posInt += step;
            }
            break;
        }
        case FIRWIN_NUTTALL:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.355768 +
                        0.487396 * cos(2.0 * M_PI * posInt +
                                       0.144232 * cos(4.0 * M_PI * posInt) +
                                       0.012604 * cos(6.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_OGG:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                double innercos = cos(M_PI * posInt);
                g[ii] = sin(M_PI / 2.0 * innercos * innercos);
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL20:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] =
                    (3.0 +
                     4.0 * cos(2.0 * M_PI * posInt + cos(4.0 * M_PI * posInt)) /
                         8.0);
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL11:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.40897 + 0.5 * cos(2.0 * M_PI * posInt +
                                            0.09103 * cos(4.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL02:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.4243801 +
                        0.4973406 * cos(2.0 * M_PI * posInt +
                                        0.0782793 * cos(4.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL30:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 10.0 + 15.0 * cos(2.0 * M_PI * posInt +
                                          6.0 * cos(4.0 * M_PI * posInt) +
                                          cos(6.0 * M_PI * posInt));
                g[ii] /= 32.0;
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL21:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.338946 +
                        0.481973 * cos(2.0 * M_PI * posInt +
                                       0.161054 * cos(4.0 * M_PI * posInt) +
                                       0.018027 * cos(6.0 * M_PI * posInt));
                posInt += step;
            }
            break;

        case FIRWIN_NUTTALL03:
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = 0.3635819 +
                        0.4891775 * cos(2.0 * M_PI * posInt +
                                        0.1365995 * cos(4.0 * M_PI * posInt) +
                                        0.0106411 * cos(6.0 * M_PI * posInt));
                posInt += step;
            }
            break;
        case FIRWIN_TRUNCGAUSS01: {
            double gamma = 4.0 * log(0.01);
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = exp(posInt * posInt * gamma);
                posInt += step;
            }
            break;
        }
        case FIRWIN_TRUNCGAUSS005: {
            double gamma = 4.0 * log(0.005);
            for (int ii = 0; ii < gl; ++ii) {
                FIRWIN_RESETCOUNTER();
                g[ii] = exp(posInt * posInt * gamma);
                posInt += step;
            }
            break;
        }
        default:
            fprintf(stderr, "Unknown window\n");
            abort();
    };

    // Fix symmetry of windows which are not zero at -0.5
    if (!domod.rem) g[domod.quot + domod.rem] = 0.0;
}

void mtgauss(int a, int M, double thr, double *g) {
    double     step, startInt, posInt, gamma;
    std::div_t domod;
    int        gl = mtgausslength(a, M, thr);

    step = 1.0 / gl;
    startInt = -0.5;
    domod = std::div(gl, 2);

    if (domod.rem) startInt = -0.5 + step / 2.0;

    posInt = 0;
    gamma = -M_PI * ((double)(gl * gl)) / ((double)(a * M));
    for (int ii = 0; ii < gl; ++ii) {
        FIRWIN_RESETCOUNTER();
        g[ii] = exp(posInt * posInt * gamma);
        posInt += step;
    }

    // Fix symmetry of windows which are not zero at -0.5
    if (!domod.rem) g[domod.quot + domod.rem] = 0.0;
}

firwin_t str2firwin(const char *win) {
    if (!strcmp("hann", win) || !strcmp("hanning", win) ||
        !strcmp("nuttall10", win))
        return FIRWIN_HANN;
    if (!strcmp("hamming", win))
        return FIRWIN_HAMMING;
    else if (!strcmp("sqrthann", win) || !strcmp("cosine", win) ||
             !strcmp("sine", win))
        return FIRWIN_SQRTHANN;
    else if (!strcmp("nuttall01", win))
        return FIRWIN_NUTTALL01;
    else if (!strcmp("square", win) || !strcmp("rect", win))
        return FIRWIN_SQUARE;
    else if (!strcmp("tria", win) || !strcmp("triangular", win) ||
             !strcmp("bartlett", win))
        return FIRWIN_TRIA;
    else if (!strcmp("sqrttria", win))
        return FIRWIN_SQRTTRIA;
    else if (!strcmp("blackman", win))
        return FIRWIN_BLACKMAN;
    else if (!strcmp("blackman2", win))
        return FIRWIN_BLACKMAN2;
    else if (!strcmp("nuttall", win) || !strcmp("nuttall12", win))
        return FIRWIN_NUTTALL;
    else if (!strcmp("ogg", win) || !strcmp("itersine", win))
        return FIRWIN_OGG;
    else if (!strcmp("nuttall20", win))
        return FIRWIN_NUTTALL20;
    else if (!strcmp("nuttall11", win))
        return FIRWIN_NUTTALL11;
    else if (!strcmp("nuttall02", win))
        return FIRWIN_NUTTALL02;
    else if (!strcmp("nuttall30", win))
        return FIRWIN_NUTTALL30;
    else if (!strcmp("nuttall21", win))
        return FIRWIN_NUTTALL21;
    else if (!strcmp("nuttall03", win))
        return FIRWIN_NUTTALL03;
    else if (!strcmp("truncgauss01", win))
        return FIRWIN_TRUNCGAUSS01;
    else if (!strcmp("truncgauss005", win))
        return FIRWIN_TRUNCGAUSS005;

    return (firwin_t)(-1);
}

int mtgausslength(int a, int M, double thr) {
    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(thr > 0.0 && thr <= 1.0, "thr must be in the range ]0,1[");

    return 2 * (int)round(sqrt(-a * M * log(thr) / M_PI));
}

#undef FIRWIN_RESETCOUNTER