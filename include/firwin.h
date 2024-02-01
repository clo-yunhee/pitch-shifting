#ifndef FIRWIN_H__
#define FIRWIN_H__

/**
 * Codes for finite support windows
 */
enum firwin_t {
    FIRWIN_HANN,
    FIRWIN_HANNING = FIRWIN_HANN,
    FIRWIN_NUTTALL10 = FIRWIN_HANN,
    FIRWIN_SQRTHANN,
    FIRWIN_COSINE = FIRWIN_SQRTHANN,
    FIRWIN_SINE = FIRWIN_SQRTHANN,
    FIRWIN_HAMMING,
    FIRWIN_NUTTALL01,
    FIRWIN_SQUARE,
    FIRWIN_RECT = FIRWIN_SQUARE,
    FIRWIN_TRIA,
    FIRWIN_TRIANGULAR = FIRWIN_TRIA,
    FIRWIN_BARTLETT = FIRWIN_TRIA,
    FIRWIN_SQRTTRIA,
    FIRWIN_BLACKMAN,
    FIRWIN_BLACKMAN2,
    FIRWIN_NUTTALL,
    FIRWIN_NUTTALL12 = FIRWIN_NUTTALL,
    FIRWIN_OGG,
    FIRWIN_ITERSINE = FIRWIN_OGG,
    FIRWIN_NUTTALL20,
    FIRWIN_NUTTALL11,
    FIRWIN_NUTTALL02,
    FIRWIN_NUTTALL30,
    FIRWIN_NUTTALL21,
    FIRWIN_NUTTALL03,
    FIRWIN_TRUNCGAUSS01,
    FIRWIN_TRUNCGAUSS005,
};

/**
 * Convert lowercase string to a firwin enum
 * E.g. "hann" returns LTFAT_HANN etc.
 */
firwin_t str2firwin(const char *name);

/**
 * Get the array length for mtgauss()
 *
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in]   thr   Threshold ]0,1[ where to truncate
 */
int mtgausslength(int a, int M, double thr);

/** Creates real, whole-point symmetric, zero delay window.
 *
 * \param[in]   win  Window type
 * \param[in]   gl   Window length
 * \param[out]  g    Window
 */
void firwin(firwin_t win, int gl, double *g);

/** Truncated Gaussian window optimally matched to the lattice
 *
 * Computes peak-normalized Gaussian window scaled such that the time-frequency
 * support is circular with respect to the lattice given by parameters a and M.
 *
 * \param[in]     a   Time hop factor
 * \param[in]     M   Number of frequency channels
 * \param[in]   thr   Threshold ]0,1[ where to truncate
 * \param[out]    g   Window. ltfat_mtgausslength
 */
void mtgauss(int a, int M, double thr, double *g);

#endif  // FIRWIN_H__