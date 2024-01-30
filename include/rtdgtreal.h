#ifndef RTDGTREAL_H__
#define RTDGTREAL_H__

#include <complex>

#include "circularbuf.h"

enum rtdgt_phase_t {
    RTDGTPHASE_ZERO,
    RTDGTPHASE_HALFSHIFT,
};

enum dgt_transformdirection {
    DGT_FORWARD,
    DGT_INVERSE,
};

class rtdgtreal_priv;

// forward
class rtdgtreal_t final {
   public:
    rtdgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype);
    ~rtdgtreal_t();

    void execute(const double *f, int W, std::complex<double> *c);

   private:
    rtdgtreal_priv *_p;
};

// inverse
class rtidgtreal_t final {
   public:
    rtidgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype);
    ~rtidgtreal_t();

    void execute(const std::complex<double> *c, int W, double *f);

   private:
    rtdgtreal_priv *_p;
};

/**
 * \param[in]  userdata   User defined data
 * \param[in]        in   Input coefficients, M2 x W array
 * \param[in]        M2   Length of the arrays;
 *                            number of unique FFT channels; equals to M/2 + 1
 * \param[in]         W   Number of channels
 * \param[out]      out   Output coefficients, M2 x W array
 */
using rtdgtreal_processor_callback = void(void                       *userdata,
                                          const std::complex<double> *in,
                                          int M2, int W,
                                          std::complex<double> *out);

class rtdgtreal_processor_t {
    rtdgtreal_processor_callback *_callback;  //!< Custom processor callback
    void                         *userDdata;  //!< Callback data
    analysis_fifo_t              *fwdfifo;
    synthesis_fifo_t             *backfifo;
    rtdgtreal_t                  *fwdplan;
    rtidgtreal_t                 *backplan;
    double                       *buf;
    std::complex<double>         *fftBufIn;
    std::complex<double>         *fftBufOut;
    int                           bufLenMax;
    void                        **garbageBin;
    int                           garbageBinSize;
    const double                **inTmp;
    double                      **outTmp;
};

#endif  // RTDGTREAL_H__