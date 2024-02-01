#ifndef RTDGTREALPROC_H__
#define RTDGTREALPROC_H__

#include <complex>

#include "firwin.h"

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

class rtdgtreal_processor_priv;

class rtdgtreal_processor_t final {
   public:
    rtdgtreal_processor_t(firwin_t win, int gl, int a, int M, int numChans,
                          int bufLenMax, int procDelay);

    ~rtdgtreal_processor_t();

    void reset();
    void set_anaa(int a);
    void set_syna(int a);
    void set_callback(rtdgtreal_processor_callback *callback, void *userdata);

    void execute_compact(const double *in, int len, int chanNo, double *out);

    void execute_gen_compact(const double *in, int inLen, int chanNo,
                             int outLen, double *out);

    void execute(const double **in, int len, int chanNo, double **out);

    void execute_gen(const double **in, int inLen, int chanNo, int outLen,
                     double **out);

   private:
    rtdgtreal_processor_priv *_p;
};

#endif  // RTDGTREALPROC_H__