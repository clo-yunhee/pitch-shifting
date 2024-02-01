#ifndef RTPGHI_P_H__
#define RTPGHI_P_H__

#include <fftw3.h>

#include "rtdgtreal.h"

enum dgt_transformdirection_t {
    DGT_FORWARD,
    DGT_INVERSE,
};

class rtdgtreal_priv final {
   public:
    rtdgtreal_priv(const double *g, int gl, int M, const rtdgt_phase_t ptype,
                   const dgt_transformdirection_t tradir);

    ~rtdgtreal_priv();

    // forward
    void execute_fwd(const double *f, int W, std::complex<double> *c);

    // inverse
    void execute_inv(const std::complex<double> *c, int W, double *f);

   private:
    double               *_g;           //!< Window
    int                   _gl;          //!< Window length
    int                   _M;           //!< Number of FFT channels
    rtdgt_phase_t         _ptype;       //!< Phase convention
    double               *_fftBuf;      //!< Internal buffer
    std::complex<double> *_fftBuf_cpx;  //!< Internal buffer
    int                   _fftBufLen;   //!< Internal buffer length
    fftw_plan             _pfft;        //!< FFTW plan
};

#endif  // RTPGHI_P_H__