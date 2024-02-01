#ifndef RTDGTREALPROC_P_H__
#define RTDGTREALPROC_P_H__

#include <memory>
#include <vector>

#include "circularbuf.h"
#include "rtdgtreal.h"
#include "rtdgtrealproc.h"

class rtdgtreal_processor_priv final {
   public:
    rtdgtreal_processor_priv(firwin_t win, int gl, int a, int M, int numChans,
                             int bufLenMax, int procDelay);

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
    void init(const double *ga, int gal, const double *gs, int gsl, int a,
              int M, int numChans, int bufLenMax, int procDelay);

    std::unique_ptr<double[]>               _g;
    std::unique_ptr<double[]>               _gd;
    std::unique_ptr<std::complex<double>[]> _fftBufIn;
    std::unique_ptr<std::complex<double>[]> _fftBufOut;
    std::unique_ptr<double[]>               _buf;
    std::vector<const double *>             _inTmp;
    std::vector<double *>                   _outTmp;
    std::unique_ptr<analysis_fifo_t>        _fwdfifo;
    std::unique_ptr<synthesis_fifo_t>       _backfifo;
    std::unique_ptr<rtdgtreal_t>            _fwdplan;
    std::unique_ptr<rtidgtreal_t>           _backplan;
    int                                     _bufLenMax;

    rtdgtreal_processor_callback *_callback;  //!< Custom processor callback
    void                         *_userdata;  //!< Callback data
};

#endif  // RTDGTREALPROC_P_H__