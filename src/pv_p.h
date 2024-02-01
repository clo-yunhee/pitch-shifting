#ifndef PV_P_H__
#define PV_P_H__

#include <memory>

#include "rtdgtrealproc.h"
#include "rtpghi.h"

class pv_priv final {
   public:
    pv_priv(double stretchmax, int Wmax, int bufLenMax);

    void print_pos() const;

    size_t next_inlen(size_t Lout) const;

    size_t next_outlen(size_t Lin) const;

    void advance_by(size_t Lin, size_t Lout);

    void set_stretch(double stretch);

    void execute(const double* in[], int Lin, int chan, double stretch,
                 int Lout, double* out[]);

    void execute_compact(const double* in, int Lin, int chan, double stretch,
                         int Lout, double* out);

   private:
    std::unique_ptr<rtdgtreal_processor_t> _proc;
    std::unique_ptr<rtpghi_t>              _rtpghi;
    double                                 _stretch;
    size_t                                 _in_pos;
    size_t                                 _out_pos;
    double                                 _in_in_out_offset;
    double                                 _out_in_in_offset;
    int                                    _aana;
    int                                    _asyn;

    friend void rtpghi_processor_callback(void*                       userdata,
                                          const std::complex<double>* in,
                                          int M2, int W,
                                          std::complex<double>* out);
};

#endif  // PV_P_H__