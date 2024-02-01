#ifndef RTDGTREAL_H__
#define RTDGTREAL_H__

#include <complex>

enum rtdgt_phase_t {
    RTDGTPHASE_ZERO,
    RTDGTPHASE_HALFSHIFT,
};

class rtdgtreal_priv;

// forward
class rtdgtreal_t final {
   public:
    rtdgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype);
    ~rtdgtreal_t();

    void execute(const double *f, int W, std::complex<double> *c);

    const int M;

   private:
    rtdgtreal_priv *_p;
};

// inverse
class rtidgtreal_t final {
   public:
    rtidgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype);
    ~rtidgtreal_t();

    void execute(const std::complex<double> *c, int W, double *f);

    const int M;

   private:
    rtdgtreal_priv *_p;
};

#endif  // RTDGTREAL_H__