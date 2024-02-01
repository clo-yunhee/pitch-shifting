#include "rtdgtreal.h"

#include "rtdgtreal_p.h"

rtdgtreal_t::rtdgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype)
    : M(M) {
    _p = new rtdgtreal_priv(g, gl, M, ptype, DGT_FORWARD);
}

rtdgtreal_t::~rtdgtreal_t() { delete _p; }

void rtdgtreal_t::execute(const double *f, int W, std::complex<double> *c) {
    _p->execute_fwd(f, W, c);
}

rtidgtreal_t::rtidgtreal_t(const double *g, int gl, int M, rtdgt_phase_t ptype)
    : M(M) {
    _p = new rtdgtreal_priv(g, gl, M, ptype, DGT_INVERSE);
}

rtidgtreal_t::~rtidgtreal_t() { delete _p; }

void rtidgtreal_t::execute(const std::complex<double> *c, int W, double *f) {
    _p->execute_inv(c, W, f);
}