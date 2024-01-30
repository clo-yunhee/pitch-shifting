#include "rtpghi.h"

#include "rtpghi_p.h"

rtpghi_t::rtpghi_t(int W, int a, int M, double tol) {
    _p = new rtpghi_priv(W, a, M, tol);
}

rtpghi_t::~rtpghi_t() { delete _p; }

void rtpghi_t::reset(const double** sinit) { _p->reset(sinit); }

double rtpghi_t::get_stretch() const { return _p->get_stretch(); }

void rtpghi_t::set_tolerance(double tol) { _p->set_tolerance(tol); }

void rtpghi_t::execute(const std::complex<double>* s, double stretch,
                       std::complex<double>* c) {
    _p->execute(s, stretch, c);
}
