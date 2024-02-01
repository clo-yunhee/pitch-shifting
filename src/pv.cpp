#include "pv.h"

#include "pv_p.h"

pv_t::pv_t(double stretchmax, int Wmax, int bufLenMax) {
    _p = new pv_priv(stretchmax, Wmax, bufLenMax);
}

pv_t::~pv_t() { delete _p; }

void pv_t::print_pos() const { _p->print_pos(); }

size_t pv_t::next_inlen(size_t Lout) const { return _p->next_inlen(Lout); }

size_t pv_t::next_outlen(size_t Lin) const { return _p->next_outlen(Lin); }

void pv_t::advance_by(size_t Lin, size_t Lout) { _p->advance_by(Lin, Lout); }

void pv_t::set_stretch(double stretch) { _p->set_stretch(stretch); }

void pv_t::execute(const double* in[], int Lin, int chan, double stretch,
                   int Lout, double* out[]) {
    _p->execute(in, Lin, chan, stretch, Lout, out);
}

void pv_t::execute_compact(const double* in, int Lin, int chan, double stretch,
                           int Lout, double* out) {
    _p->execute_compact(in, Lin, chan, stretch, Lout, out);
}