#include "rtdgtrealproc.h"

#include "rtdgtrealproc_p.h"

rtdgtreal_processor_t::rtdgtreal_processor_t(firwin_t win, int gl, int a, int M,
                                             int numChans, int bufLenMax,
                                             int procDelay) {
    _p = new rtdgtreal_processor_priv(win, gl, a, M, numChans, bufLenMax,
                                      procDelay);
}

rtdgtreal_processor_t::~rtdgtreal_processor_t() { delete _p; }

void rtdgtreal_processor_t::reset() { _p->reset(); }

void rtdgtreal_processor_t::set_anaa(int a) { _p->set_anaa(a); }

void rtdgtreal_processor_t::set_syna(int a) { _p->set_syna(a); }

void rtdgtreal_processor_t::set_callback(rtdgtreal_processor_callback *callback,
                                         void *userdata) {
    _p->set_callback(callback, userdata);
}

void rtdgtreal_processor_t::execute_compact(const double *in, int len,
                                            int chanNo, double *out) {
    _p->execute_compact(in, len, chanNo, out);
}

void rtdgtreal_processor_t::execute_gen_compact(const double *in, int inLen,
                                                int chanNo, int outLen,
                                                double *out) {
    _p->execute_gen_compact(in, inLen, chanNo, outLen, out);
}

void rtdgtreal_processor_t::execute(const double **in, int len, int chanNo,
                                    double **out) {
    _p->execute(in, len, chanNo, out);
}

void rtdgtreal_processor_t::execute_gen(const double **in, int inLen,
                                        int chanNo, int outLen, double **out) {
    _p->execute_gen(in, inLen, chanNo, outLen, out);
}