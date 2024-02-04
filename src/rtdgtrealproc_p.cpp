#include "rtdgtrealproc_p.h"

#include <fftw3.h>

#include "circularbuf.h"
#include "gabdual_painless.h"
#include "rtdgtreal.h"
#include "rtpghi.h"

rtdgtreal_processor_priv::rtdgtreal_processor_priv(firwin_t win, int gl, int a,
                                                   int M, int numChans,
                                                   int bufLenMax,
                                                   int procDelay) {
    rtpghi_assert(gl > 0, "gl must be positive");
    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(numChans > 0, "numChans must be positive");

    _g = std::make_unique<double[]>(gl);
    _gd = std::make_unique<double[]>(gl);

    firwin(win, gl, _g.get());
    gabdual_painless(_g.get(), gl, a, M, _gd.get());
    init(_g.get(), gl, _gd.get(), gl, a, M, numChans, bufLenMax, procDelay);
}

void rtdgtreal_processor_priv::init(const double *ga, int gal, const double *gs,
                                    int gsl, int a, int M, int numChans,
                                    int bufLenMax, int procDelay) {
    int glmax;

    rtpghi_assert(gal > 0, "gal must be positive");
    rtpghi_assert(gsl > 0, "gsl must be positive");

    glmax = gal > gsl ? gal - 1 : gsl - 1;

    rtpghi_assert(
        procDelay >= glmax && procDelay <= glmax + bufLenMax,
        "procDelay must be at least the window length at most the bufLenMax");

    rtpghi_assert(a > 0, "a must be positive");
    rtpghi_assert(M > 0, "M must be positive");
    rtpghi_assert(numChans > 0, "numChans must be positive");
    rtpghi_assert(bufLenMax > 0, "bufLenMax, must be positive");

    _fftBufIn =
        std::make_unique<std::complex<double>[]>(numChans * (M / 2 + 1));
    _fftBufOut =
        std::make_unique<std::complex<double>[]>(numChans * (M / 2 + 1));

    _buf = std::make_unique<double[]>(numChans * gal);
    _inTmp.resize(numChans);
    _outTmp.resize(numChans);

    _fwdfifo = std::make_unique<analysis_fifo_t>(bufLenMax + gal, procDelay,
                                                 gal, a, numChans);
    _backfifo =
        std::make_unique<synthesis_fifo_t>(bufLenMax + gsl, gsl, a, numChans);

    _fwdplan = std::make_unique<rtdgtreal_t>(ga, gal, M, RTDGTPHASE_ZERO);
    _backplan = std::make_unique<rtidgtreal_t>(gs, gsl, M, RTDGTPHASE_ZERO);

    _bufLenMax = bufLenMax;
}

void rtdgtreal_processor_priv::reset() {
    _fwdfifo->reset();
    _backfifo->reset();
}

void rtdgtreal_processor_priv::set_anaa(int a) { _fwdfifo->set_hop(a); }

void rtdgtreal_processor_priv::set_syna(int a) { _backfifo->set_hop(a); }

void rtdgtreal_processor_priv::set_callback(
    rtdgtreal_processor_callback *callback, void *userdata) {
    _callback = callback;
    _userdata = userdata;
}

void rtdgtreal_processor_priv::execute_compact(const double *in, int len,
                                               int chanNo, double *out) {
    return execute_gen_compact(in, len, chanNo, len, out);
}

void rtdgtreal_processor_priv::execute_gen_compact(const double *in, int inLen,
                                                   int chanNo, int outLen,
                                                   double *out) {
    int chanLoc;

    chanLoc =
        chanNo > _fwdfifo->get_numchans() ? _fwdfifo->get_numchans() : chanNo;

    for (int w = 0; w < chanLoc; ++w) {
        _inTmp[w] = &in[w * inLen];
        _outTmp[w] = &out[w * outLen];
    }

    // Clear superfluous channels
    if (chanNo > chanLoc) {
        std::fill(out + chanLoc * outLen, out + chanNo * outLen, 0);
    }

    execute_gen(_inTmp.data(), inLen, chanLoc, outLen, _outTmp.data());
}

void rtdgtreal_processor_priv::execute(const double **in, int len, int chanNo,
                                       double **out) {
    execute_gen(in, len, chanNo, len, out);
}

void rtdgtreal_processor_priv::execute_gen(const double **in, int inLen,
                                           int chanNo, int outLen,
                                           double **out) {
    int samplesWritten = 0;
    int samplesRead = 0;

    rtdgtreal_processor_callback *callback = _callback;

    rtpghi_assert(inLen >= 0 && outLen >= 0, "len must be nonnegative");
    rtpghi_assert(chanNo >= 0, "chanNo must be nonnegative");

    if (chanNo == 0 || (inLen == 0 && outLen == 0)) return;

    if (chanNo > _fwdfifo->get_numchans()) {
        for (int w = _fwdfifo->get_numchans(); w < chanNo; ++w) {
            std::fill(out[w], out[w] + outLen, 0);
        }
        chanNo = _fwdfifo->get_numchans();
    }

    if (inLen > _bufLenMax) {
        inLen = _bufLenMax;
    }

    if (outLen > _bufLenMax) {
        for (int w = 0; w < chanNo; ++w) {
            std::fill(out[w] + _bufLenMax, out[w] + outLen, 0);
        }
        outLen = _bufLenMax;
    }

    // Get default processor if none was set
    if (!callback) {
        callback = [](void *, const std::complex<double> *in, int M2, int W,
                      std::complex<double> *out) {
            std::copy(in, in + W * M2, out);
        };
    }

    // Write new data
    samplesWritten = _fwdfifo->write(in, inLen, chanNo);

    // While there is new data in the input fifo
    while (_fwdfifo->read(_buf.get()) > 0) {
        // Transform
        _fwdplan->execute(_buf.get(), _fwdfifo->get_numchans(),
                          _fftBufIn.get());

        // Process
        callback(_userdata, _fftBufIn.get(), _fwdplan->M / 2 + 1,
                 _fwdfifo->get_numchans(), _fftBufOut.get());

        // Reconstruct
        _backplan->execute(_fftBufOut.get(), _backfifo->get_numchans(),
                           _buf.get());

        // Write (and overlap) to out fifo
        _backfifo->write(_buf.get());
    }

    // Read samples for output
    samplesRead = _backfifo->read(outLen, chanNo, out);

    if (samplesWritten != inLen) {
        fprintf(stderr, "Samples written != input length\n");
    }
    if (samplesRead != outLen) {
        fprintf(stderr, "Samples read != output length\n");
    }
}