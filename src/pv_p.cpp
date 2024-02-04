#include "pv_p.h"

#include <limits>

#include "firwin.h"
#include "pv.h"

void rtpghi_processor_callback(void *userdata, const std::complex<double> *in,
                               int M2, int W, std::complex<double> *out) {
    (void)M2;
    (void)W;
    auto state = static_cast<pv_priv *>(userdata);
    state->_rtpghi->execute(in, state->_stretch, out);
}

pv_priv::pv_priv(double stretchmax, int Wmax, int bufLenMax) {
    int asyn = 1024;
    int M = 8192;
    int gl = 4096;
    int fifoSize = (bufLenMax + asyn) * stretchmax;

    _procdelay = gl < fifoSize ? fifoSize : gl;

    _in_pos = 0;
    _out_pos = 0;
    _in_in_out_offset = 0.0;
    _out_in_in_offset = 0.0;

    _asyn = asyn;

    _proc = std::make_unique<rtdgtreal_processor_t>(FIRWIN_HANN, gl, asyn, M,
                                                    Wmax, fifoSize, _procdelay);

    _rtpghi = std::make_unique<rtpghi_t>(Wmax, asyn, M, 1e-6);

    _proc->set_callback(rtpghi_processor_callback, this);

    set_stretch(1.0);
}

int pv_priv::get_procdelay() const { return _procdelay; }

void pv_priv::print_pos() const {
    printf(
        "in_pos: %zu, out_pos: %zu, in_out_offset: %.3f, out_in_offset: %.3f, "
        "stretch: %.3f\n",
        _in_pos, _out_pos, _in_in_out_offset, _out_in_in_offset, _stretch);
}

size_t pv_priv::next_inlen(size_t Lout) const {
    double stretch = _rtpghi->get_stretch();
    size_t in_pos_end = (size_t)std::round(Lout / stretch + _out_in_in_offset);
    return in_pos_end;
}

size_t pv_priv::next_outlen(size_t Lin) const {
    double stretch = _rtpghi->get_stretch();
    size_t out_pos_end =
        (size_t)std::round((double)Lin * stretch + _in_in_out_offset);
    printf("stretch:%f, in_in_out_offset:%f\n", stretch, _in_in_out_offset);
    return out_pos_end;
}

void pv_priv::advance_by(size_t Lin, size_t Lout) {
    double stretch = _rtpghi->get_stretch();

    _in_pos += Lin;
    _out_pos += Lout;

    _in_in_out_offset += Lin * stretch;
    _in_in_out_offset -= Lout;

    _out_in_in_offset += Lout / stretch;
    _out_in_in_offset -= Lin;
}

void pv_priv::set_stretch(double stretch) {
    int    newaana = std::round(_asyn / stretch);
    double truestretch = ((double)_asyn) / newaana;

    if (fabs(truestretch - _stretch) > std::numeric_limits<double>::epsilon()) {
        _aana = newaana;
        _stretch = truestretch;
        _proc->set_anaa(newaana);
    }
}

void pv_priv::execute(const double *in[], int Lin, int chan, double stretch,
                      int Lout, double *out[]) {
    advance_by(Lin, Lout);
    set_stretch(stretch);
    _proc->execute_gen(in, Lin, chan, Lout, out);
}

void pv_priv::execute_compact(const double *in, int Lin, int chan,
                              double stretch, int Lout, double *out) {
    advance_by(Lin, Lout);
    set_stretch(stretch);
    _proc->execute_gen_compact(in, Lin, chan, Lout, out);
}