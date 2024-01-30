#include "rtpghi_heap.h"

rtpghi_heap_t::rtpghi_heap_t(int initmaxsize, const double *s) : _s(s) {
    _h.reserve(initmaxsize);
}

const double *rtpghi_heap_t::get_dataptr() const { return _s; }

void rtpghi_heap_t::reset(const double *news) { _s = news; }

void rtpghi_heap_t::push(int key) {
    int pos, pos2;

    pos = _h.size();
    _h.push_back(0);  // Let std::vector handle dynamic growth.

    double val = _s[key];

    while (pos > 0) {
        pos2 = (pos - 1) >> 1;

        if (_s[_h[pos2]] < val)
            _h[pos] = _h[pos2];
        else
            break;

        pos = pos2;
    }

    _h[pos] = key;
}

int rtpghi_heap_t::peek() const {
    if (_h.empty()) return -1;
    return _h[0];
}

int rtpghi_heap_t::pop() {
    int    pos, pos2, retkey, key;
    double maxchildkey, val;

    if (_h.empty()) return -1;

    /* Extract first element */
    retkey = _h[0];
    key = _h.back();
    val = _s[key];

    _h.pop_back();

    pos = 0;
    pos2 = 1;

    while (pos2 < _h.size()) {
        if ((pos2 + 2 > _h.size()) || (_s[_h[pos2]] >= _s[_h[pos2 + 1]])) {
            maxchildkey = _s[_h[pos2]];
        } else {
            maxchildkey = _s[_h[++pos2]];
        }

        if (maxchildkey > val)
            _h[pos] = _h[pos2];
        else
            break;

        pos = pos2;
        pos2 = (pos << 1) + 1;
    }

    _h[pos] = key;

    return retkey;
}