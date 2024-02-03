#include "circularbuf.h"

#include "rtpghi.h"

analysis_fifo_t::analysis_fifo_t(int fifoLen, int procDelay, int winLen,
                                 int hop, int numChans) {
    rtpghi_assert(fifoLen > 0, "fifoLen must be positive");
    rtpghi_assert(winLen > 0, "winLen must be positive");
    rtpghi_assert(hop > 0, "hop must be positive");
    rtpghi_assert(numChans > 0, "numChans must be positive");
    rtpghi_assert(procDelay >= winLen - 1,
                  "procDelay must be greater or equal than winLen-1");
    rtpghi_assert(fifoLen > winLen + 1,
                  "fifoLen must be greater than winLen+1");

    _winLen = winLen;
    _readchanstride = winLen;
    _hop = hop;
    _buf = std::make_unique<double[]>(numChans * (fifoLen + 1));
    _bufLen = fifoLen + 1;
    _readIdx = fifoLen + 1 - procDelay;
    _writeIdx=0;
    _numChans = numChans;
}

int analysis_fifo_t::get_numchans() const { return _numChans; }

void analysis_fifo_t::reset() {
    std::fill(_buf.get(), _buf.get() + _numChans * _bufLen, 0);
}

void analysis_fifo_t::set_hop(int hop) {
    rtpghi_assert(hop > 0, "hop must be positive");
    _hop = hop;
}

void analysis_fifo_t::set_readchanstride(int stride) {
    rtpghi_assert(stride > 0, "stride must be positive");
    _readchanstride = stride;
}

int analysis_fifo_t::write(const double** buf, int bufLen, int W) {
    int Wact, freeSpace, toWrite, valid, over, endWriteIdx;

    rtpghi_assert(bufLen >= 0, "bufLen must be positive");
    rtpghi_assert(W > 0, "W must be positive");

    if (bufLen == 0) return 0;

    freeSpace = _readIdx - _writeIdx - 1;
    if (freeSpace < 0) freeSpace += _bufLen;

    Wact = _numChans < W ? _numChans : W;

    toWrite = bufLen > freeSpace ? freeSpace : bufLen;
    valid = toWrite;
    over = 0;

    endWriteIdx = _writeIdx + toWrite;

    if (endWriteIdx > _bufLen) {
        valid = _bufLen - _writeIdx;
        over = endWriteIdx - _bufLen;
    }

    if (valid > 0) {
        for (int w = 0; w < _numChans; ++w) {
            double* pbufchan = _buf.get() + w * _bufLen + _writeIdx;
            if (w < Wact)
                std::copy(buf[w], buf[w] + valid, pbufchan);
            else
                std::fill(pbufchan, pbufchan + valid, 0);
        }
    }
    if (over > 0) {
        for (int w = 0; w < Wact; ++w) {
            double* pbufchan = _buf.get() + w * _bufLen;
            if (w < Wact)
                std::copy(buf[w] + valid, buf[w] + valid + over, pbufchan);
            else
                std::fill(pbufchan, pbufchan + over, 0);
        }
    }
    _writeIdx = (_writeIdx + toWrite) % _bufLen;

    return toWrite;
}

int analysis_fifo_t::read(double* buf) {
    int available, toRead, valid, over, endReadIdx;

    available = _writeIdx - _readIdx;
    if (available < 0) available += _bufLen;

    if (available < _winLen || available < _hop) return 0;

    toRead = _winLen;

    valid = toRead;
    over = 0;

    endReadIdx = _readIdx + valid;

    if (endReadIdx > _bufLen) {
        valid = _bufLen - _readIdx;
        over = endReadIdx - _bufLen;
    }

    if (valid > 0) {
        for (int w = 0; w < _numChans; ++w) {
            double* pbufchan = _buf.get() + w * _bufLen + _readIdx;
            std::copy(pbufchan, pbufchan + valid,
                      _buf.get() + w * _readchanstride);
        }
    }
    if (over > 0) {
        for (int w = 0; w < _numChans; ++w) {
            std::copy(_buf.get() + w * _bufLen, _buf.get() + w * _bufLen + over,
                      _buf.get() + valid + w * _readchanstride);
        }
    }

    // Only advance by hop
    _readIdx = (_readIdx + _hop) % _bufLen;

    return toRead;
}

synthesis_fifo_t::synthesis_fifo_t(int fifoLen, int winLen, int hop,
                                   int numChans) {
    rtpghi_assert(fifoLen > 0, "fifoLen must be positive");
    rtpghi_assert(winLen > 0, "winLen must be positive");
    rtpghi_assert(hop > 0, "hop must be positive");
    rtpghi_assert(numChans > 0, "numChans must be positive");
    rtpghi_assert(fifoLen > winLen + 1,
                  "fifoLen must be greater than winLen+1");

    _winLen = winLen;
    _writechanstride = winLen;
    _hop = hop;
    _buf = std::make_unique<double[]>(numChans * (fifoLen + winLen + 1));
    _bufLen = fifoLen + winLen + 1;
    _readIdx=0;
    _writeIdx=0;
    _numChans = numChans;
}

int synthesis_fifo_t::get_numchans() const { return _numChans; }

void synthesis_fifo_t::reset() {
    std::fill(_buf.get(), _buf.get() + _numChans * _bufLen, 0);
}

void synthesis_fifo_t::set_hop(int hop) {
    rtpghi_assert(hop > 0, "hop must be positive");
    _hop = hop;
}

void synthesis_fifo_t::set_writechanstride(int stride) {
    rtpghi_assert(stride > 0, "stride must be positive");
    _writechanstride = stride;
}

int synthesis_fifo_t::write(const double* buf) {
    int freeSpace, toWrite, valid, over, endWriteIdx;

    freeSpace = _readIdx - _writeIdx - 1;
    if (freeSpace < 0) freeSpace += _bufLen;

    if (freeSpace < _winLen) return 0;

    toWrite = _winLen;
    valid = toWrite;
    over = 0;

    endWriteIdx = _writeIdx + toWrite;

    if (endWriteIdx > _bufLen) {
        valid = _bufLen - _writeIdx;
        over = endWriteIdx - _bufLen;
    }

    if (valid > 0) {
        for (int w = 0; w < _numChans; ++w) {
            double* pbufchan = _buf.get() + _writeIdx + w * _bufLen;
            std::copy(pbufchan, pbufchan + valid,
                      _buf.get() + w * _writechanstride);
        }
    }
    if (over > 0) {
        for (int w = 0; w < _numChans; ++w) {
            double* pbufchan = _buf.get() + w * _bufLen;
            std::copy(pbufchan, pbufchan + over,
                      _buf.get() + valid + w * _writechanstride);
        }
    }

    _writeIdx = (_writeIdx + _hop) % _bufLen;

    return toWrite;
}

int synthesis_fifo_t::read(int bufLen, int W, double** buf) {
    int available, toRead, valid, over, endReadIdx;

    rtpghi_assert(W > 0, "W must be positive");
    rtpghi_assert(bufLen >= 0, "bufLen must be nonnegative");

    if (bufLen == 0) return 0;

    available = _writeIdx - _readIdx;
    if (available < 0) available += _bufLen;

    toRead = available < bufLen ? available : bufLen;

    valid = toRead;
    over = 0;

    endReadIdx = _readIdx + valid;

    if (endReadIdx > _bufLen) {
        valid = _bufLen - _readIdx;
        over = endReadIdx - _bufLen;
    }

    // Set the just read samples to zero so that the values
    // are not used in write again
    if (valid > 0) {
        for (int w = 0; w < W; ++w) {
            double* pbufchan = _buf.get() + _readIdx + w * _bufLen; 
            std::copy(pbufchan, pbufchan + valid, buf[w]);
            std::fill(pbufchan, pbufchan + valid, 0);
        }
    }
    if (over > 0) {
        for (int w = 0; w < W; ++w) {
            double* pbufchan = _buf.get() + w * _bufLen;
            std::copy(pbufchan, pbufchan + over, buf[w] + valid);
            std::fill(pbufchan, pbufchan + over, 0);
        }
    }

    _readIdx = (_readIdx + toRead) % _bufLen;

    return toRead;
}