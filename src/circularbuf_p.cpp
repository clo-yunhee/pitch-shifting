#include "circularbuf_p.h"

#include "rtpghi.h"

analysis_fifo_priv::analysis_fifo_priv(int fifoLen, int procDelay, int winLen,
                                       int hop, int numChans) {
    rtpghi_assert(fifoLen > 0, "fifoLen must be positive");
    rtpghi_assert(winLen > 0, "winLen must be positive");
    rtpghi_assert(hop > 0, "hop must be positive");
    rtpghi_assert(numChans > 0, "numChans must be positive");
    rtpghi_assert(procDelay >= winLen - 1,
                  "procDelay must be greater or equal than winLen-1");
    rtpghi_assert(fifoLen > winLen + 1,
                  "fifoLen must be greater than winLen+1");

    _buf = std::make_unique<double[]>(numChans * (fifoLen + 1));
    _bufLen = fifoLen + 1;
    _hop = hop;
    _winLen = winLen;
    _readIdx = fifoLen + 1 - procDelay;
    _numChans = numChans;
    _readchanstride = winLen;
}

void analysis_fifo_priv::reset() {
    std::fill(_buf.get(), _buf.get() + _numChans * _bufLen, 0);
}

void analysis_fifo_priv::set_hop(int hop) {
    rtpghi_assert(hop > 0, "hop must be positive");
    _hop = hop;
}

void analysis_fifo_priv::set_readchanstride(int stride) {
    rtpghi_assert(stride > 0, "stride must be positive");
    _readchanstride = stride;
}

int analysis_fifo_priv::write(const double** buf, int bufLen, int W) {
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

int analysis_fifo_priv::read(double* buf) {
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