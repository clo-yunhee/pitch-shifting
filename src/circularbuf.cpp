#include "circularbuf.h"

#include "circularbuf_p.h"

analysis_fifo_t::analysis_fifo_t(int fifoLen, int procDelay, int winLen,
                                 int hop, int numChans) {
    _p = new analysis_fifo_priv(fifoLen, procDelay, winLen, hop, numChans);
}