#ifndef CIRCULAR_BUF_P_H__
#define CIRCULAR_BUF_P_H__

#include <memory>

class analysis_fifo_priv {
   public:
    analysis_fifo_priv(int fifoLen, int procDelay, int winLen, int hop,
                       int numChans);

    void reset();
    void set_hop(int hop);
    void set_readchanstride(int stride);

    int write(const double** buf, int bufLen, int W);
    int read(double* buf);

   private:
    int                       _winLen;          //!< Window length
    int                       _readchanstride;  //!< Window length
    int                       _hop;             //!< Hop size
    std::unique_ptr<double[]> _buf;             //!< Ring buffer array
    int                       _bufLen;          //!< Length of the previous
    int                       _readIdx;         //!< Read pos.
    int                       _writeIdx;        //!< Write pos.
    int                       _numChans;
};

class synthesis_fifo_priv {
   public:
    synthesis_fifo_priv(int fifoLen, int winLen, int hop, int numChans);

    void reset();
    void set_hop(int hop);
    void set_writechanstride(int stride);

    int write(const double* buf);
    int read(int bufLen, int W, double** buf);

   private:
    int                       _winLen;           //!< Window length
    int                       _writechanstride;  //!< Window length
    int                       _hop;              //!< Hop size
    std::unique_ptr<double[]> _buf;              //!< Ring buffer array
    int                       _bufLen;           //!< Length of the previous
    int                       _readIdx;          //!< Read pos.
    int                       _writeIdx;         //!< Write pos.
    int                       _numChans;
};

#endif  // CIRCULAR_BUF_P_H__