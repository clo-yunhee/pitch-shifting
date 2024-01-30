#ifndef CIRCULAR_BUF_H__
#define CIRCULAR_BUF_H__

class analysis_fifo_priv;
class synthesis_fifo_priv;

class analysis_fifo_t final {
   public:
    /** Create constant size output ring buffer
     *
     * The ring buffer works as usual when written to, but only constant
     * size (winLen) chunks can be read from it and the read pointer is
     * only advanced by hop after read.
     *
     * The buffer read and write pointers are initialized such that they
     * reflect the processing delay.
     *
     * \param[in]  fifoLen  Ring buffer size. This should be at least winLen +
     * max. expected buffer length. One more slot is actually allocated for the
     * "one slot open" implementation.
     * \param[in]  winLen   Window length
     * \param[in]  hop      Hop factor
     * \param[in]  numChans Maximum number of channels
     *
     * \returns RTDGTREAL_FIFO struct pointer
     */
    analysis_fifo_t(int fifoLen, int procDelay, int winLen, int hop,
                    int numChans);

    ~analysis_fifo_t();

    void reset();
    void set_hop(int hop);
    void set_read_chan_stride(int stride);

    /** Write bufLen samples to the analysis ring buffer
     *
     * The function returns number of samples written and a negative number if
     * something went wrong. If there is not enough space for all bufLen
     * samples, only available space is used and the number of actually written
     * samples is returned.
     *
     * \param[in]  p        Analysis ring buffer struct
     * \param[in]  buf      Channels to be written.
     * \param[in]  bufLen   Number of samples to be written
     * \param[in]  W        Number of channels
     *
     * \returns Number of samples written
     */
    int write(const double *buf[], int bufLen, int W);

    /** Read p->winLen samples from the analysis ring buffer
     *
     * The function attempts to read p->winLen samples from the buffer.
     *
     * The function does mothing and returns 0 if there is less than p->winLen
     * samples available.
     *
     * \param[in]   p        Analysis ring buffer struct
     * \param[out]  buf      Output array, it is expected to be able
     *                           to hold p->winLen * p->numChans samples.
     *
     * \returns Number of samples read
     */
    int read(double buf[]);

   private:
    analysis_fifo_priv *_p;
};

class synthesis_fifo_t final {
   public:
    /** Create constant size input ring buffer
     *
     * The ring buffer behaves as usual when read from, except it sets the read
     * samples to zero.
     * Only chunks of size winLen can be written to it and the write pointer is
     * advanced by hop. The samples are added to the existing values instead of
     * the usual overwrite.
     *
     * The buffer read and write pointers are both initialized to the same
     * value.
     *
     * \param[in]  fifoLen  Ring buffer size. This should be at least
                            winLen + max. expected buffer length.
                            (winLen+1) more slots are actually allocated
     *                      to accomodate the overlaps.
     * \param[in]  winLen   Window length
     * \param[in]  hop      Hop factor
     * \param[in]  numChans Maximum number of channels
     *
     * \returns RTIDGTREAL_FIFO struct pointer
     */
    synthesis_fifo_t(int fifoLen, int winLen, int hop, int numChans);

    ~synthesis_fifo_t();

    void reset();
    void set_hop(int hop);
    void set_write_chan_stride(int stride);

    /** Write p->winLen samples to DGT synthesis ring buffer
     *
     * The function returns 0 if there is not enough space to write all
     * p->winLen samples.
     *
     * \param[in]  p        Synthesis ring buffer struct
     * \param[in]  buf      Samples to be written
     *
     * \returns Number of samples written
     */
    int write(const double buf[]);

    /** Read bufLen samples from DGT analysis ring buffer
     *
     * The function attempts to read bufLen samples from the buffer.
     *
     * \param[in]   p        Analysis ring buffer struct
     * \param[in]   bufLen   Number of samples to be read
     * \param[in]   W        Number of channels
     * \param[out]  buf      Output channels, each channel is expected to be
     * able to hold bufLen samples.
     *
     * \returns Number of samples read
     */
    int read(int bufLen, int W, double *buf[]);

   private:
    synthesis_fifo_priv *_p;
};

#endif  // CIRCULAR_BUF_H__
