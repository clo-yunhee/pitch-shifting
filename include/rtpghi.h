#ifndef RTPGHI_H__
#define RTPGHI_H__

#include <complex>

/**
 * Implementation class of RTPGHI.
 */
class rtpghi_priv;

/**
 * Interface for using RTPGHI
 */
class rtpghi_t {
public:
    /**
     * Create a RTPGHI state.
     *
     * \param[in]   W       Number of channels
     * \param[in]   a       Hop size
     * \param[in]   M       Number of frequency channels (FFT length)
     * \param[in]   tol     Relative coefficient tolerance
     */
    rtpghi(int W, int a, int M, double tol);

    /**
     * Destroys a RTPGHI plan.
     */
    ~rtpghi();

    /**
     * Reset RTPGHI state to the inital state.
     */
    void reset();

    /**
     * Change tolerance.
     *
     * \param[in]   tol     Relative coefficient tolerance
     */
    void set_tolerance(double tol);

    double get_stretch() const;
    
    /**
     * Execute RTPGHI plan for a single frame.
     * 
     * The function is intended to be called for consecutive stream of frames as it
     * reuses some data from the previous frames stored in the plan.
     *
     * \param[in]   s       Target magnitude
     * \param[in]   stretch Stretch factor
     * \param[out]  c       Reconstructed coefficients
     */
    void execute(const std::complex<double> *s, double stretch, std::complex<double> *c);

private:
    rtpghi_priv *_p;
    
};

#endif RTPGHI_H__
