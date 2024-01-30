#ifndef RTPGHI_H__
#define RTPGHI_H__

#include <complex>

/**
 * Assert macro.
 */
#ifndef NDEBUG
    #define rtpghi_assert(cond, msg) \
        __rtpghi_assert(#cond, cond, __FILE__, __LINE__, msg);
#else
    #define rtpghi_assert(cond, msg)
#endif

/**
 * Implementation class of RTPGHI.
 */
class rtpghi_priv;

/**
 * Interface for using RTPGHI.
 */
class rtpghi_t final {
   public:
    /**
     * Create a RTPGHI state.
     *
     * \param[in]   W       Number of channels
     * \param[in]   a       Hop size
     * \param[in]   M       Number of frequency channels (FFT length)
     * \param[in]   tol     Relative coefficient tolerance
     */
    rtpghi_t(int W, int a, int M, double tol);

    /**
     * Destroys a RTPGHI plan.
     */
    ~rtpghi_t();

    /**
     * Reset RTPGHI state to the inital state.
     */
    void reset(const double** sinit);

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
     * The function is intended to be called for consecutive stream of frames as
     * it reuses some data from the previous frames stored in the plan.
     *
     * \param[in]   s       Target magnitude
     * \param[in]   stretch Stretch factor
     * \param[out]  c       Reconstructed coefficients
     */
    void execute(const std::complex<double>* s, double stretch,
                 std::complex<double>* c);

   private:
    rtpghi_priv* _p;
};

#ifndef NDEBUG
void __rtpghi_assert(const char* expr_str, bool expr, const char* file,
                     int line, const char* msg);
#endif

#endif  // RTPGHI_H__
