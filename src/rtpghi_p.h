#ifndef RTPGHI_P_H__
#define RTPGHI_P_H__

#define _USE_MATH_DEFINES 
#include <cmath>
#include <cstdint>

#include <memory>
#include <random>
#include <vector>

class rtpghi_update_plan;

class rtpghi_priv final {
public:
    rtpghi_priv(int W, int a, int M, double tol);
    ~rtpghi_priv();

    double get_stretch() const;
    void set_tolerance(double tol);

    void reset(const double **sinit);

    void execute(const std::complex<double> *cin, double stretch,
                 std::complex<double> *cout);


private:
    rtpghi_update_plan *_p;
    int _M;
    int _a;
    int _W;
    std::vector<double> _s;
    std::vector<double> _tgrad; //!< Time gradient buffer
    std::vector<double> _fgrad; //!< Frequency gradient buffer
    std::vector<double> _phase;
    std::vector<double> _phasein;
    double _stretch;
};

class rtpghi_update_plan {
public:
    rtpghi_update_plan(int M, int W, double tol);

    void execute_with_mask(const double *s,
                           const double *tgrad,
                           const double *fgrad,
                           const double *startphase,
                           const int8_t *mask,
                           double *phase);

    void execute(const double *s,
                 const double *tgrad,
                 const double *fgrad,
                 const double *startphase,
                 double *phase);

private:
    void execute_common(const double *s,
                        const double *tgrad,
                        const double *fgrad,
                        const double *startphase,
                        double *phase);

    double random_phase();

    std::random_device _rd;
    std::mt19937 _rand;
    std::uniform_real_distribution<> _dist;

    rtpghi_heap *_h;
    std::vector<int8_t> _donemask;
    double _tol;
    int _M;
};

#endif // RTPGHI_P_H__
