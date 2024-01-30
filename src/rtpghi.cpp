#include "rtpghi.h"
#include "rtpghi_p.h"

rtpghi::rtpghi(int W, int a, int M, double tol)
{
    assert(W > 0, "W must be positive");
    assert(a > 0, "a must be positive");
    assert(M > 0, "M must be positive");
    assert(tol > 0 && tol < 1, "tol must be in range ]0,1[");

    _p = new rtpghi_priv(W, a, M, tol);
}

double rtpghi::get_stretch() const
{
    return _p->get_stretch();
}

void rtpghi::set_tolerance(double tol)
{
    assert(tol > 0 && tol < 1, "tol must be in range ]0,1[");
    _p->set_tolerance(tol);
}