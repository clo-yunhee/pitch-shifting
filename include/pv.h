#ifndef PV_H__
#define PV_H__

#include <cstddef>

/**
 * Implementation class of PV.
 */
class pv_priv;

/**
 * Interface for using PV.
 */
class pv_t final {
   public:
    pv_t(double stretchmax, int Wmax, int buflenMax);

    ~pv_t();

    void print_pos() const;

    size_t next_inlen(size_t Lout) const;

    size_t next_outlen(size_t Lin) const;

    void advance_by(size_t Lin, size_t Lout);

    void set_stretch(double stretch);

    void execute(const double* in[], int Lin, int chan, double stretch,
                 int Lout, double* out[]);

    void execute_compact(const double* in, int Lin, int chan, double stretch,
                         int Lout, double* out);

   private:
    pv_priv* _p;
};

#endif  // PV_H__
