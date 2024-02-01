#include "arrayutils.h"

#include "rtpghi.h"

void circshift(const double *in, int L, int shift, double *out) {
    int p;
    rtpghi_assert(L > 0, "L must be positive");

    // Fix shift
    p = (L - shift) % L;

    if (p < 0) p += L;

    if (in == out) {
        if (p) {  // Do nothing if no shift is needed
            int m, count, i, j;

            for (m = 0, count = 0; count != L; ++m) {
                double t = in[m];

                for (i = m, j = m + p; j != m;
                     i = j, j = j + p < L ? j + p : j + p - L, ++count)
                    out[i] = out[j];

                out[i] = t;
                ++count;
            }
        }
    } else {
        // Still ok if p == 0
        std::copy(in + p, in + p + (L - p), out);
        std::copy(in, in + p, out + L - p);
    }
}

void fftshift(const double *in, int L, double *out) {
    circshift(in, L, (L / 2), out);
}

int positiverem(int a, int b) {
    int c = a % b;
    return (c < 0 ? c + b : c);
}

void clear_array(double *in, int L) {
    rtpghi_assert(L >= 0, "L must be nonnegative");
    if (L > 0) std::fill(in, in + L, 0);
}

void fold_array(const double *in, int Lin, int offset, int Lfold, double *out) {
    int startIdx;

    rtpghi_assert(Lin > 0, "Lin must be positive");
    rtpghi_assert(Lfold > 0, "Lfold must be positive");

    // Sanitize offset.
    startIdx = positiverem(offset, Lfold);

    // Clear output, we will use it as an accumulator
    if (in != out)
        clear_array(out, Lfold);
    else if (Lfold > Lin)
        clear_array(out + Lin, Lfold - Lin);

    if (!startIdx) {
        // Common code for no offset
        int startAt = (in == out) ? Lfold : 0;

        for (int ii = startAt; ii < Lin;)
            for (int kk = 0; ii < Lin && kk < Lfold; ++ii, ++kk)
                out[kk] += in[ii];
    } else {
        if (in == out) {
            // We cannot avoid the (slow) inplace circshift anyway.
            // Do it after the folding.
            fold_array(in, Lin, 0, Lfold, out);
            circshift(in, Lfold, startIdx, out);
        } else {
            // We avoid the inplace circshift by effectively
            // doing circshift of all blocks.
            for (int ii = 0; ii < Lin;) {
                int kk = startIdx;
                for (; kk < Lfold && ii < Lin; ++ii, ++kk) out[kk] += in[ii];

                for (kk = 0; kk < startIdx && ii < Lin; ++ii, ++kk)
                    out[kk] += in[ii];
            }
        }
    }
}

void periodize_array(const double *in, int Lin, int Lout, double *out) {
    rtpghi_assert(Lin > 0, "Lin must be positive");
    rtpghi_assert(Lout > 0, "Lout must be positive");

    // Do nothing if there is no place where to put periodized samples
    if (Lout <= Lin) {
        if (in != out) {
            std::copy(in, in + Lout, out);
        }
    } else {
        int periods = Lout / Lin;
        int lastL = Lout - periods * Lin;
        int startPer = (in == out) ? 1 : 0;

        for (int ii = startPer; ii < periods; ++ii)
            std::copy(in, in + Lin, out + ii * Lin);

        std::copy(in, in + lastL, out + periods * Lin);
    }
}