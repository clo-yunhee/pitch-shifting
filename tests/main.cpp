#include <sndfile-64.h>
#include <sndfile.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>

#include "pv.h"

#define DEFAULT_RATIO 0.5
#define MAX_RATIO 10
#define BUFFER_LEN 1024
#define INOUT_BUFFER_LEN (BUFFER_LEN * MAX_RATIO)
#define MAX_CHANNELS 1

int main(int argc, char **argv) {
    static double indata[INOUT_BUFFER_LEN];
    static double outdata[INOUT_BUFFER_LEN];

    SNDFILE *infile, *outfile;

    SF_INFO sfinfo;
    int     readcount;
    int     writecount;

    char  *end{};
    double ratio = (argc >= 2 ? std::strtod(argv[1], &end) : DEFAULT_RATIO);
    if (errno == ERANGE) {
        ratio = DEFAULT_RATIO;
    }

    const char *infilename = (argc >= 3 ? argv[1] : "input.wav");
    const char *outfilename = (argc >= 4 ? argv[2] : "output.wav");

    printf("ratio:%.2f, infile:%s, outfile%s\n", ratio, infilename,
           outfilename);

    memset(&sfinfo, 0, sizeof(sfinfo));

    if (!(infile = sf_open(infilename, SFM_READ, &sfinfo))) {
        printf("Not able to open input file %s.\n", infilename);
        puts(sf_strerror(nullptr));
        return 1;
    }

    if (sfinfo.channels > MAX_CHANNELS) {
        printf("Not able to process more than %d channels.\n", MAX_CHANNELS);
        sf_close(infile);
        return 1;
    }

    if (!(outfile = sf_open(outfilename, SFM_WRITE, &sfinfo))) {
        printf("Not able to open output file %s.\n", outfilename);
        puts(sf_strerror(nullptr));
        sf_close(infile);
        return 1;
    }

    int Wmax = sfinfo.channels;
    int bufLenMax = BUFFER_LEN;

    pv_t pv(MAX_RATIO, Wmax, bufLenMax);

    do {
        int outlen = BUFFER_LEN;
        int inlen = pv.next_inlen(outlen);
        int channels = sfinfo.channels;

        readcount = (int)sf_read_double(infile, indata, inlen);

        pv.execute_compact(indata, readcount, channels, ratio, outlen, outdata);

        sf_write_double(outfile, outdata, outlen);
    } while (readcount != 0);

    sf_close(infile);
    sf_close(outfile);

    return 0;
}