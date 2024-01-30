#ifndef RTPGHI_HEAP_H__
#define RTPGHI_HEAP_H__

#include <vector>

class rtpghi_heap_t {
   public:
    rtpghi_heap_t(int initmaxsize, const double *s);

    const double *get_dataptr() const;

    void reset(const double *news);

    void push(int key);
    int  peek() const;
    int  pop();

   private:
    std::vector<int> _h;
    const double    *_s;
};

#endif  // RTPGHI_HEAP_H__