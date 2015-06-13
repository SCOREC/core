#include "parma_monitor.h"
#include <assert.h>
#include <stdlib.h>
namespace {
  const unsigned int order = 2;
  const double c[]={-3./2, 2., -1./2}; // 2nd order
}

namespace parma {
  CircBuffer::CircBuffer(unsigned int l) {
    len = l;
    q = static_cast<double*>(calloc(len, sizeof(double)));
    sz = next = 0;
  }
  CircBuffer::~CircBuffer() { free(q); }
  bool CircBuffer::full() { return (sz==len); }
  unsigned int CircBuffer::length() { return len; }
  unsigned int CircBuffer::size() { return sz; }
  double CircBuffer::get(unsigned int item) { 
    assert(item < sz);
    unsigned int tail = next;
    if( sz < len )
      tail = 0;
    unsigned int itemIdx = (tail + item) % len;
    return q[itemIdx];
  }
  void CircBuffer::push(double v) {
    q[next] = v;
    next = (next+1)%len;
    if( sz != len ) sz++;
  }
  Slope::Slope() : CircBuffer(order+1) {}
  double Slope::slope() {
    assert( full() );
    double s = 0;
    for(unsigned int i=0; i<length(); i++)
      s += c[i]*get(i);
    return s;
  }
  Average::Average(unsigned int l) : CircBuffer(l) {}
  double Average::avg() {
    double a = 0;
    for(unsigned int i=0; i<size(); i++)
      a += get(i);
    return a /= size();
  }
}
