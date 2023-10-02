#ifndef PARMA_MONITOR_H
#define PARMA_MONITOR_H
namespace parma {
  class CircBuffer {
    private:
      unsigned int len; /* length of queue */
      unsigned int next; /* currently available */
      unsigned int sz; /* stored entries */
      double* q;
    public:
      CircBuffer(unsigned int l);
      ~CircBuffer();
      void push(double v);
      double get(unsigned int item);
      unsigned int length();
      unsigned int size();
      bool full();
      double slope();
      double avg();
  };
  class Slope : public CircBuffer {
    public:
      Slope();
      double slope();
  };
  class Average : public CircBuffer {
    public:
      Average(unsigned int l);
      double avg();
  };
}
#endif
