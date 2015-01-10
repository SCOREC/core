#ifndef PARMA_STOP_H
#define PARMA_STOP_H

#include "parma_monitor.h"

namespace parma {
  class Stop {
    public:
      virtual ~Stop() {}
      virtual bool stop(double imb, double maxImb)=0;
  };
  class Less : public Stop { 
    public:
      ~Less() {}
      bool stop(double imb, double maxImb) {
        return imb < maxImb;
      }
  };
  class BalOrStall : public Stop {
    public:
      BalOrStall(Average* imb, Average* sides, double sidesTol);
      ~BalOrStall() {}
      bool stop(double imb, double maxImb);
    private:
      Average* i;
      Average* s;
      double sTol;
  };
}
#endif
