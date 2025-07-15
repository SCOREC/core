#ifndef PARMA_STOP_H
#define PARMA_STOP_H

#include "parma_monitor.h"
#include "PCU.h"

namespace parma {
  class Stop {
    public:
      virtual ~Stop() {}
      virtual bool stop(double imb, double maxImb, pcu::PCU*)=0;
  };
  class Less : public Stop { 
    public:
      ~Less() {}
      bool stop(double imb, double maxImb, pcu::PCU*) {
        return imb < maxImb;
      }
  };
  class BalOrStall : public Stop {
    public:
      BalOrStall(Average* imb, Average* sides, double sidesTol, int verbose=0);
      ~BalOrStall() {}
      bool stop(double imb, double maxImb, pcu::PCU *PCUObj);
    private:
      Average* i;
      Average* s;
      double sTol;
      int verbose;
  };
}
#endif
