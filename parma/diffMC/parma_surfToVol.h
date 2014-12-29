#ifndef PARMA_SURFTOVOL_H
#define PARMA_SURFTOVOL_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides;
  class SurfToVol : public Associative<double> {
    public:
      SurfToVol() {};
      virtual ~SurfToVol() {}
      virtual double self() { return surfToVolImb; }
    protected:
      double surfToVolImb;
  };
  SurfToVol* makeSidesToElements(apf::Mesh* m, Sides* s);
}

#endif

