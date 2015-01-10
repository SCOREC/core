#ifndef PARMA_SIDES_H
#define PARMA_SIDES_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides : public Associative<int> {
    public:
      Sides(apf::Mesh*) {}
      virtual ~Sides() {}
      virtual int total() { return totalSides;}
    protected:
      int totalSides;
  };
  Sides* makeElmBdrySides(apf::Mesh* m);
  Sides* makeElmSideSides(apf::Mesh* m);
  Sides* makeVtxSides(apf::Mesh* m);

  double avgSharedSides(Sides* s);
}

#endif
