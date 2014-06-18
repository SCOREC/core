#ifndef PARMA_SIDES_H
#define PARMA_SIDES_H
#include <apfMesh.h>
#include "parma_associative.h"

namespace parma {
  class Sides : public Associative<int> {
    public:
      Sides(apf::Mesh* m);
      virtual int total()=0;
    private:
      int totalSides;
  };
  Sides* makeElmBdrySides(apf::Mesh* m);
};

#endif
