#ifndef PARMA_BALANCER_H
#define PARMA_BALANCER_H

#include <apfPartition.h>

namespace parma {
  class Slope;
  class Average;
  class Balancer : public apf::Balancer {
    public:
      Balancer(apf::Mesh* m, double f, int v, const char* n);
      ~Balancer();
      virtual bool runStep(apf::MeshTag* wtag, double tolerance)=0;
      virtual void balance(apf::MeshTag* wtag, double tolerance);
      void monitorUpdate(double v, Slope* s, Average* a);
      apf::Mesh* mesh;
      double factor;
      int verbose;
      const char* name;
      int maxStep;
    protected:
      Slope* iS;
      Average* iA;
      Slope* sS;
      Average* sA;
  };
}
#endif
