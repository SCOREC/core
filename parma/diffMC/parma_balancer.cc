#include <PCU.h>
#include "parma_balancer.h"
#include "parma_monitor.h"
#include "parma_graphDist.h"

namespace {
  void printTiming(const char* type, int steps, double tol, double time) {
    if (!PCU_Comm_Self())
      printf("%s balanced in %d steps to %f in %f seconds\n",
          type, steps, tol, time);
  }
}

namespace parma {
  Balancer::Balancer(apf::Mesh* m, double f, int v, const char* n)
    : mesh(m), factor(f), verbose(v), name(n) {
      maxStep = 300;
      iS = new parma::Slope();
      iA = new parma::Average(8);
      sS = new parma::Slope();
      sA = new parma::Average(8);
  }
  Balancer::~Balancer() {
    delete iA;
    delete iS;
    delete sA;
    delete sS;
    apf::MeshTag* dist = parma::getDistTag(mesh);
    if(dist) {
      apf::removeTagFromDimension(mesh,dist,0);
      mesh->destroyTag(dist);
    }
  }
  void Balancer::balance(apf::MeshTag* wtag, double tolerance) {
    int step = 0;
    double t0 = PCU_Time();
    while (runStep(wtag,tolerance) && step++ < maxStep);
    printTiming(name, step, tolerance, PCU_Time()-t0);
  }
  void Balancer::monitorUpdate(double v, Slope* s, Average* a) {
    s->push(v);
    const double slope = (s->full()) ? s->slope() : 1.0;
    a->push(slope);
  }
}
