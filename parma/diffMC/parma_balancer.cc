#include "parma_balancer.h"
#include "parma_monitor.h"
#include "parma_graphDist.h"
#include "parma_commons.h"

namespace {
  void printTiming(const char* type, int steps, double tol, double time, pcu::PCU *PCUObj) {
    if (!PCUObj->Self())
      parmaCommons::status("%s balanced in %d steps to %f in %f seconds\n",
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
    if( 1 == mesh->getPCU()->Peers() ) return;
    int step = 0;
    double t0 = pcu::Time();
    while (runStep(wtag,tolerance) && step++ < maxStep);
    printTiming(name, step, tolerance, pcu::Time()-t0, mesh->getPCU());
  }
  void Balancer::monitorUpdate(double v, Slope* s, Average* a) {
    s->push(v);
    const double slope = (s->full()) ? s->slope() : 1.0;
    a->push(slope);
  }
}
