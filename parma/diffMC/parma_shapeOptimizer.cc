#include <PCU.h>
#include <stdio.h>
#include <limits.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_centroids.h"
#include "parma_targets.h"
#include "parma_selector.h"

namespace {
  int getSmallestSide(parma::Sides* s) {
    int minSides = INT_MAX;
    s->begin();
    const parma::Sides::Item* side;
    while( (side = s->iterate()) ) 
      if(side->second < minSides ) {
	minSides = side->second;
      }
    s->end();
    PCU_Min_Ints(&minSides,1);
    return minSides;
  }
  
  double getAvgSides(parma::Sides* s) {
    double tot = s->total();
    PCU_Add_Doubles(&tot, 1);
    int cnt = static_cast<int>(s->size());
    PCU_Add_Ints(&cnt, 1);
    return tot/cnt;
  }
  class ShapeOptimizer : public parma::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "gap") { }
      static parma::Sides* s;
      static double avgSide;
      static bool greater(double imb, double maxImb) {
	int small = getSmallestSide(s);
	if (!PCU_Comm_Self())
	  fprintf(stdout,"Smallest Side: %d, endPoint: %f\n",small,avgSide*0.7);
        return imb > maxImb||small>0.7*avgSide;
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        s = parma::makeVtxSides(mesh);
        if (!avgSide) 
          avgSide=getAvgSides(s);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        parma::Targets* t = parma::makeShapeTargets(mesh, s, w, factor);
        PCU_Debug_Print("%s\n", t->print("targets").c_str());
        parma::Centroids c(mesh, wtag, s);
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag, &c);
        parma::Stepper b(mesh, factor, s, w, t, sel, greater);
        return b.step(tolerance, verbose);
      }
  };
  double ShapeOptimizer::avgSide=0;
  parma::Sides* ShapeOptimizer::s;
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}
