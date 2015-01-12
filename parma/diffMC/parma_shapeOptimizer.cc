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

  class ImbOrLong : public parma::Stop {
    public:
      ImbOrLong(parma::Sides* s, double tol)
        : sides(s), sideTol(tol) {}
      bool stop(double imb, double maxImb) {
        const double small = static_cast<double>(getSmallestSide(sides));
        if (!PCU_Comm_Self())
          fprintf(stdout,"Smallest Side: %f, endPoint: %f\n", small, sideTol);
        return imb > maxImb || small > sideTol;
      }
    private:
      parma::Sides* sides;
      double sideTol;
  };
  
  class ShapeOptimizer : public parma::Balancer {
    public:
      ShapeOptimizer(apf::Mesh* m, double f, int v)
        : Balancer(m, f, v, "gap") {
          parma::Sides* s = parma::makeVtxSides(mesh);
          avgSide=getAvgSides(s);
	  avgSideMult=0.4;
	  iter=0;
	  
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
	if (iter==0) {
	  Parma_ProcessDisconnectedParts(mesh);
	  const double t1 = PCU_Time();
	  misNumber = Parma_MisNumbering(mesh,0);
	  double elapsedTime = PCU_Time() - t1;
	  PCU_Max_Doubles(&elapsedTime, 1);
	  if( !PCU_Comm_Self() )
	    fprintf(stdout,"mis completed in %f (seconds)\n", elapsedTime);
	  maxMis = misNumber;
          PCU_Max_Ints(&maxMis,1);
	  avgSideMult+=.1;
        }
        parma::Targets* t = 
          parma::makeShapeTargets(mesh, s, w, factor, avgSideMult, misNumber==iter);
	iter++;
	if (iter>maxMis)
	  iter=0;
        PCU_Debug_Print("%s\n", t->print("targets").c_str());
        parma::Centroids c(mesh, wtag, s);
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag, &c);
        ImbOrLong* stopper = new ImbOrLong(s, avgSide*0.7);
        parma::Stepper b(mesh, factor, s, w, t, sel, stopper);
        return b.step(tolerance, verbose);
      }
    private:
      double avgSide;
      double avgSideMult;
      int iter;
      int misNumber;
      int maxMis;
    
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}
