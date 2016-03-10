#include <PCU.h>
#include <stdio.h>
#include <limits.h>
#include "parma.h"
#include "parma_balancer.h"
#include "parma_step.h"
#include "parma_sides.h"
#include "parma_weights.h"
#include "parma_targets.h"
#include "parma_selector.h"
#include "parma_commons.h"
#include "parma_convert.h"

namespace {
  using parmaCommons::status;

  int getMinSide(parma::Sides* s) {
    int minSides = INT_MAX;
    s->begin();
    const parma::Sides::Item* side;
    while( (side = s->iterate()) ) 
      if(side->second < minSides && side->second>0) {
        minSides = side->second;
      }
    s->end();
    PCU_Debug_Print("minside %d\n", minSides);
    return minSides;
  }
  int getSmallestSide(parma::Sides* s) {
    int minSides=getMinSide(s);
    minSides = PCU_Min_Int(minSides);
    return minSides;
  }

  int getMisNumber(apf::Mesh* m) {
    const double t0 = PCU_Time();
    int misNumber = Parma_MisNumbering(m,0);
    double elapsedTime = PCU_Max_Double(PCU_Time()-t0);
    if( !PCU_Comm_Self() )
      status("mis completed in %f (seconds)\n", elapsedTime);
    return misNumber;
  }

  class ImbOrLong : public parma::Stop {
    public:
      ImbOrLong(parma::Sides* s, int tol)
        : sides(s), sideTol(tol) {}
      bool stop(double imb, double maxImb) {
        const int small = getSmallestSide(sides);
        if (!PCU_Comm_Self())
          status("Smallest Side %d, Target Side %f\n", small, sideTol);
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
        smallestTgtSide = 10;
        iter=0;
        misNumber = 0;
        maxMis = 0;
        if (!PCU_Comm_Self())
          status("Factor %f Smallest target side %d\n",f, smallestTgtSide);
      }
      bool runStep(apf::MeshTag* wtag, double tolerance) {
        if (!PCU_Comm_Self())
          status("Iteration: %d\n",iter);
        parma::Sides* s = parma::makeVtxSides(mesh);
        parma::Weights* w =
          parma::makeEntWeights(mesh, wtag, s, mesh->getDimension());
        if( !iter ) {
          misNumber = getMisNumber(mesh);
          maxMis = PCU_Max_Int(misNumber);
          if (!PCU_Comm_Self())
            status("mis maxNum %d\n", maxMis);
        }
        const int small = getSmallestSide(s);
        parma::Targets* t = 
          parma::makeShapeTargets(s, small, misNumber==iter);
        iter++;
        if (iter>maxMis)
          iter=0;
        parma::Selector* sel = parma::makeShapeSelector(mesh, wtag);
        ImbOrLong* stopper = new ImbOrLong(s, smallestTgtSide);
        parma::Stepper b(mesh, factor, s, w, t, sel, "elm", stopper);
        return b.step(tolerance, verbose);
      }
    private:
      int smallestTgtSide;
      int iter;
      int misNumber;
      int maxMis;
  };
}

apf::Balancer* Parma_MakeShapeOptimizer(apf::Mesh* m,
    double stepFactor, int verbosity) {
  return new ShapeOptimizer(m, stepFactor, verbosity);
}
