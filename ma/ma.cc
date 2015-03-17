/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include <PCU.h>
#include "ma.h"
#include "maAdapt.h"
#include "maCoarsen.h"
#include "maRefine.h"
#include "maSnap.h"
#include "maShape.h"
#include "maBalance.h"
#include "maLayer.h"

namespace ma {

void adapt(Input* in)
{
  print("version 2.0 !");
  double t0 = PCU_Time();
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  for (int i=0; i < in->maximumIterations; ++i)
  {
    print("iteration %d",i);
    coarsen(a);
    coarsenLayer(a);
    midBalance(a);
    refine(a);
  }
  allowSplitCollapseOutsideLayer(a);
  snap(a);
  fixElementShapes(a);
  cleanupLayer(a);
  tetrahedronize(a);
  printQuality(a);
  postBalance(a);
  Mesh* m = a->mesh;
  delete a;
  delete in;
  double t1 = PCU_Time();
  print("mesh adapted in %f seconds",t1-t0);
  apf::printStats(m);
}

void adapt(Mesh* m, AnisotropicFunction* f, SolutionTransfer* s)
{
  adapt(configure(m,f,s));
}

void adapt(Mesh* m, IsotropicFunction* f, SolutionTransfer* s)
{
  adapt(configure(m,f,s));
}

void runUniformRefinement(Mesh* m, int n, SolutionTransfer* s)
{
  adapt(configureUniformRefine(m,n,s));
}

void adaptMatching(Mesh* m, int n, SolutionTransfer* s)
{
  adapt(configureMatching(m,n,s));
}

}
