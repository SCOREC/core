/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "ma.h"
#include "maAdapt.h"
#include "maCoarsen.h"
#include "maRefine.h"
#include "maShape.h"
#include "maBalance.h"
#include "maLayer.h"
#include <PCU.h>

namespace ma {

void adapt(Input* in)
{
  print("version 2.0 !");
  double t0 = MPI_Wtime();
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  preventChangesToLayer(a);
  allowSplitInLayer(a);
  for (int i=0; i < in->maximumIterations; ++i)
  {
    print("iteration %d",i);
    coarsen(a);
    midBalance(a);
    refine(a);
  }
  preventChangesToLayer(a);
  allowSplitCollapseOutsideLayer(a);
  fixElementShapes(a);
  turnLayerToTets(a);
  postBalance(a);
  Mesh* m = a->mesh;
  delete a;
  delete in;
  double t1 = MPI_Wtime();
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

void adapt(Mesh* m, apf::Field* f, SolutionTransfer* s)
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
