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
#include "maDBG.h"
#include <pcu_util.h>

namespace ma {

void adapt(Input* in, bool verbose)
{
  print("version 2.0 !");
  if (verbose)
    PCU_ALWAYS_ASSERT_VERBOSE(PCU_Comm_Peers() == 1,
    	"Verbose mesh adapt only works for 1 part meshes. Aborting!");
  double t0 = PCU_Time();
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    print("iteration %d",i);
    coarsen(a);
    if (verbose) ma_dbg::dumpMeshWithQualities(a,i,"after_coarsening");
    coarsenLayer(a);
    midBalance(a);
    refine(a);
    if (verbose) ma_dbg::dumpMeshWithQualities(a,i,"after_refining");
    snap(a);
    if (verbose) ma_dbg::dumpMeshWithQualities(a,i,"after_snapping");
    alignElements(a, i, true);
  }
  allowSplitCollapseOutsideLayer(a);
  fixElementShapes(a,verbose);
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
