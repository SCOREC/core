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
#include "maSnap.h"
#include "maShape.h"
#include "maBalance.h"
#include "maLayer.h"
#include "maDBG.h"
#include <pcu_util.h>

namespace ma {

void adapt(Input* in)
{
  double t0 = pcu::Time();
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  print("version 2.0 !", a->mesh->getPCU());
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    print("iteration %d", a->mesh->getPCU(), i);
    coarsen(a);
    coarsenLayer(a);
    midBalance(a);
    refine(a);
    snap(a);
  }
  allowSplitCollapseOutsideLayer(a);
  fixElementShapes(a);
  cleanupLayer(a);
  tetrahedronize(a);
  printQuality(a);
  postBalance(a);
  Mesh* m = a->mesh;
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  double t1 = pcu::Time();
  print("mesh adapted in %f seconds", m->getPCU(), t1-t0);
  apf::printStats(m);
}

void adapt(const Input* in)
{
  adapt(makeAdvanced(in));
}

void adaptVerbose(Input* in, bool verbose)
{
  double t0 = pcu::Time();
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  print("version 2.0 - dev !", a->mesh->getPCU());
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    print("iteration %d",a->mesh->getPCU(), i);
    coarsen(a);
    if (verbose && in->shouldCoarsen)
      ma_dbg::dumpMeshWithQualities(a,i,"after_coarsen");
    coarsenLayer(a);
    midBalance(a);
    refine(a);
    if (verbose)
      ma_dbg::dumpMeshWithQualities(a,i,"after_refine");
    snap(a);
    if (verbose && in->shouldSnap)
      ma_dbg::dumpMeshWithQualities(a,i,"after_snap");
    fixElementShapes(a);
    if (verbose && in->shouldFixShape)
      ma_dbg::dumpMeshWithQualities(a,i,"after_fix");
  }
  allowSplitCollapseOutsideLayer(a);
  if (verbose) ma_dbg::dumpMeshWithQualities(a,999,"after_final_fix");
  // The following is applied to 2D surface meshes only and has no effect
  // on 3D meshes
  improveQualities(a);
  if (verbose)
    ma_dbg::dumpMeshWithQualities(a,999,"after_improveQualities");
  /* The following loop ensures that no long edges are left in
   * the mesh. Note that at this point all elements are of "good"
   * quality and a few refinement iterations will not depreciate
   * the overall quality of the mesh, significantly.
   */
  int count = 0;
  double lMax = ma::getMaximumEdgeLength(a->mesh, a->sizeField);
  print("Maximum (metric) edge length in the mesh is %f", a->mesh->getPCU(), lMax);
  while (lMax > 1.5) {
    print("%dth additional refine-snap call", a->mesh->getPCU(), count);
    refine(a);
    snap(a);
    lMax = ma::getMaximumEdgeLength(a->mesh, a->sizeField);
    count++;
    print("Maximum (metric) edge length in the mesh is %f", a->mesh->getPCU(), lMax);
    if (count > 5) break;
  }
  if (verbose)
    ma_dbg::dumpMeshWithQualities(a,999,"after_final_refine_snap_loop");
  printQuality(a);
  cleanupLayer(a);
  tetrahedronize(a);
  printQuality(a);
  postBalance(a);
  Mesh* m = a->mesh;
  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  double t1 = pcu::Time();
  print("mesh adapted in %f seconds", m->getPCU(), t1-t0);
  apf::printStats(m);
}

void adaptVerbose(const Input* in, bool verbose)
{
  adaptVerbose(makeAdvanced(in), verbose);
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
