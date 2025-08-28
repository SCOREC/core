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
#include "maShapeNew.h"
#include <pcu_util.h>
#include <iostream>

namespace ma {

void printHistogramData(std::string name, std::vector<double> input, double min, double max, Mesh* m)
{
  const int nbins = 10;
  int count[nbins] = {0};
  const double bin_size = (max-min)/(nbins*1.0);
  double inputMax = 0;
  double inputMin = max;

  for (size_t i = 0; i < input.size(); ++i) {
    if (input[i] > inputMax)
      inputMax = input[i];
    if (input[i] < inputMin)
      inputMin = input[i];
    int bin = (int)std::round((input[i] - min)/bin_size);
    count[bin] += 1;
  }

  std::cout << name << "\n";
  printf("Min: %f, Max: %f\n", inputMin, inputMax);
  for (int i = 0; i < nbins; ++i)
    fprintf(stderr, "%d\n", count[i]);
}

void printHistogramStats(Adapt* a)
{
  std::vector<double> lengths;
  std::vector<double> qualities;
  ma::stats(a->mesh, a->input->sizeField, lengths, qualities, true);
  printHistogramData("quality", qualities, 0, 1, a->mesh);
  printHistogramData("lengths", lengths, 0, MAXLENGTH+1, a->mesh);
}

void adapt(Input* in)
{
  double t0 = pcu::Time();
  print(in->mesh->getPCU(), "version 2.0 !");
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  coarsenMultiple(a);
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    print(a->mesh->getPCU(), "iteration %d", i);
    midBalance(a);
    refine(a);
    snap(a);
    coarsenMultiple(a);
    coarsenLayer(a);
  }
  allowSplitCollapseOutsideLayer(a);
  fixElementShapesNew(a);
  cleanupLayer(a);
  tetrahedronize(a);
  printQuality(a);
  postBalance(a);
  Mesh* m = a->mesh;

  double t1 = pcu::Time();
  print(m->getPCU(), "mesh adapted in %f seconds", t1-t0);
  printHistogramStats(a);

  delete a;
  // cleanup input object and associated sizefield and solutiontransfer objects
  if (in->ownsSizeField)
    delete in->sizeField;
  if (in->ownsSolutionTransfer)
    delete in->solutionTransfer;
  delete in;
  apf::printStats(m);
}

void adapt(const Input* in)
{
  adapt(makeAdvanced(in));
}

void adaptVerbose(Input* in, bool verbose)
{
  double t0 = pcu::Time();
  print(in->mesh->getPCU(), "version 2.0 - dev !");
  validateInput(in);
  Adapt* a = new Adapt(in);
  preBalance(a);
  for (int i = 0; i < in->maximumIterations; ++i)
  {
    print(a->mesh->getPCU(), "iteration %d", i);
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
  print(a->mesh->getPCU(), "Maximum (metric) edge length in the mesh is %f", lMax);
  while (lMax > MAXLENGTH) {
    print(a->mesh->getPCU(), "%dth additional refine-snap call", count);
    refine(a);
    snap(a);
    lMax = ma::getMaximumEdgeLength(a->mesh, a->sizeField);
    count++;
    print(a->mesh->getPCU(), "Maximum (metric) edge length in the mesh is %f", lMax);
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
  print(m->getPCU(), "mesh adapted in %f seconds", t1-t0);
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
