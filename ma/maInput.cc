/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maInput.h"
#include "maAdapt.h"
#include <lionPrint.h>
#include <apfShape.h>
#include <cstdio>
#include <pcu_util.h>
#include <cstdlib>

namespace ma {

void setDefaultValues(Input* in)
{
  in->ownsSizeField = true;
  in->maximumIterations = 3;
  in->shouldCoarsen = true;
  in->shouldSnap = in->mesh->canSnap();
  in->shouldTransferParametric = in->mesh->canSnap();
  in->shouldTransferToClosestPoint = false;
  in->shouldHandleMatching = in->mesh->hasMatching();
  in->shouldFixShape = true;
  in->shouldForceAdaptation = false;
  in->shouldPrintQuality = true;
  if (in->mesh->getDimension()==3)
  {
    in->goodQuality = 0.027;
    in->maximumEdgeRatio = 2.0;
  }
  else
  { PCU_ALWAYS_ASSERT(in->mesh->getDimension()==2);
    //old MA says .04, but that rarely kicks in.
    //.2 is strict, but at least quality goes up
    in->goodQuality = 0.2;
    //this basically turns off short edge removal...
    in->maximumEdgeRatio = 100.0;
    //2D mesh adapt performs better if forceAdapt is on
    in->shouldForceAdaptation = true;
  }
  in->shouldCheckQualityForDoubleSplits = false;
  in->validQuality = 1e-10;
  in->maximumImbalance = 1.10;
  in->shouldRunPreZoltan = false;
  in->shouldRunPreZoltanRib = false;
  in->shouldRunPreParma = false;
  in->shouldRunMidZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostZoltan = false;
  in->shouldRunPostZoltanRib = false;
  in->shouldRunPostParma = false;
  in->shouldTurnLayerToTets = false;
  in->shouldCleanupLayer = false;
  in->shouldRefineLayer = false;
  in->shouldCoarsenLayer = false;
  in->splitAllLayerEdges = false;
  in->userDefinedLayerTagName = "";
  in->shapeHandler = 0;
  in->debugFolder = nullptr;
}

void rejectInput(const char* str, pcu::PCU *PCUObj)
{
  if (PCUObj->Self() != 0)
    return;
  lion_eprint(1,"MeshAdapt input error:\n");
  lion_eprint(1,"%s\n",str);
  abort();
}

// if more than 1 option is true return true
static bool moreThanOneOptionIsTrue(bool op1, bool op2, bool op3)
{
  int cnt = 0;
  if (op1) cnt++;
  if (op2) cnt++;
  if (op3) cnt++;
  if (cnt > 1)
    return true;
  return false;
}

void validateInput(Input* in)
{
  if ( ! in->sizeField)
    rejectInput("no size field", in->mesh->getPCU());
  if ( ! in->solutionTransfer)
    rejectInput("no solution transfer object", in->mesh->getPCU());
  if (in->maximumIterations < 0)
    rejectInput("negative maximum iteration count", in->mesh->getPCU());
  if (in->maximumIterations > 10)
    rejectInput("unusually high maximum iteration count", in->mesh->getPCU());
  if (in->shouldSnap
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested snapping "
                "but the geometric model does not support it", in->mesh->getPCU());
  if (in->shouldTransferParametric
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested parametric coordinate transfer "
                "but the geometric model does not support it", in->mesh->getPCU());
  if (in->shouldTransferToClosestPoint
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested transfer to closest point on model"
                "but the geometric model does not support it", in->mesh->getPCU());
  if (in->shouldSnap && ( ! (in->shouldTransferParametric ||
			     in->shouldTransferToClosestPoint)))
    rejectInput("snapping requires parametric coordinate transfer or transfer to closest point", in->mesh->getPCU());
  if ((in->mesh->hasMatching())
    &&( ! in->shouldHandleMatching))
    rejectInput("the mesh has matching entities but matched support is off", in->mesh->getPCU());
  if (in->shouldHandleMatching
    && in->shouldFixShape)
    rejectInput("user requested matched mesh handling and shape correction "
        "but shape correction does not support matching yet", in->mesh->getPCU());
  if (in->goodQuality < 0.0)
    rejectInput("negative desired element quality", in->mesh->getPCU());
  if (in->goodQuality > 1.0)
    rejectInput("desired element quality greater than one", in->mesh->getPCU());
  if (in->validQuality < 0.0)
    rejectInput("negative minimum element quality", in->mesh->getPCU());
  if (in->maximumImbalance < 1.0)
    rejectInput("maximum imbalance less than 1.0", in->mesh->getPCU());
  if (in->maximumEdgeRatio < 1.0)
    rejectInput("maximum tet edge ratio less than one", in->mesh->getPCU());
  if (moreThanOneOptionIsTrue(
  	in->shouldRunPreZoltan,
  	in->shouldRunPreZoltanRib,
  	in->shouldRunPreParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PreBalance options can be set to true!", in->mesh->getPCU());
  if (moreThanOneOptionIsTrue(
  	in->shouldRunPostZoltan,
  	in->shouldRunPostZoltanRib,
  	in->shouldRunPostParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PostBalance options can be set to true!", in->mesh->getPCU());
  if (in->shouldRunMidZoltan && in->shouldRunMidParma)
    rejectInput("only one of Zoltan and Parma MidBalance options can be set to true!", in->mesh->getPCU());
#ifndef PUMI_HAS_ZOLTAN
  if (in->shouldRunPreZoltan ||
      in->shouldRunPreZoltanRib ||
      in->shouldRunMidZoltan)
    rejectInput("core is not compiled with Zoltan. Use a different balancer or compile core with ENABLE_ZOLTAN=ON!", in->mesh->getPCU());
#endif
}

static void updateMaxIterBasedOnSize(Mesh* m, Input* in)
{
  // number of iterations
  double maxMetricLength = getMaximumEdgeLength(m, in->sizeField);
  int iter = std::ceil(std::log2(maxMetricLength));
  if (iter >= 10) {
    print(m->getPCU(), "ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations,\n"
    	"           which is equal to or larger than the maximum of 10 allowed.\n"
    	"           Setting the number of iteration to 10!", iter);
    in->maximumIterations = 10;
  }
  else {
    print(m->getPCU(), "ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations.\n"
    	"           Setting the number of iteration to %d!", iter, iter+1);
    in->maximumIterations = iter+1;
  }
}

void setSolutionTransfer(Input* in, SolutionTransfer* s)
{
  if (s)
  {
    in->solutionTransfer = s;
    in->ownsSolutionTransfer = false;
  }
  else
  {
    in->solutionTransfer = new AutoSolutionTransfer(in->mesh);
    in->ownsSolutionTransfer = true;
  }
}

static Input* configure(
    Mesh* m,
    SolutionTransfer* s)
{
  Input* in = new Input;
  in->mesh = m;
  setDefaultValues(in);
  setSolutionTransfer(in,s);
  return in;
}

const Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s,
    bool logInterpolation)
{
/* AutoSolutionTransfer object has to be created
   before the MetricSizeField, to avoid
   AutoSolutionTransfer taking ownership of
   the metric field, which has its own built-in
   solution transfer */
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f, logInterpolation);
  updateMaxIterBasedOnSize(m, in);
  return in;
}

const Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  updateMaxIterBasedOnSize(m, in);
  return in;
}

const Input* configure(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  updateMaxIterBasedOnSize(m, in);
  return in;
}

const Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, sizes, frames, logInterpolation);
  updateMaxIterBasedOnSize(m, in);
  return in;
}

const Input* configureUniformRefine(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = new UniformRefiner(m);
  in->maximumIterations = n;
  in->shouldRefineLayer = true;
  in->splitAllLayerEdges = true;
  return in;
}

const Input* configureMatching(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = makeAdvanced(configureUniformRefine(m,n,s));
  in->shouldHandleMatching = true;
  in->shouldFixShape = false;
  return in;
}

const Input* configureIdentity(Mesh* m, SizeField* f, SolutionTransfer* s)
{
  Input* in = configure(m,s);
  if (f)
  {
    in->sizeField = f;
    in->ownsSizeField = false;
  }
  else
  {
    in->sizeField = new IdentitySizeField(m);
    in->ownsSizeField = true;
  }
  in->maximumIterations = 0;
  in->shouldFixShape = false;
  in->shouldSnap = false;
  return in;
}

Input* makeAdvanced(const Input* in)
{
  Input* in2 = new Input(*in);
  delete in;
  return in2;
}

}
