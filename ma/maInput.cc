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
#include <PCU.h>
#include <pcu_util.h>
#include <cstdlib>

namespace ma {

Input::~Input()
{
  if (ownsSizeField)
    delete sizeField;
  if (ownsSolutionTransfer)
    delete solutionTransfer;
}

void setDefaultValues(Input* in)
{
  in->ownsSizeField = true;
  in->_maximumIterations = 3;
  in->_shouldCoarsen = true;
  in->_shouldSnap = in->mesh->canSnap();
  in->_shouldTransferParametric = in->mesh->canSnap();
  in->_shouldTransferToClosestPoint = false;
  in->_shouldHandleMatching = in->mesh->hasMatching();
  in->_shouldFixShape = true;
  in->_shouldForceAdaptation = false;
  in->_shouldPrintQuality = true;
  if (in->mesh->getDimension()==3)
  {
    in->_goodQuality = 0.027;
    in->_maximumEdgeRatio = 2.0;
  }
  else
  { PCU_ALWAYS_ASSERT(in->mesh->getDimension()==2);
    //old MA says .04, but that rarely kicks in.
    //.2 is strict, but at least quality goes up
    in->_goodQuality = 0.2;
    //this basically turns off short edge removal...
    in->_maximumEdgeRatio = 100.0;
    //2D mesh adapt performs better if forceAdapt is on
    in->_shouldForceAdaptation = true;
  }
  in->_shouldCheckQualityForDoubleSplits = false;
  in->_validQuality = 1e-10;
  in->_maximumImbalance = 1.10;
  in->_shouldRunPreZoltan = false;
  in->_shouldRunPreZoltanRib = false;
  in->_shouldRunPreParma = false;
  in->_shouldRunMidZoltan = false;
  in->_shouldRunMidParma = false;
  in->_shouldRunPostZoltan = false;
  in->_shouldRunPostZoltanRib = false;
  in->_shouldRunPostParma = false;
  in->_shouldTurnLayerToTets = false;
  in->_shouldCleanupLayer = false;
  in->_shouldRefineLayer = false;
  in->_shouldCoarsenLayer = false;
  in->_splitAllLayerEdges = false;
  in->_userDefinedLayerTagName = "";
  in->shapeHandler = 0;
}

void rejectInput(const char* str)
{
  if (PCU_Comm_Self() != 0)
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
    rejectInput("no size field");
  if ( ! in->solutionTransfer)
    rejectInput("no solution transfer object");
  if (in->_maximumIterations < 0)
    rejectInput("negative maximum iteration count");
  if (in->_maximumIterations > 10)
    rejectInput("unusually high maximum iteration count");
  if (in->_shouldSnap
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested snapping "
                "but the geometric model does not support it");
  if (in->_shouldTransferParametric
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested parametric coordinate transfer "
                "but the geometric model does not support it");
  if (in->_shouldTransferToClosestPoint
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested transfer to closest point on model"
                "but the geometric model does not support it");
  if (in->_shouldSnap && ( ! (in->_shouldTransferParametric ||
			     in->_shouldTransferToClosestPoint)))
    rejectInput("snapping requires parametric coordinate transfer or transfer to closest point");
  if ((in->mesh->hasMatching())
    &&( ! in->_shouldHandleMatching))
    rejectInput("the mesh has matching entities but matched support is off");
  if (in->_shouldHandleMatching
    && in->_shouldFixShape)
    rejectInput("user requested matched mesh handling and shape correction "
        "but shape correction does not support matching yet");
  if (in->_goodQuality < 0.0)
    rejectInput("negative desired element quality");
  if (in->_goodQuality > 1.0)
    rejectInput("desired element quality greater than one");
  if (in->_validQuality < 0.0)
    rejectInput("negative minimum element quality");
  if (in->_maximumImbalance < 1.0)
    rejectInput("maximum imbalance less than 1.0");
  if (in->_maximumEdgeRatio < 1.0)
    rejectInput("maximum tet edge ratio less than one");
  if (moreThanOneOptionIsTrue(
  	in->_shouldRunPreZoltan,
  	in->_shouldRunPreZoltanRib,
  	in->_shouldRunPreParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PreBalance options can be set to true!");
  if (moreThanOneOptionIsTrue(
  	in->_shouldRunPostZoltan,
  	in->_shouldRunPostZoltanRib,
  	in->_shouldRunPostParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PostBalance options can be set to true!");
  if (in->_shouldRunMidZoltan && in->_shouldRunMidParma)
    rejectInput("only one of Zoltan and Parma MidBalance options can be set to true!");
#ifndef PUMI_HAS_ZOLTAN
  if (in->_shouldRunPreZoltan ||
      in->_shouldRunPreZoltanRib ||
      in->_shouldRunMidZoltan)
    rejectInput("core is not compiled with Zoltan. Use a different balancer or compile core with ENABLE_ZOLTAN=ON!");
#endif
}

void updateMaxIterBasedOnSize(Input* in)
{
  Mesh* m = in->mesh;
  double maxMetricLength = getMaximumEdgeLength(m, in->sizeField);
  int iter = std::ceil(std::log2(maxMetricLength));
  if (iter >= 10) {
    print("ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations,\n"
    	"           which is equal to or larger than the maximum of 10 allowed.\n"
    	"           Setting the number of iteration to 10!", iter);
    in->_maximumIterations = 10;
  }
  else {
    print("ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations.\n"
    	"           Setting the number of iteration to %d!", iter, iter+1);
    in->_maximumIterations = iter+1;
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

Input* configure(
    Mesh* m,
    SolutionTransfer* s)
{
  Input* in = new Input;
  in->mesh = m;
  setDefaultValues(in);
  setSolutionTransfer(in,s);
  return in;
}

Input* configure(
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
  updateMaxIterBasedOnSize(in);
  return in;
}

Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  updateMaxIterBasedOnSize(in);
  return in;
}

Input* configure(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  if (f) {
    in->sizeField = makeSizeField(m, f);
    updateMaxIterBasedOnSize(in);
  }
  return in;
}

Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, sizes, frames, logInterpolation);
  updateMaxIterBasedOnSize(in);
  return in;
}

Input* configureUniformRefine(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = new UniformRefiner(m);
  in->_maximumIterations = n;
  in->_shouldRefineLayer = true;
  in->_splitAllLayerEdges = true;
  return in;
}

Input* configureMatching(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configureUniformRefine(m,n,s);
  in->_shouldHandleMatching = true;
  in->_shouldFixShape = false;
  return in;
}

Input* configureIdentity(Mesh* m, SizeField* f, SolutionTransfer* s)
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
  in->_maximumIterations = 0;
  in->_shouldFixShape = false;
  in->_shouldSnap = false;
  return in;
}

}
