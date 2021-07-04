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

Input::Input(const Input& in)
{
  this->ops = in.ops;
  this->mesh = in.mesh;
  this->sizeField = in.sizeField;
  this->ownsSizeField = in.ownsSizeField;
  this->solutionTransfer = in.solutionTransfer;
  this->ownsSolutionTransfer = in.ownsSolutionTransfer;
  this->shapeHandler = in.shapeHandler;
}

void Input::setDefaultValues()
{
  ownsSizeField = true;
  shapeHandler = 0;
  ops.maximumIterations = 3;
  ops.shouldCoarsen = true;
  ops.shouldSnap = mesh->canSnap();
  ops.shouldTransferParametric = mesh->canSnap();
  ops.shouldTransferToClosestPoint = false;
  ops.shouldHandleMatching = mesh->hasMatching();
  ops.shouldFixShape = true;
  ops.shouldForceAdaptation = false;
  ops.shouldPrintQuality = true;
  if (mesh->getDimension()==3)
  {
    ops.goodQuality = 0.027;
    ops.maximumEdgeRatio = 2.0;
  }
  else
  { PCU_ALWAYS_ASSERT(mesh->getDimension()==2);
    //old MA says .04, but that rarely kicks in.
    //.2 is strict, but at least quality goes up
    ops.goodQuality = 0.2;
    //this basically turns off short edge removal...
    ops.maximumEdgeRatio = 100.0;
    //2D mesh adapt performs better if forceAdapt is on
    ops.shouldForceAdaptation = true;
  }
  ops.shouldCheckQualityForDoubleSplits = false;
  ops.validQuality = 1e-10;
  ops.maximumImbalance = 1.10;
  ops.shouldRunPreZoltan = false;
  ops.shouldRunPreZoltanRib = false;
  ops.shouldRunPreParma = false;
  ops.shouldRunMidZoltan = false;
  ops.shouldRunMidParma = false;
  ops.shouldRunPostZoltan = false;
  ops.shouldRunPostZoltanRib = false;
  ops.shouldRunPostParma = false;
  ops.shouldTurnLayerToTets = false;
  ops.shouldCleanupLayer = false;
  ops.shouldRefineLayer = false;
  ops.shouldCoarsenLayer = false;
  ops.splitAllLayerEdges = false;
  ops.userDefinedLayerTagName = "";
}

void rejectInput(const char* str)
{
  if (PCU_Comm_Self() != 0)
    return;
  lion_eprint(1,"MeshAdapt input error:\n");
  lion_eprint(1,"%s\n",str);
  abort();
}

void rejectSetter(const char* str)
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
void Input::validateInput()
{
  if ( ! sizeField)
    rejectInput("no size field");
  if ( ! solutionTransfer)
    rejectInput("no solution transfer object");
  if (ops.maximumIterations < 0)
    rejectInput("negative maximum iteration count");
  if (ops.maximumIterations > 10)
    rejectInput("unusually high maximum iteration count");
  if (ops.shouldSnap
    &&( ! mesh->canSnap()))
    rejectInput("user requested snapping "
                "but the geometric model does not support it");
  if (ops.shouldTransferParametric
    &&( ! mesh->canSnap()))
    rejectInput("user requested parametric coordinate transfer "
                "but the geometric model does not support it");
  if (ops.shouldTransferToClosestPoint
    &&( ! mesh->canSnap()))
    rejectInput("user requested transfer to closest point on model"
                "but the geometric model does not support it");
  if (ops.shouldSnap && ( ! (ops.shouldTransferParametric ||
			     ops.shouldTransferToClosestPoint)))
    rejectInput("snapping requires parametric coordinate transfer or transfer to closest point");
  if ((mesh->hasMatching())
    &&( ! ops.shouldHandleMatching))
    rejectInput("the mesh has matching entities but matched support is off");
  if (ops.shouldHandleMatching
    && ops.shouldFixShape)
    rejectInput("user requested matched mesh handling and shape correction "
        "but shape correction does not support matching yet");
  if (ops.goodQuality < 0.0)
    rejectInput("negative desired element quality");
  if (ops.goodQuality > 1.0)
    rejectInput("desired element quality greater than one");
  if (ops.validQuality < 0.0)
    rejectInput("negative minimum element quality");
  if (ops.maximumImbalance < 1.0)
    rejectInput("maximum imbalance less than 1.0");
  if (ops.maximumEdgeRatio < 1.0)
    rejectInput("maximum tet edge ratio less than one");
  if (moreThanOneOptionIsTrue(
  	ops.shouldRunPreZoltan,
  	ops.shouldRunPreZoltanRib,
  	ops.shouldRunPreParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PreBalance options can be set to true!");
  if (moreThanOneOptionIsTrue(
  	ops.shouldRunPostZoltan,
  	ops.shouldRunPostZoltanRib,
  	ops.shouldRunPostParma))
    rejectInput("only one of Zoltan, ZoltanRib, and Parma PostBalance options can be set to true!");
  if (ops.shouldRunMidZoltan && ops.shouldRunMidParma)
    rejectInput("only one of Zoltan and Parma MidBalance options can be set to true!");
#ifndef PUMI_HAS_ZOLTAN
  if (ops.shouldRunPreZoltan ||
      ops.shouldRunPreZoltanRib ||
      ops.shouldRunMidZoltan)
    rejectInput("core is not compiled with Zoltan. Use a different balancer or compile core with ENABLE_ZOLTAN=ON!");
#endif
}

void Input::updateMaxIterBasedOnSize()
{
  double maxMetricLength = getMaximumEdgeLength(mesh, sizeField);
  int iter = std::ceil(std::log2(maxMetricLength));
  if (iter >= 10) {
    print("ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations,\n"
    	"           which is equal to or larger than the maximum of 10 allowed.\n"
    	"           Setting the number of iteration to 10!", iter);
    ops.maximumIterations = 10;
  }
  else {
    print("ma::configure:  Based on requested sizefield, MeshAdapt requires at least %d iterations.\n"
    	"           Setting the number of iteration to %d!", iter, iter+1);
    ops.maximumIterations = iter+1;
  }
}

void Input::init(
    Mesh* m,
    SolutionTransfer* s)
{
  /* Input* in = new Input; */
  mesh = m;
  setDefaultValues();
  if (s)
  {
    solutionTransfer = s;
    ownsSolutionTransfer = false;
  }
  else
  {
    solutionTransfer = new AutoSolutionTransfer(mesh);
    ownsSolutionTransfer = true;
  }
}

Input::Input(
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
  init(m,s);
  sizeField = makeSizeField(m, f, logInterpolation);
  updateMaxIterBasedOnSize();
}

Input::Input(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  init(m,s);
  sizeField = makeSizeField(m, f);
  updateMaxIterBasedOnSize();
}

Input::Input(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  init(m,s);
  if (f) {
    sizeField = makeSizeField(m, f);
    updateMaxIterBasedOnSize();
  }
}

Input::Input(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  init(m,s);
  sizeField = makeSizeField(m, sizes, frames, logInterpolation);
  updateMaxIterBasedOnSize();
}

Input::Input(Mesh* m, int n, SolutionTransfer* s)
{
  init(m,s);
  sizeField = new UniformRefiner(m);
  ops.maximumIterations = n;
  ops.shouldRefineLayer = true;
  ops.splitAllLayerEdges = true;
}

Input::Input(Mesh* m, SizeField* f, SolutionTransfer* s)
{
  init(m,s);
  if (f)
  {
    sizeField = f;
    ownsSizeField = false;
  }
  else
  {
    sizeField = new IdentitySizeField(m);
    ownsSizeField = true;
  }
  ops.maximumIterations = 0;
  ops.shouldFixShape = false;
  ops.shouldSnap = false;
}

Input* makeAdvanced(const Input* in)
{
  return const_cast<Input*>(in);
}

const Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new Input(m, f, s, logInterpolation);
}

const Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  return new Input(m, f, s);
}

const Input* configure(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  return new Input(m, f, s);
}

const Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new Input(m, sizes, frames, s, logInterpolation);
}

const Input* configureUniformRefine(Mesh* m, int n, SolutionTransfer* s)
{
  return new Input(m, n, s);
}

const Input* configureMatching(Mesh* m, int n, SolutionTransfer* s)
{
  const Input* in = configureUniformRefine(m,n,s);
  Input* inAdv = makeAdvanced(in);
  inAdv->shouldHandleMatching(true);
  inAdv->shouldFixShape(false);
  return (const Input*) inAdv;
}

const Input* configureIdentity(Mesh* m, SizeField* f, SolutionTransfer* s)
{
  return new Input(m, f, s);
}

}
