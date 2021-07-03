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

int Input::maximumIterations() { return ops.maximumIterations; }
void Input::maximumIterations(int i)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.maximumIterations = i;
}

bool Input::shouldCoarsen() { return ops.shouldCoarsen; }
void Input::shouldCoarsen(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldCoarsen = b;
}

bool Input::shouldSnap() { return ops.shouldSnap; }
void Input::shouldSnap(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldSnap = b;
}

bool Input::shouldTransferParametric() { return ops.shouldTransferParametric; }
void Input::shouldTransferParametric(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldTransferParametric = b;
}

bool Input::shouldTransferToClosestPoint() { return ops.shouldTransferToClosestPoint; }
void Input::shouldTransferToClosestPoint(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldTransferToClosestPoint = b;
}

bool Input::shouldHandleMatching() { return ops.shouldHandleMatching; }
void Input::shouldHandleMatching(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldHandleMatching = b;
}

bool Input::shouldFixShape() { return ops.shouldFixShape; }
void Input::shouldFixShape(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldFixShape = b;
}

bool Input::shouldForceAdaptation() { return ops.shouldForceAdaptation; }
void Input::shouldForceAdaptation(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldForceAdaptation = b;
}
/* void Input::shouldForceAdaptation(bool b) { ops.shouldFixShape = b; } */

bool Input::shouldPrintQuality() { return ops.shouldPrintQuality; }
void Input::shouldPrintQuality(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldPrintQuality = b;
}

double Input::goodQuality() { return ops.goodQuality; }
void Input::goodQuality(double d)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.goodQuality = d;
}

bool Input::shouldCheckQualityForDoubleSplits() { return ops.shouldCheckQualityForDoubleSplits; }
void Input::shouldCheckQualityForDoubleSplits(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldCheckQualityForDoubleSplits = b;
}

double Input::validQuality() { return ops.validQuality; }
void Input::validQuality(double d)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.validQuality = d;
}

double Input::maximumImbalance() { return ops.maximumImbalance; }
void Input::maximumImbalance(double d)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.maximumImbalance = d;
}

bool Input::shouldRunPreZoltan() { return ops.shouldRunPreZoltan; }
void Input::shouldRunPreZoltan(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPreZoltan = b;
}

bool Input::shouldRunPreZoltanRib() { return ops.shouldRunPreZoltanRib; }
void Input::shouldRunPreZoltanRib(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPreZoltanRib = b;
}

bool Input::shouldRunPreParma() { return ops.shouldRunPreParma; }
void Input::shouldRunPreParma(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPreParma = b;
}


bool Input::shouldRunMidZoltan() { return ops.shouldRunMidZoltan; }
void Input::shouldRunMidZoltan(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunMidZoltan = b;
}

bool Input::shouldRunMidParma() { return ops.shouldRunMidParma; }
void Input::shouldRunMidParma(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunMidParma = b;
}

bool Input::shouldRunPostZoltan() { return ops.shouldRunPostZoltan; }
void Input::shouldRunPostZoltan(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPostZoltan = b;
}

bool Input::shouldRunPostZoltanRib() { return ops.shouldRunPostZoltanRib; }
void Input::shouldRunPostZoltanRib(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPostZoltanRib = b;
}

bool Input::shouldRunPostParma() { return ops.shouldRunPostParma; }
void Input::shouldRunPostParma(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRunPostParma = b;
}

double Input::maximumEdgeRatio() { return ops.maximumEdgeRatio; }
void Input::maximumEdgeRatio(double d)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.maximumEdgeRatio = d;
}

bool Input::shouldTurnLayerToTets() { return ops.shouldTurnLayerToTets; }
void Input::shouldTurnLayerToTets(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldTurnLayerToTets = b;
}

bool Input::shouldCleanupLayer() { return ops.shouldCleanupLayer; }
void Input::shouldCleanupLayer(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldCleanupLayer = b;
}

bool Input::shouldRefineLayer() { return ops.shouldRefineLayer; }
void Input::shouldRefineLayer(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldRefineLayer = b;
}

bool Input::shouldCoarsenLayer() { return ops.shouldCoarsenLayer; }
void Input::shouldCoarsenLayer(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.shouldCoarsenLayer = b;
}

bool Input::splitAllLayerEdges() { return ops.splitAllLayerEdges; }
void Input::splitAllLayerEdges(bool b)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.splitAllLayerEdges = b;
}

const char* Input::userDefinedLayerTagName() { return ops.userDefinedLayerTagName; }
void Input::userDefinedLayerTagName(const char* c)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.userDefinedLayerTagName = c;
}

const char* Input::debugFolder() { return ops.debugFolder; }
void Input::debugFolder(const char* c)
{
  if (!modifiable())
    rejectSetter("Cannot change mesh adapt options for basic input!");
  else
    ops.debugFolder = c;
}



Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new InputBasic(m, f, s, logInterpolation);
}
Input* configureAdvanced(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new InputAdvanced(m, f, s, logInterpolation);
}

Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  return new InputBasic(m, f, s);
}
Input* configureAdvanced(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  return new InputAdvanced(m, f, s);
}

Input* configure(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  return new InputBasic(m, f, s);
}
Input* configureAdvanced(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  return new InputAdvanced(m, f, s);
}

Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new InputBasic(m, sizes, frames, s, logInterpolation);
}
Input* configureAdvanced(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s,
    bool logInterpolation)
{
  return new InputAdvanced(m, sizes, frames, s, logInterpolation);
}

Input* configureUniformRefine(Mesh* m, int n, SolutionTransfer* s)
{
  return new InputBasic(m, n, s);
}
Input* configureUniformRefineAdvanced(Mesh* m, int n, SolutionTransfer* s)
{
  return new InputAdvanced(m, n, s);
}

Input* configureMatching(Mesh* m, int n, SolutionTransfer* s)
{
  Input* inAdv = configureUniformRefineAdvanced(m,n,s);
  inAdv->shouldHandleMatching(true);
  inAdv->shouldFixShape(false);
  // make an un-modifiable version
  Input* in = new InputBasic(inAdv);
  delete inAdv;
  return in;
}
Input* configureMatchingAdvanced(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configureUniformRefineAdvanced(m,n,s);
  in->shouldHandleMatching(true);
  in->shouldFixShape(false);
  return in;
}


Input* configureIdentity(Mesh* m, SizeField* f, SolutionTransfer* s)
{
  return new InputBasic(m, f, s);
}
Input* configureIdentityAdvanced(Mesh* m, SizeField* f, SolutionTransfer* s)
{
  return new InputAdvanced(m, f, s);
}

}
