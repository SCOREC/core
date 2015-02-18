/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maInput.h"
#include <apfShape.h>
#include <cstdio>

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
  in->maximumIterations = 3;
  in->shouldSnap = in->mesh->canSnap();
  in->shouldTransferParametric = in->mesh->canSnap();
  in->shouldHandleMatching = in->mesh->hasMatching();
  in->shouldFixShape = true;
  in->shouldPrintQuality = true;
  if (in->mesh->getDimension()==3)
  {
    in->goodQuality = 0.027;
    in->maximumEdgeRatio = 2.0;
  }
  else
  { assert(in->mesh->getDimension()==2);
    //old MA says .04, but that rarely kicks in.
    //.2 is strict, but at least quality goes up
    in->goodQuality = 0.2;
    //this basically turns off short edge removal...
    in->maximumEdgeRatio = 100.0;
  }
  in->validQuality = 1e-10;
  in->maximumImbalance = 1.10;
  in->shouldRunPreZoltan = false;
  in->shouldRunPreParma = false;
  in->shouldRunMidZoltan = false;
  in->shouldRunMidParma = false;
  in->shouldRunPostZoltan = false;
  in->shouldRunPostParma = false;
  in->shouldTurnLayerToTets = false;
  in->shouldCleanupLayer = false;
  in->shouldRefineLayer = false;
  in->shouldCoarsenLayer = false;
  in->isUniform = false;
}

void rejectInput(const char* str)
{
  fprintf(stderr,"MeshAdapt input error:\n");
  fprintf(stderr,"%s\n",str);
  abort();
}

void validateInput(Input* in)
{
  if ( ! in->sizeField)
    rejectInput("no size field");
  if ( ! in->solutionTransfer)
    rejectInput("no solution transfer object");
  if (in->maximumIterations < 0)
    rejectInput("negative maximum iteration count");
  if (in->maximumIterations > 10)
    rejectInput("unusually high maximum iteration count");
  if (in->shouldSnap
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested snapping "
                "but the geometric model does not support it");
  if (in->shouldTransferParametric
    &&( ! in->mesh->canSnap()))
    rejectInput("user requested parametric coordinate transfer "
                "but the geometric model does not support it");
  if (in->shouldSnap && ( ! in->shouldTransferParametric))
    rejectInput("snapping requires parametric coordinate transfer");
  if ((in->mesh->hasMatching())
    &&( ! in->shouldHandleMatching))
    rejectInput("the mesh has matching entities but matched support is off");
  if (in->shouldHandleMatching
    && in->shouldFixShape)
    rejectInput("user requested matched mesh handling and shape correction "
        "but shape correction does not support matching yet");
  if (in->goodQuality < 0.0)
    rejectInput("negative desired element quality");
  if (in->goodQuality > 1.0)
    rejectInput("desired element quality greater than one");
  if (in->validQuality < 0.0)
    rejectInput("negative minimum element quality");
  if (in->maximumImbalance < 1.0)
    rejectInput("maximum imbalance less than 1.0");
  if (in->maximumEdgeRatio < 1.0)
    rejectInput("maximum tet edge ratio less than one");
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
    SolutionTransfer* s)
{
/* AutoSolutionTransfer object has to be created
   before the MetricSizeField, to avoid
   AutoSolutionTransfer taking ownership of
   the metric field, which has its own built-in
   solution transfer */
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  return in;
}

Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  return in;
}

Input* configure(
    Mesh* m,
    apf::Field* f,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, f);
  return in;
}

Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = makeSizeField(m, sizes, frames);
  return in;
}

Input* configureUniformRefine(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configure(m,s);
  in->sizeField = new UniformRefiner(m);
  in->maximumIterations = n;
  in->shouldRefineLayer = true;
  in->isUniform = true;
  return in;
}

Input* configureMatching(Mesh* m, int n, SolutionTransfer* s)
{
  Input* in = configureUniformRefine(m,n,s);
  in->shouldHandleMatching = true;
  in->shouldFixShape = false;
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
  in->maximumIterations = 0;
  in->shouldFixShape = false;
  in->shouldSnap = false;
  return in;
}

}
