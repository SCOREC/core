/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_INPUT_H
#define MA_INPUT_H

#include "maMesh.h"
#include "maSize.h"
#include "maSolutionTransfer.h"

namespace ma {

class Input
{
  public:
    ~Input();
    Mesh* mesh;
    SizeField* sizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSizeField;
/* whether or not we are resposible for deleting solutionTrasfer */
    bool ownsSolutionTransfer;
/* number of refine/coarsen iterations to run */
    int maximumIterations;
/* whether to snap new vertices to the model surface (requires modeler support) */
    bool shouldSnap;
/* whether to transfer parametric coordinates (requires modeler support) */
    bool shouldTransferParametric;
/* whether to update matched entity info (limited support) */
    bool shouldHandleMatching;
/* whether to run shape correction */
    bool shouldFixShape;
/* minimum desired mean ratio cubed for simplex elements
   (a different measure is used for curved elements) */
    double goodQuality;
/* minimum valid mean ratio cubed for simplex elements
   (a different measure is used for curved elements),
   used to define inside-out tetrahedra */
    double validQuality;
/* imbalance target for all load balancing tools */
    double maximumImbalance;
/* whether or not to run zoltan predictive load balancing before adapting */
    bool shouldRunPreZoltan;
/* whether or not to run parma predictive load balancing before adapting */
    bool shouldRunPreParma;
    bool shouldRunPreDiffusion;
/* whether or not to run zoltan load balancing during adaptation */
    bool shouldRunMidZoltan;
/* whether or not to run parma iterative load balancing during adaptation */
    bool shouldRunMidParma;
    bool shouldRunMidDiffusion;
/* whether or not to run zoltan load balancing after adapting */
    bool shouldRunPostZoltan;
/* whether or not to run parma cleanup load balancing after adapting */
    bool shouldRunPostParma;
    bool shouldRunPostDiffusion;
/* maximum iterations for parma's diffuser */
    int diffuseIterations;
/* the ratio between longest and shortest edges that differentiates a
   "short edge" element from a "large angle" element. */
    double maximumEdgeRatio;
/* whether to tetrahedronize the boundary layer */
    bool shouldTurnLayerToTets;
/* whether to tetrahedronize abnormal pyramids */
    bool shouldCleanupLayer;
/* whether to allow layer refinement */
    bool shouldRefineLayer;
/* hack to enable boundary layer uniform refinement. do not touch */
    bool isUniform;
};

Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0);
Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s=0);
Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s=0);
Input* configureUniformRefine(Mesh* m, int n=1, SolutionTransfer* s=0);
Input* configureMatching(Mesh* m, int n=1, SolutionTransfer* s=0);
Input* configureIdentity(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);

void validateInput(Input* in);

}

#endif
