/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_INPUT_H
#define MA_INPUT_H

/** \file maInput.h
  \brief MeshAdapt user configuration
  \details Usually, users of MeshAdapt need to tweak its behavior slightly
           from the default configuration. For this purpose, the ma::Input
           object is the fully configurable way to interact with MeshAdapt.
           Users should call one of the configure functions in this file,
           then edit the options in the resulting ma::Input object,
           then give this object to the adapt function. */

#include "maMesh.h"
#include "maSize.h"
#include "maSolutionTransfer.h"

namespace ma {

/** \brief User configuration for a MeshAdapt run */
class Input
{
  public:
    ~Input();
    Mesh* mesh;
    SizeField* sizeField;
    bool ownsSizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSolutionTransfer;
/** \brief number of refine/coarsen iterations to run (default 3) */
    int maximumIterations;
/** \brief whether to snap new vertices to the model surface
    \details requires modeler support, see gmi_can_eval */
    bool shouldSnap;
/** \brief whether to transfer parametric coordinates
  \details requires modeler support, see gmi_reparam */
    bool shouldTransferParametric;
/** \brief whether to update matched entity info (limited support) */
    bool shouldHandleMatching;
/** \brief whether to run shape correction (default true) */
    bool shouldFixShape;
/** \brief whether to print the worst shape quality */
    bool shouldPrintQuality;
/** \brief minimum desired mean ratio cubed for simplex elements
   \details a different measure is used for curved elements */
    double goodQuality;
/** \brief minimum valid mean ratio cubed for simplex elements (default 1e-10)
   \details used to define inside-out tetrahedra.
   a different measure is used for curved elements */
    double validQuality;
/** \brief imbalance target for all load balancing tools (default 1.10) */
    double maximumImbalance;
/** \brief whether to run zoltan predictive load balancing (default false) */
    bool shouldRunPreZoltan;
/** \brief whether to run parma predictive load balancing (default false) */
    bool shouldRunPreParma;
/** \brief whether to run zoltan during adaptation (default false) */
    bool shouldRunMidZoltan;
/** \brief whether to run parma during adaptation (default false)*/
    bool shouldRunMidParma;
/** \brief whether to run zoltan after adapting (default false) */
    bool shouldRunPostZoltan;
/** \brief whether to run parma after adapting (default false) */
    bool shouldRunPostParma;
/** \brief the ratio between longest and shortest edges that differentiates a
   "short edge" element from a "large angle" element. */
    double maximumEdgeRatio;
/** \brief whether to tetrahedronize the boundary layer (default false)  */
    bool shouldTurnLayerToTets;
/** \brief whether to tetrahedronize abnormal pyramids (default false) */
    bool shouldCleanupLayer;
/** \brief whether to allow layer refinement (default false) */
    bool shouldRefineLayer;
/** \brief whether to allow layer coarsening (default false) */
    bool shouldCoarsenLayer;
/** \brief hack to enable boundary layer uniform refinement (do not touch!) */
    bool isUniform;
};

/** \brief generate a configuration based on an anisotropic function.
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on an isotropic function
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on anisotropic fields
 \param sizes a vector field of desired element sizes along the
              axes of the anisotropic frame
 \param frames a matrix field containing anisotropic frames
               for each vertex
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on an isotropic field
 \param size a scalar field of desired element size
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configure(
    Mesh* m,
    apf::Field* size,
    SolutionTransfer* s=0);

/** \brief generate a uniform refinement configuration */
Input* configureUniformRefine(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a matched uniform refinement configuration */
Input* configureMatching(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a no-op configuration */
Input* configureIdentity(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);

/** \brief used internally, but users can call this if they want */
void validateInput(Input* in);

}

#endif
