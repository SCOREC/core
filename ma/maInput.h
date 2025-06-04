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

class ShapeHandler;
class Adapt;

typedef ShapeHandler* (*ShapeHandlerFunction)(Adapt* a);

/** \brief User configuration for a MeshAdapt run */
class Input
{
  public:
    ~Input() {}
    Input() {} // default empty c-tor
    Input(const Input& in) = default; // copy c-tor
    Mesh* mesh;
    SizeField* sizeField;
    bool ownsSizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSolutionTransfer;
    ShapeHandlerFunction shapeHandler;
/** \brief number of refine/coarsen iterations to run (default 3)
    \details this value will be set to the minimum required iterations
    inside the call to ma::configure in cases where there is a size information.
    Users can override this by setting in->maximumIterations after the
    call to ma::configure and before the call to ma::adapt routine.*/
    int maximumIterations;
/** \brief whether to perform the collapse step */
    bool shouldCoarsen;
/** \brief whether to snap new vertices to the model surface
    \details requires modeler support, see gmi_can_eval */
    bool shouldSnap;
/** \brief whether to transfer parametric coordinates
    \details requires modeler support, see gmi_reparam */
    bool shouldTransferParametric;
/** \brief whether to transfer to the parametric coords of the closest point
    \details requires modeler support, see gmi_closest_point */
    bool shouldTransferToClosestPoint;
/** \brief whether to update matched entity info (limited support) */
    bool shouldHandleMatching;
/** \brief whether to run shape correction (default true) */
    bool shouldFixShape;
/** \brief whether to adapt if it makes local quality worse (default false) */
    bool shouldForceAdaptation;
/** \brief whether to print the worst shape quality */
    bool shouldPrintQuality;
/** \brief minimum desired mean ratio cubed for simplex elements
    \details a different measure is used for curved elements */
    double goodQuality;
/** \brief whether to check the quality of split elems in DoubleSplitsCollapse
    (default false) */
    double shouldCheckQualityForDoubleSplits;
/** \brief minimum valid mean ratio cubed for simplex elements (default 1e-10)
    \details used to define inside-out tetrahedra.
    a different measure is used for curved elements */
    double validQuality;
/** \brief imbalance target for all load balancing tools (default 1.10) */
    double maximumImbalance;
/** \brief whether to run zoltan predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreZoltan;
/** \brief whether to run zoltan predictive load balancing using RIB (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreZoltanRib;
    /**
     * \brief Whether to run METIS before adaptation (default false)
     *
     * If all pre-balance options are false, pre-balancing only occurs if the
     * estimated imbalance is greater than in->maximumImbalance.
     */
    bool shouldRunPreMetis;
/** \brief whether to run parma predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreParma;
/** \brief whether to run zoltan during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunMidZoltan;
    /**
     * \brief Whether to run METIS during adaptation (default false)
     *
     * If all mid-balance options are false, mid-balancing only occurs if the
     * estimated imbalance is greater than in->maximumImbalance.
     */
    bool shouldRunMidMetis;
/** \brief whether to run parma during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunMidParma;
/** \brief whether to run zoltan after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPostZoltan;
/** \brief whether to run zoltan RIB after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPostZoltanRib;
    /**
     * \brief Whether to run METIS after adaptation (default false)
     *
     * If all post-balance options are false, post-balancing only occurs if the
     * estimated imbalance is greater than in->maximumImbalance.
     */
    bool shouldRunPostMetis;
/** \brief whether to run parma after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
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
/** \brief set to true during UR to get splits in the normal direction */
    bool splitAllLayerEdges;
/** \brief the name of the (user defined) INT tag specifying the boundary
    layer elements. Use the value of 0 for non-layer elements and a non-zero value
    for layer elements. (default "") */
    const char* userDefinedLayerTagName;
/** \brief this a folder that debugging meshes will be written to, if provided! */
    const char* debugFolder;
};

/** \brief generate a configuration based on an anisotropic function.
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
const Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
/** \brief generate a configuration based on an isotropic function
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
const Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on anisotropic fields
 \param sizes a vector field of desired element sizes along the
              axes of the anisotropic frame
 \param frames a matrix field containing anisotropic frames
               for each vertex along the columns
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
const Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
/** \brief generate a configuration based on an isotropic field
 \param size a scalar field of desired element size
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
const Input* configure(
    Mesh* m,
    apf::Field* size,
    SolutionTransfer* s=0);

/** \brief generate a uniform refinement configuration */
const Input* configureUniformRefine(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a matched uniform refinement configuration */
const Input* configureMatching(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a no-op configuration */
const Input* configureIdentity(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);

/** \brief used internally, but users can call this if they want */
void validateInput(Input* in);

/** \brief creates a new non-constant Input for advanced users */
Input* makeAdvanced(const Input* in);

}

#endif
