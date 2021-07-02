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
/** \brief number of refine/coarsen iterations to run (default 3)
    \details this value will be set to the minimum required iterations
    inside the call to ma::configure in cases where there is a size information.
    Users can override this by setting in->maximumIterations after the
    call to ma::configure and before the call to ma::adapt routine.*/
    int _maximumIterations;
/** \brief whether to perform the collapse step */
    bool _shouldCoarsen;
/** \brief whether to snap new vertices to the model surface
    \details requires modeler support, see gmi_can_eval */
    bool _shouldSnap;
/** \brief whether to transfer parametric coordinates
    \details requires modeler support, see gmi_reparam */
    bool _shouldTransferParametric;
/** \brief whether to transfer to the parametric coords of the closest point
    \details requires modeler support, see gmi_closest_point */
    bool _shouldTransferToClosestPoint;
/** \brief whether to update matched entity info (limited support) */
    bool _shouldHandleMatching;
/** \brief whether to run shape correction (default true) */
    bool _shouldFixShape;
/** \brief whether to adapt if it makes local quality worse (default false) */
    bool _shouldForceAdaptation;
/** \brief whether to print the worst shape quality */
    bool _shouldPrintQuality;
/** \brief minimum desired mean ratio cubed for simplex elements
    \details a different measure is used for curved elements */
    double _goodQuality;
/** \brief whether to check the quality of split elems in DoubleSplitsCollapse
    (default false) */
    double _shouldCheckQualityForDoubleSplits;
/** \brief minimum valid mean ratio cubed for simplex elements (default 1e-10)
    \details used to define inside-out tetrahedra.
    a different measure is used for curved elements */
    double _validQuality;
/** \brief imbalance target for all load balancing tools (default 1.10) */
    double _maximumImbalance;
/** \brief whether to run zoltan predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPreZoltan;
/** \brief whether to run zoltan predictive load balancing using RIB (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPreZoltanRib;
/** \brief whether to run parma predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPreParma;
/** \brief whether to run zoltan during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunMidZoltan;
/** \brief whether to run parma during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunMidParma;
/** \brief whether to run zoltan after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPostZoltan;
/** \brief whether to run zoltan RIB after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPostZoltanRib;
/** \brief whether to run parma after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool _shouldRunPostParma;
/** \brief the ratio between longest and shortest edges that differentiates a
    "short edge" element from a "large angle" element. */
    double _maximumEdgeRatio;
/** \brief whether to tetrahedronize the boundary layer (default false)  */
    bool _shouldTurnLayerToTets;
/** \brief whether to tetrahedronize abnormal pyramids (default false) */
    bool _shouldCleanupLayer;
/** \brief whether to allow layer refinement (default false) */
    bool _shouldRefineLayer;
/** \brief whether to allow layer coarsening (default false) */
    bool _shouldCoarsenLayer;
/** \brief set to true during UR to get splits in the normal direction */
    bool _splitAllLayerEdges;
/** \brief the name of the (user defined) INT tag specifying the boundary
    layer elements. Use the value of 0 for non-layer elements and a non-zero value
    for layer elements. (default "") */
    const char* _userDefinedLayerTagName;
/** \brief this a folder that debugging meshes will be written to, if provided! */
    const char* _debugFolder;
  public:
    ~Input();
    Mesh* mesh;
    SizeField* sizeField;
    bool ownsSizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSolutionTransfer;
    ShapeHandlerFunction shapeHandler;
    int maximumIterations() { return _maximumIterations; }
    bool shouldFixShape() { return _shouldFixShape; }
    bool shouldTransferParametric() { return _shouldTransferParametric; }
    bool shouldTransferToClosestPoint() { return _shouldTransferToClosestPoint; }
    bool shouldRefineLayer() { return _shouldRefineLayer; }
    bool shouldCoarsen() { return _shouldCoarsen; }
    bool shouldCoarsenLayer() { return _shouldCoarsenLayer; }
    bool shouldRunPreZoltan() { return _shouldRunPreZoltan; }
    bool shouldRunPreZoltanRib() { return _shouldRunPreZoltanRib; }
    bool shouldRunPreParma() { return _shouldRunPreParma; }
    bool shouldRunMidZoltan() { return _shouldRunMidZoltan; }
    bool shouldRunMidParma() { return _shouldRunMidParma; }
    bool shouldRunPostZoltan() { return _shouldRunPostZoltan; }
    bool shouldRunPostZoltanRib() { return _shouldRunPostZoltanRib; }
    bool shouldRunPostParma() { return _shouldRunPostParma; }
    int maximumImbalance() { return _maximumImbalance; }
    bool shouldSnap() { return _shouldSnap; }
    bool shouldForceAdaptation() { return _shouldForceAdaptation; }
    void setShouldForceAdaptation(bool b) { _shouldForceAdaptation = b; }
    double goodQuality() { return _goodQuality; }
    double validQuality() { return _validQuality; }
    bool shouldHandleMatching() { return _shouldHandleMatching; }
    bool splitAllLayerEdges() { return _splitAllLayerEdges; }
    bool shouldCheckQualityForDoubleSplits() { return _shouldCheckQualityForDoubleSplits; }
    bool shouldCleanupLayer() { return _shouldCleanupLayer; }
    bool shouldPrintQuality() { return _shouldPrintQuality; }
    double maximumEdgeRatio() { return _maximumEdgeRatio; }
    bool shouldTurnLayerToTets() { return _shouldTurnLayerToTets; }
    const char* userDefinedLayerTagName() { return _userDefinedLayerTagName; }
    const char* debugFolder() { return _debugFolder; }

    // friend classes are defined here so they can change the private members above
    friend Input* configure(
	Mesh*, AnisotropicFunction*, SolutionTransfer*, bool);
    friend Input* configure(
	Mesh*, IsotropicFunction*, SolutionTransfer*);
    friend Input* configure(
	Mesh*, apf::Field*, apf::Field*, SolutionTransfer*, bool);
    friend Input* configure(Mesh*, apf::Field*, SolutionTransfer*);
    friend Input* configureUniformRefine(Mesh*, int, SolutionTransfer*);
    friend Input* configureMatching(Mesh*, int, SolutionTransfer*);
    friend Input* configureIdentity(Mesh*, SizeField*, SolutionTransfer*);
    friend void setDefaultValues(Input* in);
    friend void validateInput(Input*);
    friend void updateMaxIterBasedOnSize(Input*);
    /* friend void setSolutionTransfer(Input*, SolutionTransfer*); */
    friend void setAdvancedInputOptions(Input*);
};

/** \brief generate a configuration based on an anisotropic function.
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
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
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
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

/** \brief gives advanced users the option to manually set input adapt options
 \param in pointer to input object. This function is defined as a friend of Input
           class but is never implemented in the code. The users can implement it
           in their application code and override the default adapt options by
           setting them using statements like "in->shouldRunPreZoltan = true".
           see core/test/torus_ma_test.cc for and example of how this is done.  */
void setAdvancedInputOptions(Input* in);

}

#endif
