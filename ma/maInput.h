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
  \details Most MeshAdapt options are set internally via the configure APIs.
           Occasionally (advanced) users of MeshAdapt need to tweak its behavior
           slightly from the default configuration. Such users should see
	   ma::makeAdvanced documentation for more information. */

#include "maMesh.h"
#include "maSize.h"
#include "maSolutionTransfer.h"

namespace ma {

class ShapeHandler;
class Adapt;

typedef ShapeHandler* (*ShapeHandlerFunction)(Adapt* a);


/** \brief Configuration options for a MeshAdapt run */
struct InputOptions
{
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
/** \brief whether to run parma predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreParma;
/** \brief whether to run zoltan during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunMidZoltan;
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


/** \brief User configuration for a MeshAdapt run */
class Input
{
  public:
    virtual ~Input();
    /** \brief generate an Input based on an anisotropic function.
    \param s if non-zero, use that to transfer all fields. otherwise,
	      transfer any associated fields with default algorithms
    \param logInterpolation if true uses logarithmic interpolation */
    Input(
      Mesh* m,
      AnisotropicFunction* f,
      SolutionTransfer* s=0,
      bool logInterpolation=true);
    /** \brief generate an Input based on an isotropic function
    \param s if non-zero, use that to transfer all fields. otherwise,
	      transfer any associated fields with default algorithms */
    Input(
      Mesh* m,
      IsotropicFunction* f,
      SolutionTransfer* s=0);
    /** \brief generate an Input based on anisotropic fields
    \param sizes a vector field of desired element sizes along the
		  axes of the anisotropic frame
    \param frames a matrix field containing anisotropic frames
		  for each vertex
    \param s if non-zero, use that to transfer all fields. otherwise,
	      transfer any associated fields with default algorithms
    \param logInterpolation if true uses logarithmic interpolation */
    Input(
      Mesh* m,
      apf::Field* sizes,
      apf::Field* frames,
      SolutionTransfer* s=0,
      bool logInterpolation=true);
    /** \brief generate an Input based on an isotropic field
    \param size a scalar field of desired element size
    \param s if non-zero, use that to transfer all fields. otherwise,
	      transfer any associated fields with default algorithms */
    Input(
      Mesh* m,
      apf::Field* size,
      SolutionTransfer* s=0);
    /** \brief generate a uniform refinement Input */
    Input(Mesh* m, int n=1, SolutionTransfer* s=0);
    /** \brief generate a no-op Input */
    Input(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);
    /** \brief getter for maximumIterations
        \details see InputOptions */
    int maximumIterations() { return ops.maximumIterations; }
    /** \brief setter for maximumIterations */
    void maximumIterations(int i) { ops.maximumIterations = i; }
    /** \brief getter for shouldCoarsen
        \details see InputOptions */
    bool shouldCoarsen() { return ops.shouldCoarsen; }
    /** \brief setter for shouldCoarsen */
    void shouldCoarsen(bool b) { ops.shouldCoarsen = b; }
    /** \brief getter for shouldSnap
        \details see InputOptions */
    bool shouldSnap() { return ops.shouldSnap; }
    /** \brief setter for shouldSnap */
    void shouldSnap(bool b) { ops.shouldSnap = b; }
    /** \brief getter for shouldTransferParametric
        \details see InputOptions */
    bool shouldTransferParametric() { return ops.shouldTransferParametric; }
    /** \brief setter for shouldTransferParametric */
    void shouldTransferParametric(bool b) { ops.shouldTransferParametric = b; }
    /** \brief getter for shouldTransferToClosestPoint
        \details see InputOptions */
    bool shouldTransferToClosestPoint() { return ops.shouldTransferToClosestPoint; }
    /** \brief setter for shouldTransferToClosestPoint */
    void shouldTransferToClosestPoint(bool b) { ops.shouldTransferToClosestPoint = b; }
    /** \brief getter for shouldHandleMatching
        \details see InputOptions */
    bool shouldHandleMatching() { return ops.shouldHandleMatching; }
    /** \brief setter for shouldHandleMatching */
    void shouldHandleMatching(bool b) { ops.shouldHandleMatching = b; }
    /** \brief getter for shouldFixShape
        \details see InputOptions */
    bool shouldFixShape() { return ops.shouldFixShape; }
    /** \brief setter for shouldFixShape */
    void shouldFixShape(bool b) { ops.shouldFixShape = b; }
    /** \brief getter for shouldForceAdaptation
        \details see InputOptions */
    bool shouldForceAdaptation() { return ops.shouldForceAdaptation; }
    /** \brief setter for shouldForceAdaptation */
    void shouldForceAdaptation(bool b) { ops.shouldForceAdaptation = b; }
    /** \brief getter for shouldPrintQuality
        \details see InputOptions */
    bool shouldPrintQuality() { return ops.shouldPrintQuality; }
    /** \brief setter for shouldPrintQuality */
    void shouldPrintQuality(bool b) { ops.shouldPrintQuality = b; }
    /** \brief getter for goodQuality
        \details see InputOptions */
    double goodQuality() { return ops.goodQuality; }
    /** \brief setter for goodQuality */
    void goodQuality(double d) { ops.goodQuality = d; }
    /** \brief getter for shouldCheckQualityForDoubleSplits
        \details see InputOptions */
    bool shouldCheckQualityForDoubleSplits() { return ops.shouldCheckQualityForDoubleSplits; }
    /** \brief setter for shouldCheckQualityForDoubleSplits */
    void shouldCheckQualityForDoubleSplits(bool b) { ops.shouldCheckQualityForDoubleSplits = b; }
    /** \brief getter for validQuality
        \details see InputOptions */
    double validQuality() { return ops.validQuality; }
    /** \brief setter for validQuality */
    void validQuality(double d) { ops.validQuality = d; }
    /** \brief getter for maximumImbalance
        \details see InputOptions */
    double maximumImbalance() { return ops.maximumImbalance; }
    /** \brief setter for maximumImbalance */
    void maximumImbalance(double d) { ops.maximumImbalance = d; }
    /** \brief getter for shouldRunPreZoltan
        \details see InputOptions */
    bool shouldRunPreZoltan() { return ops.shouldRunPreZoltan; }
    /** \brief setter for shouldRunPreZoltan */
    void shouldRunPreZoltan(bool b) { ops.shouldRunPreZoltan = b; }
    /** \brief getter for shouldRunPostZoltanRib
        \details see InputOptions */
    bool shouldRunPreZoltanRib() { return ops.shouldRunPreZoltanRib; }
    /** \brief setter for shouldRunPreZoltanRib */
    void shouldRunPreZoltanRib(bool b) { ops.shouldRunPreZoltanRib = b; }
    /** \brief getter for shouldRunPreParma
        \details see InputOptions */
    bool shouldRunPreParma() { return ops.shouldRunPreParma; }
    /** \brief setter for shouldRunPreParma */
    void shouldRunPreParma(bool b) { ops.shouldRunPreParma = b; }
    /** \brief getter for shouldRunMidZoltan
        \details see InputOptions */
    bool shouldRunMidZoltan() { return ops.shouldRunMidZoltan; }
    /** \brief setter for shouldRunMidZoltan */
    void shouldRunMidZoltan(bool b) { ops.shouldRunMidZoltan = b; }
    /** \brief getter for shouldRunMidParma
        \details see InputOptions */
    bool shouldRunMidParma() { return ops.shouldRunMidParma; }
    /** \brief setter for shouldRunMidParma */
    void shouldRunMidParma(bool b) { ops.shouldRunMidParma = b; }
    /** \brief getter for shouldRunPostZoltan
        \details see InputOptions */
    bool shouldRunPostZoltan() { return ops.shouldRunPostZoltan; }
    /** \brief setter for shouldRunPostZoltan */
    void shouldRunPostZoltan(bool b) { ops.shouldRunPostZoltan = b; }
    /** \brief getter for shouldRunPostZoltanRib
        \details see InputOptions */
    bool shouldRunPostZoltanRib() { return ops.shouldRunPostZoltanRib; }
    /** \brief setter for shouldRunPostZoltanRib */
    void shouldRunPostZoltanRib(bool b) { ops.shouldRunPostZoltanRib = b; }
    /** \brief getter for shouldRunPostParma
        \details see InputOptions */
    bool shouldRunPostParma() { return ops.shouldRunPostParma; }
    /** \brief setter for shouldRunPostParma */
    void shouldRunPostParma(bool b) { ops.shouldRunPostParma = b; }
    /** \brief getter for maximumEdgeRatio
        \details see InputOptions */
    double maximumEdgeRatio() { return ops.maximumEdgeRatio; }
    /** \brief setter for maximumEdgeRatio */
    void maximumEdgeRatio(double d) { ops.maximumEdgeRatio = d; }
    /** \brief getter for shouldTurnLayerToTets
        \details see InputOptions */
    bool shouldTurnLayerToTets() { return ops.shouldTurnLayerToTets; }
    /** \brief setter for shouldTurnLayerToTets */
    void shouldTurnLayerToTets(bool b) { ops.shouldTurnLayerToTets = b; }
    /** \brief getter for shouldCleanupLayer
        \details see InputOptions */
    bool shouldCleanupLayer() { return ops.shouldCleanupLayer; }
    /** \brief setter for shouldCleanupLayer */
    void shouldCleanupLayer(bool b) { ops.shouldCleanupLayer = b; }
    /** \brief getter for shouldRefineLayer
        \details see InputOptions */
    bool shouldRefineLayer() { return ops.shouldRefineLayer; }
    /** \brief setter for shouldRefineLayer */
    void shouldRefineLayer(bool b) { ops.shouldRefineLayer = b; }
    /** \brief getter for shouldCoarsenLayer
        \details see InputOptions */
    bool shouldCoarsenLayer() { return ops.shouldCoarsenLayer; }
    /** \brief setter for shouldCoarsenLayer */
    void shouldCoarsenLayer(bool b) { ops.shouldCoarsenLayer = b; }
    /** \brief getter for splitAllLayerEdges
        \details see InputOptions */
    bool splitAllLayerEdges() { return ops.splitAllLayerEdges; }
    /** \brief setter for splitAllLayerEdges */
    void splitAllLayerEdges(bool b) { ops.splitAllLayerEdges = b; }
    /** \brief getter for userDefinedLayerTagName
        \details see InputOptions */
    const char* userDefinedLayerTagName() { return ops.userDefinedLayerTagName; }
    /** \brief setter for userDefinedLayerTagName */
    void userDefinedLayerTagName(const char* c) { ops.userDefinedLayerTagName = c; }
    /** \brief getter for debugFolder
        \details see InputOptions */
    const char* debugFolder() { return ops.debugFolder; }
    /** \brief setter for debugFolder */
    void debugFolder(const char* c) { ops.debugFolder = c; }

    void validateInput();
    Mesh* mesh;
    SizeField* sizeField;
    bool ownsSizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSolutionTransfer;
    ShapeHandlerFunction shapeHandler;
  private:
    Input(const Input& in);
    void init(Mesh* m, SolutionTransfer* s);
    void setDefaultValues();
    void updateMaxIterBasedOnSize();
    InputOptions ops;
};

/** \brief removes the constantness of in
 \details users can call this on a const Input* to make an Input* object that can be modified
          after the call to configure functions.
 \param in input const Input pointer */
Input* makeAdvanced(const Input* in);

/** \brief generate a configuration based on an anisotropic function.
 \details the return Input object is un-modifiable
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
const Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
/** \brief generate a configuration based on an isotropic function
 \details the return Input object is un-modifiable
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
const Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on anisotropic fields
 \details the return Input object is un-modifiable
 \param sizes a vector field of desired element sizes along the
              axes of the anisotropic frame
 \param frames a matrix field containing anisotropic frames
               for each vertex
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
 \details the return Input object is un-modifiable
 \param size a scalar field of desired element size
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
const Input* configure(
    Mesh* m,
    apf::Field* size,
    SolutionTransfer* s=0);
/** \brief generate a uniform refinement configuration (un-modifiable Input) */
const Input* configureUniformRefine(Mesh* m, int n=1, SolutionTransfer* s=0);

/** \brief generate a matched uniform refinement configuration  (un-modifiable Input) */
const Input* configureMatching(Mesh* m, int n=1, SolutionTransfer* s=0);

/** \brief generate a no-op configuration (un-modifiable Input) */
const Input* configureIdentity(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);

}

#endif
