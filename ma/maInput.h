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
    Input(Input* in) {this->ops = in->ops;}

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

    // setters and getters
/** \brief number of refine/coarsen iterations to run (default 3)
    \details this value will be set to the minimum required iterations
    inside the call to ma::configure in cases where there is a size information.
    Users can override this by setting in->maximumIterations after the
    call to ma::configure and before the call to ma::adapt routine.*/
    int maximumIterations();
    void maximumIterations(int i);
/** \brief whether to perform the collapse step */
    bool shouldCoarsen();
    void shouldCoarsen(bool b);
/** \brief whether to snap new vertices to the model surface
    \details requires modeler support, see gmi_can_eval */
    bool shouldSnap();
    void shouldSnap(bool b);
/** \brief whether to transfer parametric coordinates
    \details requires modeler support, see gmi_reparam */
    bool shouldTransferParametric();
    void shouldTransferParametric(bool b);
/** \brief whether to transfer to the parametric coords of the closest point
    \details requires modeler support, see gmi_closest_point */
    bool shouldTransferToClosestPoint();
    void shouldTransferToClosestPoint(bool b);
/** \brief whether to update matched entity info (limited support) */
    bool shouldHandleMatching();
    void shouldHandleMatching(bool b);
/** \brief whether to run shape correction (default true) */
    bool shouldFixShape();
    void shouldFixShape(bool b);
/** \brief whether to adapt if it makes local quality worse (default false) */
    bool shouldForceAdaptation();
    void shouldForceAdaptation(bool b);
/** \brief whether to print the worst shape quality */
    bool shouldPrintQuality();
    void shouldPrintQuality(bool b);
/** \brief minimum desired mean ratio cubed for simplex elements
    \details a different measure is used for curved elements */
    double goodQuality();
    void goodQuality(double d);
/** \brief whether to check the quality of split elems in DoubleSplitsCollapse
    (default false) */
    bool shouldCheckQualityForDoubleSplits();
    void shouldCheckQualityForDoubleSplits(bool b);
/** \brief minimum valid mean ratio cubed for simplex elements (default 1e-10)
    \details used to define inside-out tetrahedra.
    a different measure is used for curved elements */
    double validQuality();
    void validQuality(double d);
/** \brief imbalance target for all load balancing tools (default 1.10) */
    double maximumImbalance();
    void maximumImbalance(double d);
/** \brief whether to run zoltan predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreZoltan();
    void shouldRunPreZoltan(bool b);
/** \brief whether to run zoltan predictive load balancing using RIB (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreZoltanRib();
    void shouldRunPreZoltanRib(bool b);
/** \brief whether to run parma predictive load balancing (default false)
    \details if this and all the other PreBalance options are false, pre-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPreParma();
    void shouldRunPreParma(bool b);
/** \brief whether to run zoltan during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunMidZoltan();
    void shouldRunMidZoltan(bool b);
/** \brief whether to run parma during adaptation (default false)
    \details if this and all the other MidBalance options are false, mid-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunMidParma();
    void shouldRunMidParma(bool b);
/** \brief whether to run zoltan after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPostZoltan();
    void shouldRunPostZoltan(bool b);
/** \brief whether to run zoltan RIB after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPostZoltanRib();
    void shouldRunPostZoltanRib(bool b);
/** \brief whether to run parma after adapting (default false)
    \details if this and all the other PostBalance options are false, post-balancing
    occurs only if the imbalance is greater than in->maximumImbalance */
    bool shouldRunPostParma();
    void shouldRunPostParma(bool b);
/** \brief the ratio between longest and shortest edges that differentiates a
    "short edge" element from a "large angle" element. */
    double maximumEdgeRatio();
    void maximumEdgeRatio(double d);
/** \brief whether to tetrahedronize the boundary layer (default false)  */
    bool shouldTurnLayerToTets();
    void shouldTurnLayerToTets(bool b);
/** \brief whether to tetrahedronize abnormal pyramids (default false) */
    bool shouldCleanupLayer();
    void shouldCleanupLayer(bool b);
/** \brief whether to allow layer refinement (default false) */
    bool shouldRefineLayer();
    void shouldRefineLayer(bool b);
/** \brief whether to allow layer coarsening (default false) */
    bool shouldCoarsenLayer();
    void shouldCoarsenLayer(bool b);
/** \brief set to true during UR to get splits in the normal direction */
    bool splitAllLayerEdges();
    void splitAllLayerEdges(bool b);
/** \brief the name of the (user defined) INT tag specifying the boundary
    layer elements. Use the value of 0 for non-layer elements and a non-zero value
    for layer elements. (default "") */
    const char* userDefinedLayerTagName();
    void userDefinedLayerTagName(const char* c);
/** \brief this a folder that debugging meshes will be written to, if provided! */
    const char* debugFolder();
    void debugFolder(const char* c);

    virtual bool modifiable() = 0;
    void validateInput();

    Mesh* mesh;
    SizeField* sizeField;
    bool ownsSizeField;
    SolutionTransfer* solutionTransfer;
    bool ownsSolutionTransfer;
    ShapeHandlerFunction shapeHandler;
  private:
    void init(Mesh* m, SolutionTransfer* s);
    void setDefaultValues();
    void updateMaxIterBasedOnSize();
    InputOptions ops;
};

class InputBasic : public Input
{
  public:
    ~InputBasic() {}
    InputBasic(Input* in) : Input(in) {}
    InputBasic(
      Mesh* m,
      AnisotropicFunction* f,
      SolutionTransfer* s=0,
      bool logInterpolation=true) : Input(m,f,s,logInterpolation) {}
    InputBasic(
      Mesh* m,
      IsotropicFunction* f,
      SolutionTransfer* s=0) : Input(m,f,s) {}
    InputBasic(
      Mesh* m,
      apf::Field* sizes,
      apf::Field* frames,
      SolutionTransfer* s=0,
      bool logInterpolation=true) : Input(m,sizes,frames,s,logInterpolation) {}
    InputBasic(
      Mesh* m,
      apf::Field* size,
      SolutionTransfer* s=0) : Input(m,size,s) {}
    InputBasic(Mesh* m, int n=1, SolutionTransfer* s=0) : Input(m,n,s) {}
    InputBasic(Mesh* m, SizeField* f=0, SolutionTransfer* s=0) : Input(m,f,s) {}
    bool modifiable() { return false; }
};

class InputAdvanced : public Input
{
  public:
    ~InputAdvanced() {}
    InputAdvanced(Input* in) : Input(in) {}
    InputAdvanced(
      Mesh* m,
      AnisotropicFunction* f,
      SolutionTransfer* s=0,
      bool logInterpolation=true) : Input(m,f,s,logInterpolation) {}
    InputAdvanced(
      Mesh* m,
      IsotropicFunction* f,
      SolutionTransfer* s=0) : Input(m,f,s) {}
    InputAdvanced(
      Mesh* m,
      apf::Field* sizes,
      apf::Field* frames,
      SolutionTransfer* s=0,
      bool logInterpolation=true) : Input(m,sizes,frames,s,logInterpolation) {}
    InputAdvanced(
      Mesh* m,
      apf::Field* size,
      SolutionTransfer* s=0) : Input(m,size,s) {}
    InputAdvanced(Mesh* m, int n=1, SolutionTransfer* s=0) : Input(m,n,s) {}
    InputAdvanced(Mesh* m, SizeField* f=0, SolutionTransfer* s=0) : Input(m,f,s) {}
    bool modifiable() { return true; }
};




/** \brief generate a configuration based on an anisotropic function.
 \details the return Input object is un-modifiable
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
Input* configure(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
/** \brief generate a configuration based on an anisotropic function.
 \details the return Input object is modifiable (for advanced users)
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
Input* configureAdvanced(
    Mesh* m,
    AnisotropicFunction* f,
    SolutionTransfer* s=0,
    bool logInterpolation=true);

/** \brief generate a configuration based on an isotropic function
 \details the return Input object is un-modifiable
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configure(
    Mesh* m,
    IsotropicFunction* f,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on an isotropic function
 \details the return Input object is modifiable (for advanced users)
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configureAdvanced(
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
Input* configure(
    Mesh* m,
    apf::Field* sizes,
    apf::Field* frames,
    SolutionTransfer* s=0,
    bool logInterpolation=true);
/** \brief generate a configuration based on anisotropic fields
 \details the return Input object is modifiable (for advanced users)
 \param sizes a vector field of desired element sizes along the
              axes of the anisotropic frame
 \param frames a matrix field containing anisotropic frames
               for each vertex
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms
 \param logInterpolation if true uses logarithmic interpolation */
Input* configureAdvanced(
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
Input* configure(
    Mesh* m,
    apf::Field* size,
    SolutionTransfer* s=0);
/** \brief generate a configuration based on an isotropic field
 \details the return Input object is modifiable (for advanced users)
 \param size a scalar field of desired element size
 \param s if non-zero, use that to transfer all fields. otherwise,
          transfer any associated fields with default algorithms */
Input* configureAdvanced(
    Mesh* m,
    apf::Field* size,
    SolutionTransfer* s=0);

/** \brief generate a uniform refinement configuration (un-modifiable Input) */
Input* configureUniformRefine(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a uniform refinement configuration (modifiable Input) */
Input* configureUniformRefineAdvanced(Mesh* m, int n=1, SolutionTransfer* s=0);

/** \brief generate a matched uniform refinement configuration  (un-modifiable Input) */
Input* configureMatching(Mesh* m, int n=1, SolutionTransfer* s=0);
/** \brief generate a matched uniform refinement configuration (modifiable Input)  */
Input* configureMatchingAdvanced(Mesh* m, int n=1, SolutionTransfer* s=0);

/** \brief generate a no-op configuration (un-modifiable Input) */
Input* configureIdentity(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);
/** \brief generate a no-op configuration (modifiable Input)  */
Input* configureIdentityAdvanced(Mesh* m, SizeField* f=0, SolutionTransfer* s=0);

// /**  \brief used internally, but users can call this if they want */
// void validateInput(Input* in);

}

#endif
