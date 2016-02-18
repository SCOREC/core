/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRV_H
#define CRV_H

#include "apfMesh2.h"
#include "apfShape.h"
#include <ma.h>
#include <stdio.h>

/** \file crv.h
  * \brief main file for curved element support,
  * only one accessible from install right now */

/** \namespace crv
  * \brief the curving functions are contained in this namespace */
namespace crv {

/** \brief actually 1 greater than max order */
static unsigned const MAX_ORDER = 19;

/** \brief sets order used in bezier shape functions */
void setOrder(const int order);
/** \brief gets order used in bezier shape functions */
int getOrder();
/** \brief sets the blending order, if shape blending is used */
void setBlendingOrder(const int type, const int b);
/** \brief gets the blending order */
int getBlendingOrder(const int type);

/** \brief computes min det Jacobian / max det Jacobian */
double getQuality(apf::Mesh* m,apf::MeshEntity* e);

/** \brief change the order of a Bezier Mesh
 * \details going up in order is exact,
 * except for boundary elements, where snapping changes things
 * Going down in order is approximate everywhere
 * */
void changeMeshOrder(apf::Mesh2* m, int newOrder);

/** \brief Base Mesh curving object
  \details P is the order, S is the space dimension,
  different from the mesh dimension, used to distinguish between planar 2D
  meshes and surface meshes. */

class MeshCurver
{
  public:
    MeshCurver(apf::Mesh2* m, int P) : m_mesh(m), m_order(P) {};
    virtual ~MeshCurver() {};
    virtual bool run() = 0;

    /** \brief snaps points to interpolating locations */
    void snapToInterpolate(int dim);

    /** \brief wrapper around synchronizeFieldData */
    void synchronize();

  protected:
    apf::Mesh2* m_mesh;
    int m_order;
};

/** \brief curves an already changed mesh
 * \details this is a bit of a hack, meant to work with
 * all interpolating shape functions
 */
class InterpolatingCurver : public MeshCurver
{
  public:
    InterpolatingCurver(apf::Mesh2* m, int P) : MeshCurver(m,P) {};
    virtual ~InterpolatingCurver() {};
    virtual bool run();

};

/** \brief this curves a mesh with Bezier shapes
 * \details converts the mesh and snaps boundary entities to geometry
 * P is the order, B is the blending order (set to 0 to use full shapes)
 */
class BezierCurver : public MeshCurver
{
  public:
    BezierCurver(apf::Mesh2* m, int P, int B) : MeshCurver(m,P)
    {
      setBlendingOrder(apf::Mesh::TYPES,B);
    };

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see crvBezier.cc */
    virtual bool run();
    /** \brief converts interpolating points to bezier control points */
    void convertInterpolatingToBezier();
};

/** \brief this curves a mesh with 4th order G1 Patches
 * \details approach is in crvBezier.cc
 */
class GregoryCurver : public BezierCurver
{
  public:
    GregoryCurver(apf::Mesh2* m, int P, int B)
    : BezierCurver(m,P,B) {};
    /** \brief curves a mesh using G1 gregory surfaces, see crvBezier.cc */
    virtual bool run();
    /** \brief sets cubic edge points using normals */
    void setCubicEdgePointsUsingNormals();
    /** \brief sets internal points locally */
    void setInternalPointsLocally();
};

/** \brief Get the Bezier Curve or Shape of some order
 \details goes from first to sixth order */
apf::FieldShape* getBezier(int order);
/** \brief Get the 4th order Gregory Surface*/
apf::FieldShape* getGregory();

/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling
   n points per dimension of the entity through uniform
   sampling locations in parameter space */
double interpolationError(apf::Mesh* m, apf::MeshEntity* e, int n);

/** \brief Visualization, writes file for specified type */
void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix);

/** \brief Visualization, writes file of control nodes for each entity */
void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix);
/** \brief Visualization, writes file of shapes evaluated at node xi
    for each entity */
void writeInterpolationPointVtuFiles(apf::Mesh* m, const char* prefix);

/** \brief publically accessible functions */
int getTriNodeIndex(int P, int i, int j);
int getTetNodeIndex(int P, int i, int j, int k);

/** \brief check the validity (det(Jacobian) > eps) of an element
 * \details entities is a container of invalid downward entities
 * algorithm is an integer corresponding to what method to use
 * 0 - subdivision, without first check
 * 1 - elevation, without first check
 * 2 - subdivision
 * 3 - elevation
 * 4 - subdivision, using matrices
 * methods 2 and 3 exist because the first check tends to catch everything
 * without actually using subdivision and elevation, and giving this option
 * is easier for debugging and verifying the efficacy of those procedures
 * */
int checkValidity(apf::Mesh* m, apf::MeshEntity* e,
    int algorithm = 4);
/** \brief count invalid elements of the mesh */
int countNumberInvalidElements(apf::Mesh2* m);

/** \brief crv fail function */
void fail(const char* why) __attribute__((noreturn));

} //namespace crv

#endif
