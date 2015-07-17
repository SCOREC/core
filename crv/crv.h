/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRV_H
#define CRV_H

#include "apfMesh2.h"
#include "apfShape.h"

namespace crv {

/** \brief sets the blending order, if shape blending is used */
void setBlendingOrder(const int b);
/** \brief gets the blending order */
int getBlendingOrder();

class MeshCurver
{
  public:
    MeshCurver(apf::Mesh2* m, int P) : m_mesh(m), m_order(P) {};
    virtual ~MeshCurver() {};
    virtual bool run() = 0;

    /** \brief snaps points to interpolating locations */
    void snapToInterpolate(int dim);

    /** \brief these two are a per entity version of above */
    void snapToInterpolateEdge(apf::MeshEntity* e);
    void snapToInterpolateTri(apf::MeshEntity* e);

    /** \brief converts interpolating points to control points */
    void convertInterpolationPoints(apf::MeshEntity* e, int n, int ne,
      apf::NewArray<double>& c);

  protected:
    apf::Mesh2* m_mesh;
    int m_order;
};

/** \brief curves an already changed mesh */
class InterpolatingCurver : public MeshCurver
{
  public:
    InterpolatingCurver(apf::Mesh2* m, int P) : MeshCurver(m,P) {};
    virtual ~InterpolatingCurver() {};
    virtual bool run();

};

class BezierCurver : public MeshCurver
{
  public:
    BezierCurver(apf::Mesh2* m, int P, int B) : MeshCurver(m, P)
    { setBlendingOrder(B); };

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see apfBezier.cc */
    virtual bool run();
};

class GregoryCurver : public BezierCurver
{
  public:
    GregoryCurver(apf::Mesh2* m, int P, int B) : BezierCurver(m, P, B) {};
    /** \brief curves a mesh using G1 gregory surfaces, see apfBezier.cc */
    virtual bool run();
    /** \brief sets cubic edge points using normals */
    void setCubicEdgePointsUsingNormals();
    /** \brief sets internal points using neighbors (See Notes)
      \details NOT CURRENTLY FULLY IMPLEMENTED */
    void setInternalPointsUsingNeighbors();
    /** \brief sets internal points locally (4th order only) */
    void setInternalPointsLocally();
};

class SphereCurver : public MeshCurver
{
  public:
    SphereCurver(apf::Mesh2* m, int P, int B) : MeshCurver(m, P),
    m_blendOrder(B) {};

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see apfBezier.cc */
    virtual bool run();
  protected:
    int m_blendOrder;
};
/** \brief Elevate a bezier curve to a higher order
 \details This elevates from nth order to n+rth order
 requires the curve be order n+r in memory already, and
 that the first n points correspond to the lower order curve */
void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r);

/** \brief Get the Bezier Curve or Shape of some order
 \details goes from first to sixth order */
apf::FieldShape* getBezier(int dimension, int order);
/** \brief Get the Gregory Surface of some order
 \details only fourth order right now*/
apf::FieldShape* getGregory(int order);
/** \brief Get the NURBS, based off of bezier
 \details goes from first to sixth order */
apf::FieldShape* getNurbs(int orde);
/** \brief set the weights
 \details used to set these for every curved surface*/
void setNurbsEdgeWeights(apf::NewArray<double>& weights);
void setNurbsTriangleWeights(apf::NewArray<double>& weights);

/** \brief get coefficients for interpolating points to control points
 \details works only for prescribed optimal point locations */
void getTransformationCoefficients(int P, int dim, int type,
    apf::NewArray<double>& c);

/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling
   n points per dimension of the entity through uniform
   sampling locations in parameter space */
double interpolationError(apf::Mesh* m, apf::MeshEntity* e, int n);

/** \brief Visualization, writes file for specified type */
void writeCurvedVtuFiles(apf::Mesh* m, int type, int n, const char* prefix);

/** \brief Visualization, writes file of control nodes for each entity */
void writeControlPointVtuFiles(apf::Mesh* m, const char* prefix);

/** \brief binomial functions */
int binomial(int n, int i);
int trinomial(int n, int i, int j);
int quadnomial(int n, int i, int j, int k);

/** \brief crv fail function */
void fail(const char* why) __attribute__((noreturn));

} //namespace crv

#endif
