/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRV_H
#define CRV_H

#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"

namespace crv {

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
    BezierCurver(apf::Mesh2* m, int P, int B) : MeshCurver(m, P), m_blendOrder(B) {};

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see apfBezier.cc */
    virtual bool run();
  protected:
    int m_blendOrder;
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

/** \brief Get the Bezier Curve or Shape of some order
 \details goes from first to sixth order */
apf::FieldShape* getBezier(int dimension, int order, int blendOrder);
/** \brief Get the Gregory Surface of some order
 \details only fourth order right now*/
apf::FieldShape* getGregory(int order, int blendOrder);

/** \brief get coefficients for interpolating points to control points
 \details works only for prescribed optimal point locations */
void getTransformationCoefficients(int dim, int type,
    apf::NewArray<double>& c);

/** \brief binomial function */
int binomial(int n, int i);

/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling
   n points per dimension of the entity through uniform
   sampling locations in parameter space */
double interpolationError(apf::Mesh2* m, apf::MeshEntity* e, int n);

/** \brief Visualization, writes two files for edges and faces */
void writeCurvedVtuFiles(apf::Mesh* m, int n, const char* prefix);

}

#endif
