/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef MACURVEMESH_H
#define MACURVEMESH_H

#include "maMesh.h"
#include "apfShape.h"

namespace ma {

class Adapt;

class MeshCurver
{
  public:
    MeshCurver(Mesh* m, int P, int B);
    virtual ~MeshCurver() {};
    virtual bool run() = 0;

    /** \brief snaps points to interpolating locations */
    void snapToInterpolate(int dim);

    /** \brief these two are a per entity version of above */
    void snapToInterpolateEdge(Entity* e);
    void snapToInterpolateTri(Entity* e);

    /** \brief converts interpolating points to control points */
    void convertInterpolationPoints(Entity* e, int n, int ne,
      apf::NewArray<double>& c);

  protected:
    Mesh* m_mesh;
    int m_order;
    int m_blendOrder;
};

class BezierCurver : public MeshCurver
{
  public:
    BezierCurver(Mesh* m, int P, int B) : MeshCurver(m, P, B) {};

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see apfBezier.cc */
    virtual bool run();

};

class GregoryCurver : public MeshCurver
{
  public:
    GregoryCurver(Mesh* m, int P, int B) : MeshCurver(m, P, B) {};
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
/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling
   n points per dimension of the entity through uniform
   sampling locations in parameter space */
double interpolationError(Mesh* m, Entity* e, int n);

/** \brief Mostly a debugging function, writes csv file of n points per dim */
void writePointSet(Mesh* m, int d, int n, const char* prefix);

}

#endif
