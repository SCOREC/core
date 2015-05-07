#ifndef MACURVEMESH_H
#define MACURVEMESH_H

#include "maMesh.h"

namespace ma {

class Adapt;

class MeshCurver
{
  public:
    MeshCurver(Adapt* a, int order);
    virtual ~MeshCurver() {};
    virtual bool run() = 0;
    void snapToInterpolate(int dim);
    Adapt* adapt;
    int order;
};

class BezierCurver : public MeshCurver
{
  public:
    BezierCurver(Adapt* a, int o) : MeshCurver(a, o) {};

    /** \brief curves a mesh using bezier curves of chosen order
      \details finds interpolating points, then converts to control points
      see apfBezier.cc */
    virtual bool run();

  /** \brief converts interpolating points to bezier points */
    void convertInterpolationPoints(Entity* e, int n, int ne,
      apf::NewArray<double>& c);
};

/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling n points. */
double interpolationError(Mesh* m, Entity* e, int n);

/** \brief Mostly a debugging function, writes csv file of n points per dim */
void writePointSet(Mesh* m, int d, int n, const char* prefix);

}

#endif
