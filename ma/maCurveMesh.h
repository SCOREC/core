#ifndef MACURVEMESH_H
#define MACURVEMESH_H

#include "maMesh.h"

namespace ma {

/** \brief computes interpolation error of a curved entity on a mesh
  \details this computes the Hausdorff distance by sampling n points. */
double interpolationError(Mesh* m, Entity* e, int n,
    Vector &samplept, Vector &maxpt);
/** \brief computes interpolation error at nodeXi. Should always be small
  \details this computes the Hausdorff distance by sampling n points. */
double interpolationErrorAtNodeXi(Mesh* m, Entity* e, int n,
    Vector &samplept, Vector &maxpt);

/** \brief curves a mesh using bezier curves of chosen order
  \details finds interpolating points, then converts to control points
  see apfBezier.cc */
void curveMeshToBezier(Mesh* m, int order);

/** \brief Mostly a debugging function, writes csv file of n points per dim */
void writePointSet(Mesh* m, int d, int n, const char* prefix);

}

#endif
