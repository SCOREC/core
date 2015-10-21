/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crvBezier.h"
#include "crvTables.h"
#include "crvMath.h"
#include <mth_def.h>
#include <cassert>

namespace crv {

/*
 * For anything not using optimal points, just allocate as needed
 */
static void getHigherBezierNodeXi(int type, int P, int node, apf::Vector3& xi)
{
  // 19th order supported
  static apf::NewArray<double> edgePoints[20];
  static apf::NewArray<apf::Vector3> triPoints[20];
  static apf::NewArray<apf::Vector3> tetPoints[20];
  // if edges are allocated
  if (!edgePoints[P].allocated()){
    edgePoints[P].allocate(P-1);
    double dp = 2./P;
    for (int i = 1; i <= P-1; ++i)
      edgePoints[P][i-1] = -1.+dp*i;
    triPoints[P].allocate((P-1)*(P-2)/2);
    dp = 1./P;
    int start = getNumControlPoints(apf::Mesh::TRIANGLE,P)
        - getNumInternalControlPoints(apf::Mesh::TRIANGLE,P);
    for(int j = 1; j <= P-2; ++j)
      for(int i = 1; i <= P-1-j; ++i){
        int index = computeTriNodeIndex(P,i,j) - start;
        triPoints[P][index][0] = dp*j;
        triPoints[P][index][1] = dp*(P-i-j);
      }
    start = getNumControlPoints(apf::Mesh::TET,P)
        - getNumInternalControlPoints(apf::Mesh::TET,P);
    tetPoints[P].allocate((P-1)*(P-2)*(P-3)/6);
    for (int k = 1; k <= P-2; ++k)
      for (int j = 1; j <= P-k-2; ++j)
        for (int i = 1; i <= P-j-k-1; ++i){
          int index = computeTetNodeIndex(P,i,j,k) - start;
          tetPoints[P][index] = apf::Vector3(dp*j,dp*k,dp*(P-i-j-k));
        }
  }
  if(type == apf::Mesh::EDGE)
    xi[0] = edgePoints[P][node];
  if(type == apf::Mesh::TRIANGLE)
    xi = triPoints[P][node];
  if(type == apf::Mesh::TET){
    xi = tetPoints[P][node];
  }
}

void getBezierNodeXi(int type, int P, int node, apf::Vector3& xi)
{
  static double eP2[1] = {0.0};
  static double eP3[2] = {-0.4503914,0.4503914};
  static double eP4[3] = {-0.6612048,0.0,0.6612048};
  static double eP5[4] = {-0.7732854,-0.2863522,0.2863522,0.7732854};
  static double eP6[5] = {-0.8388042,-0.469821,0.0,
      0.469821,0.8388042};
  static double* edgePoints[6] =
  {eP2, eP2, eP3, eP4, eP5, eP6 };
  static apf::Vector3 triPoints5[6] =
  {apf::Vector3(0.15251715,0.15251715,0.6949657),
   apf::Vector3(0.4168658,0.1662684,0.4168658),
   apf::Vector3(0.6949657,0.15251715,0.15251715),
   apf::Vector3(0.1662684,0.4168658,0.4168658),
   apf::Vector3(0.4168658,0.4168658,0.1662684),
   apf::Vector3(0.15251715,0.6949657,0.15251715)};
  static apf::Vector3 triPoints6[10] =
  {apf::Vector3(0.10971385,0.10971385,0.7805723),
   apf::Vector3(0.3157892,0.1256031,0.5586077),
   apf::Vector3(0.5586077,0.1256031,0.3157892),
   apf::Vector3(0.7805723,0.10971385,0.10971385),
   apf::Vector3(0.1256031,0.3157892,0.5586077),
   apf::Vector3(1./3.,1./3.,1./3.),
   apf::Vector3(0.5586077,0.3157892,0.1256031),
   apf::Vector3(0.1256031,0.5586077,0.3157892),
   apf::Vector3(0.3157892,0.5586077,0.1256031),
   apf::Vector3(0.10971385,0.7805723,0.10971385)};

  switch (type) {
  case apf::Mesh::EDGE:
    if(P <= 6)
      xi[0] = edgePoints[P-1][node];
    else
      getHigherBezierNodeXi(type,P,node,xi);
    break;
  case apf::Mesh::TRIANGLE:
    // technically only two of these numbers are needed
    switch (P) {
    case 1:
    case 2:
      fail("expected P >= 3");
    case 3:
      xi = apf::Vector3(1./3.,1./3.,1./3.);
      break;
    case 4:
      xi[(node+2) % 3] = 0.5582239;
      xi[(node+0) % 3] = 0.22088805;
      xi[(node+1) % 3] = 0.22088805;
      break;
    case 5:
      xi = triPoints5[node];
      break;
    case 6:
      xi = triPoints6[node];
      break;
    default:
      getHigherBezierNodeXi(type,P,node,xi);
      break;
    }
    break;
  case apf::Mesh::TET:
    switch (P) {
    case 1:
    case 2:
    case 3:
      fail("expected P > 3");
    case 4:
      xi = apf::Vector3(0.25,0.25,0.25);
      break;
    default:
      getHigherBezierNodeXi(type,P,node,xi);
      break;
    }
    break;
    default:
      xi.zero();
      break;
  }
}

void getBezierTransformationCoefficients(apf::Mesh* m, int P, int type,
    apf::NewArray<double> & c)
{
//  if (P <= 6) return;
  int nb = getNumInternalControlPoints(type,P);
  int ni = getNumControlPoints(type,P);

  static apf::NewArray<double> transform[apf::Mesh::TYPES][20];
  if(!transform[type][P].allocated()){

    int oldB[apf::Mesh::TYPES];
    for (int i = 0; i < apf::Mesh::TYPES; ++i)
      oldB[i] = getBlendingOrder(i);

    setBlendingOrder(apf::Mesh::TYPES,0);

    transform[type][P].allocate(nb*ni);
    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[type]);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))){
      if (m->getType(e) == type)
        break;
    }
    m->end(it);
    mth::Matrix<double> A(ni,ni);
    mth::Matrix<double> Ai(ni,ni);
    getTransformationMatrix(m,e,A);
    invertMatrix(ni,A,Ai);

    for( int i = 0; i < nb; ++i)
      for( int j = 0; j < ni; ++j)
        transform[type][P][i*ni+j] = Ai(i+ni-nb,j);

    for (int i = 0; i < apf::Mesh::TYPES; ++i)
      setBlendingOrder(i,oldB[i]);
  }

  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = transform[type][P][i*ni+j];

}

void getInternalBezierTransformationCoefficients(apf::Mesh* m, int P, int blend,
    int type, apf::NewArray<double> & c)
{
  // this is fun. To compute this, the actual type must be
  // blended, but the lower entities must not be.
  // We also require the full matrix from above, call it
  // A, and the blended matrix B
  // The return is inv(A)*B

  int nb = getNumInternalControlPoints(type,P);
  int ni = getNumControlPoints(type,P);
  static apf::NewArray<double> transform[2][apf::Mesh::TYPES][20];

  if(!transform[blend-1][type][P].allocated()){

    int oldB[apf::Mesh::TYPES];
    for (int i = 0; i < apf::Mesh::TYPES; ++i)
      oldB[i] = getBlendingOrder(i);
    // get first Matrix
    setBlendingOrder(apf::Mesh::TYPES,0);
    transform[blend-1][type][P].allocate(nb*(ni-nb));

    apf::MeshIterator* it = m->begin(apf::Mesh::typeDimension[type]);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))){
      if (m->getType(e) == type &&
          m->getModelType(m->toModel(e)) == m->getDimension())
        break;
    }
    m->end(it);

    mth::Matrix<double> A(ni,ni);
    mth::Matrix<double> Ai(ni,ni);
    mth::Matrix<double> B(ni,ni);

    getTransformationMatrix(m,e,A);
    invertMatrix(ni,A,Ai);

    // now get second matrix
    setBlendingOrder(type,blend);
    getTransformationMatrix(m,e,B);

    // fill in the last few rows of B
    apf::NewArray<double> values;
    apf::Vector3 xi;
    for (int i = 0; i < nb; ++i){
      m->getShape()->getNodeXi(type,i,xi);
      m->getShape()->getEntityShape(type)->getValues(m,e, xi,values);
      for (int j = 0; j < ni; ++j)
        B(i+ni-nb,j) = values[j];
    }

    mth::multiply(Ai,B,A);

    for( int i = 0; i < nb; ++i)
      for( int j = 0; j < ni-nb; ++j)
        transform[blend-1][type][P][i*(ni-nb)+j] = A(i+ni-nb,j);

    for (int i = 0; i < apf::Mesh::TYPES; ++i)
      setBlendingOrder(i,oldB[i]);
  }
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = transform[blend-1][type][P][i*ni+j];

}

} // namespace crv
