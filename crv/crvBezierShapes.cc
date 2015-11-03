/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "crv.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"

namespace crv {

static void bezierCurve(int P, apf::Vector3 const& xi,
    apf::NewArray<double>& values)
{
  double t = 0.5*(xi[0]+1.);
  for(int i = 1; i < P; ++i)
    values[i+1] = binomial(P,i)*Bij(P-i,i,1.-t,t);
  values[0] = intpow(1-t, P);
  values[1] = intpow(t, P);
}

static void bezierCurveGrads(int P, apf::Vector3 const& xi,
    apf::NewArray<apf::Vector3>& grads)
{
  double t = 0.5*(xi[0]+1.);
  for(int i = 1; i < P; ++i)
    grads[i+1] = apf::Vector3(binomial(P,i)*(i-P*t)
        *Bij(P-1-i,i-1,1.-t,t)/2.,0,0);
  grads[0] = apf::Vector3(-P*intpow(1-t, P-1)/2.,0,0);
  grads[1] = apf::Vector3(P*intpow(t, P-1)/2.,0,0);
}

static void bezierTriangle(int P, apf::Vector3 const& xi,
    apf::NewArray<double>& values)
{
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
  for(int i = 0; i < P+1; ++i)
    for(int j = 0; j < P+1-i; ++j)
      values[getTriNodeIndex(P,i,j)] =
          trinomial(P,i,j)*Bijk(i,j,P-i-j,xii[0],xii[1],xii[2]);
}

static void bezierTriangleGrads(int P, apf::Vector3 const& xi,
    apf::NewArray<apf::Vector3>& grads)
{

  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  apf::Vector3 gxii[3] =
  {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

  for(int i = 0; i < 3; ++i)
    grads[i] = gxii[i]*P*intpow(xii[i],P-1);

  for(int i = 1; i < P+1; ++i)
    for(int j = 1; j < P-i; ++j)
      grads[getTriNodeIndex(P,i,j)] =
          gxii[0]*trinomial(P,i,j)*(i*(1.-xii[1])-(P-j)*xii[0])
          *Bijk(i-1,j,P-i-j-1,xii[0],xii[1],xii[2]) +
          gxii[1]*trinomial(P,i,j)*(j*(1.-xii[0])-(P-i)*xii[1])
          *Bijk(i,j-1,P-i-j-1,xii[0],xii[1],xii[2]);

  // i = 0
  for(int j = 1; j < P; ++j)
    grads[3+2*(P-1)-j] =
        (gxii[0]*(j-P)*xii[1] +
            gxii[1]*(j*(1.-xii[0])-P*xii[1]))
            *binomial(P,j)*Bij(j-1,P-j-1,xii[1],xii[2]);

  // j = 0
  for(int i = 1; i < P; ++i)
    grads[3+2*(P-1)-1+i] =
        (gxii[0]*(i*(1.-xii[1])-P*xii[0]) +
            gxii[1]*(i-P)*xii[0])
            *binomial(P,i)*Bij(i-1,P-i-1,xii[0],xii[2]);

  // k = 0
  for(int i = 1, j = P-1; i < P; ++i, --j)
    grads[3+(P-1)-i] =
        (gxii[0]*i*xii[1] + gxii[1]*j*xii[0])
        *trinomial(P,i,j)*Bij(i-1,j-1,xii[0],xii[1]);
}

static void bezierTet(int P, apf::Vector3 const& xi,
    apf::NewArray<double>& values)
{
  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
  for(int i = 0; i < 4; ++i)
    values[i] = intpow(xii[i],P);

  int nE = P-1;

  int const (*tev)[2] = apf::tet_edge_verts;

  for(int a = 0; a < 6; ++a)
    for(int b = 0; b < nE; ++b) // edge nodes
      values[4+a*nE+b] = binomial(P,b+1)
      *Bij(P-b-1,b+1,xii[tev[a][0]],xii[tev[a][1]]);

  // face 0, l = 0
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j)
      values[computeTetNodeIndex(P,i,j,P-i-j)] = trinomial(P,i,j)
      *Bijk(i,j,P-i-j,xii[0],xii[1],xii[2]);
  // face 1, k = 0
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j)
      values[computeTetNodeIndex(P,i,j,0)] = trinomial(P,i,j)
      *Bijk(i,j,P-i-j,xii[0],xii[1],xii[3]);
  // face 2, i = 0
  for(int j = 1; j <= P-1; ++j)
    for(int k = 1; k <= P-1-j; ++k)
      values[computeTetNodeIndex(P,0,j,k)] = trinomial(P,j,k)
      *Bijk(j,k,P-j-k,xii[1],xii[2],xii[3]);
  // face 3, j = 0
  for(int i = 1; i <= P-1; ++i)
    for(int k = 1; k <= P-1-i; ++k)
      values[computeTetNodeIndex(P,i,0,k)] = trinomial(P,i,k)
      *Bijk(i,k,P-i-k,xii[0],xii[2],xii[3]);

  // internal nodes
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j)
      for(int k = 1; k <= P-1-i-j; ++k)
        values[computeTetNodeIndex(P,i,j,k)] = quadnomial(P,i,j,k)
        *Bijkl(i,j,k,P-i-j-k,xii[0],xii[1],xii[2],xii[3]);

}

static void bezierTetGrads(int P, apf::Vector3 const& xi,
    apf::NewArray<apf::Vector3>& grads)
{
  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
  apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

  for(int i = 0; i < 4; ++i)
    grads[i] = gxii[i]*P*intpow(xii[i],P-1);

  int nE = P-1;

  int const (*tev)[2] = apf::tet_edge_verts;
  int const (*ttv)[3] = apf::tet_tri_verts;

  for(int a = 0; a < 6; ++a)
    for(int b = 0; b < nE; ++b) // edge nodes
      grads[4+a*nE+b] = gxii[tev[a][0]]*binomial(P,b+1)*(P-b-1)
                        *Bij(P-b-2,b+1,xii[tev[a][0]],xii[tev[a][1]])
                      + gxii[tev[a][1]]*binomial(P,b+1)*(b+1)
                        *Bij(P-b-1,b,xii[tev[a][0]],xii[tev[a][1]]);

  // face 0, l = 0
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j){
      int index = computeTetNodeIndex(P,i,j,P-i-j);
      grads[index].zero();
      int ijk[3] = {i,j,P-i-j};
      for(int b = 0; b < 3; ++b){
        grads[index] += gxii[ttv[0][b]]*ijk[b]
                       *Bijk(ijk[b % 3]-1,ijk[(b+1) % 3],ijk[(b+2) % 3],
                        xii[ttv[0][b % 3]],xii[ttv[0][(b+1) % 3]],
                        xii[ttv[0][(b+2) % 3]]);

      }
      grads[index] = grads[index]*trinomial(P,i,j);
    }
  // face 1, k = 0
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j){
      int index = computeTetNodeIndex(P,i,j,0);
      grads[index].zero();
      int ijk[3] = {i,j,P-i-j};
      for(int b = 0; b < 3; ++b){
        grads[index] += gxii[ttv[1][b]]*ijk[b]
                       *Bijk(ijk[b % 3]-1,ijk[(b+1) % 3],ijk[(b+2) % 3],
                        xii[ttv[1][b % 3]],xii[ttv[1][(b+1) % 3]],
                        xii[ttv[1][(b+2) % 3]]);

      }
      grads[index] = grads[index]*trinomial(P,i,j);
    }
  // face 2, i = 0
  for(int j = 1; j <= P-1; ++j)
    for(int k = 1; k <= P-1-j; ++k){
      int index = computeTetNodeIndex(P,0,j,k);
      grads[index].zero();
      int jkl[3] = {j,k,P-j-k};
      for(int b = 0; b < 3; ++b){
        grads[index] += gxii[ttv[2][b]]*jkl[b]
                       *Bijk(jkl[b % 3]-1,jkl[(b+1) % 3],jkl[(b+2) % 3],
                        xii[ttv[2][b % 3]],xii[ttv[2][(b+1) % 3]],
                        xii[ttv[2][(b+2) % 3]]);

      }
      grads[index] = grads[index]*trinomial(P,j,k);
    }
  // face 3, j = 0
  for(int i = 1; i <= P-1; ++i)
      for(int k = 1; k <= P-1-i; ++k){
        int index = computeTetNodeIndex(P,i,0,k);
        grads[index].zero();
        int ikl[3] = {i,k,P-i-k};
        for(int b = 0; b < 3; ++b){
          grads[index] += gxii[ttv[3][b]]*ikl[b]
                         *Bijk(ikl[b % 3]-1,ikl[(b+1) % 3],ikl[(b+2) % 3],
                          xii[ttv[3][b % 3]],xii[ttv[3][(b+1) % 3]],
                          xii[ttv[3][(b+2) % 3]]);

        }
        grads[index] = grads[index]*trinomial(P,i,k);
      }
  // internal nodes
  for(int i = 1; i <= P-1; ++i)
    for(int j = 1; j <= P-1-i; ++j)
      for(int k = 1; k <= P-1-i-j; ++k){
        int index = computeTetNodeIndex(P,i,j,k);
        grads[index].zero();
        int ijkl[4] = {i,j,k,P-i-j-k};
        for(int b = 0; b < 4; ++b){
          grads[index] += gxii[b]*ijkl[b]
             *Bijkl(ijkl[b % 4]-1,ijkl[(b+1) % 4],ijkl[(b+2) % 4],ijkl[(b+3) % 4],
                 xii[b],xii[(b+1) % 4],xii[(b+2) % 4],xii[(b+3) % 4]);

        }
        grads[index] = grads[index]*quadnomial(P,i,j,k);
      }
}

const bezierShape bezier[apf::Mesh::TYPES] =
{
  NULL,    //vertex
  bezierCurve,     //edge
  bezierTriangle,  //triangle
  NULL,    //quad
  bezierTet,       //tet
  NULL,    //hex
  NULL,    //prism
  NULL     //pyramid
};

const bezierShapeGrads bezierGrads[apf::Mesh::TYPES] =
{
  NULL,    //vertex
  bezierCurveGrads,     //edge
  bezierTriangleGrads,  //triangle
  NULL,    //quad
  bezierTetGrads,       //tet
  NULL,    //hex
  NULL,    //prism
  NULL     //pyramid
};

}
