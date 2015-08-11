/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvBezier.h"
#include "crvTables.h"

namespace crv {

int getTriPointIndex(int P, int i, int j)
{
  return crv::b2[P][j*(P+1)+i-j*(j-1)/2];
}

int getTetPointIndex(int P, int i, int j, int k)
{
  return crv::b3[P][i][j][k];
}

static double getTriPartialJacobianDet(apf::NewArray<apf::Vector3>& nodes,
    int P, int i1, int j1, int i2, int j2)
{
  int p00 = getTriPointIndex(P,i1+1,j1);
  int p01 = getTriPointIndex(P,i1,j1+1);
  int p10 = getTriPointIndex(P,i2+1,j2);
  int p11 = getTriPointIndex(P,i2,j2);
  return apf::cross(nodes[p01]-nodes[p00],
      nodes[p11]-nodes[p10])[2];
}

static double getTetPartialJacobianDet(apf::NewArray<apf::Vector3>& nodes,
    int P, int i1, int j1, int k1, int i2, int j2, int k2,
    int i3, int j3, int k3)
{
  int p00 = getTetPointIndex(P,i1+1,j1,k1);
  int p01 = getTetPointIndex(P,i1,j1+1,k1);
  int p10 = getTetPointIndex(P,i2+1,j2,k2);
  int p11 = getTetPointIndex(P,i2,j2,k2+1);
  int p20 = getTetPointIndex(P,i3+1,j3,k3);
  int p21 = getTetPointIndex(P,i3,j3,k3);
  return apf::cross(nodes[p01]-nodes[p00],nodes[p11]-nodes[p10])
    *(nodes[p21]-nodes[p20]);
}

static double Nijk(apf::NewArray<apf::Vector3>& nodes,
    int d, int I, int J)
{
  int CD = trinomial(2*(d-1),I,J);
  double sum = 0.;
  for(int j1 = 0; j1 <= J; ++j1){
    for(int i1 = 0; i1 <= I; ++i1){
      if(i1 > d-1 || j1 > d-1 || I-i1 > d-1 || J-j1 > d-1 ||
         i1+j1 > d-1 || I-i1 + J-j1 > d-1)
        continue;
      sum += trinomial(d-1,i1,j1)*trinomial(d-1,I-i1,J-j1)
          *getTriPartialJacobianDet(nodes,d,i1,j1,I-i1,J-j1);
    }
  }
  return sum*d*d/CD;
}

static double Nijkl(apf::NewArray<apf::Vector3>& nodes,
    int d, int I, int J, int K)
{
  double sum = 0.;
  int CD = quadnomial(3*(d-1),I,J,K);

  for(int k1 = 0; k1 <= K; ++k1){
    for(int k2 = 0; k2 <= K-k1; ++k2){
      for(int j1 = 0; j1 <= J; ++j1){
        for(int j2 = 0; j2 <= J-j1; ++j2){

          for(int i1 = 0; i1 <= I; ++i1){
            for(int i2 = 0; i2 <= I-i1; ++i2){
              int i3 = I - i1 - i2;
              int j3 = J - j1 - j2;
              int k3 = K - k1 - k2;

              int l1 = d-1 - i1 - j1 - k1;
              int l2 = d-1 - i2 - j2 - k2;
              int l3 = d-1 - i3 - j3 - k3;

              if(i1 > d-1 || j1 > d-1 || k1 > d-1 ||
                  i2 > d-1 || j2 > d-1 || k2 > d-1 ||
                  i3 > d-1 || j3 > d-1 || k3 > d-1 ||
                  l1 > d-1 || l2 > d-1 || l3 > d-1 ||
                  l1 < 0 || l2 < 0 || l3 < 0)
                continue;
              sum += quadnomial(d-1,i1,j1,k1)*quadnomial(d-1,i2,j2,k2)*quadnomial(d-1,i3,j3,k3)
                  *getTetPartialJacobianDet(nodes,d,i1,j1,k1,i2,j2,k2,i3,j3,k3);
            }
          }
        }
      }
    }
  }
  return sum*d*d*d/CD;
}

double computeAlternateTriJacobianDet(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi)
{
  double detJ = 0.;
  int P = m->getShape()->getOrder();
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);

  for(int I = 0; I <= 2*(P-1); ++I){
    for(int J = 0; J <= 2*(P-1)-I; ++J){
      detJ += trinomial(2*(P-1),I,J)
          *Bijk(I,J,2*(P-1)-I-J,1.-xi[0]-xi[1],xi[0],xi[1])
          *Nijk(nodes,P,I,J);
    }
  }
  apf::destroyElement(elem);
  return detJ;
}

double computeAlternateTetJacobianDet(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi)
{
  int P = m->getShape()->getOrder();
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);

  double detJ = 0.;
  for(int I = 0; I <= 3*(P-1); ++I){
    for(int J = 0; J <= 3*(P-1)-I; ++J){
      for(int K = 0; K <= 3*(P-1)-I-J; ++K){
        detJ += quadnomial(3*(P-1),I,J,K)*Bijkl(I,J,K,3*(P-1)-I-J-K,
            1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2])
            *Nijkl(nodes,P,I,J,K);
      }
    }
  }
  apf::destroyElement(elem);
  return detJ;
}

}
