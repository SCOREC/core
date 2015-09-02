/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvBezier.h"
#include "crvTables.h"
#include "crvQuality.h"

namespace crv {

static int maxAdaptiveIter = 5;

static double minThreshold = 0.05;
/* This work is based on the approach of Geometric Validity of high-order
 * lagrange finite elements, theory and practical guidance,
 * by George, Borouchaki, and Barral. (2014)
 *
 * The notation follows theirs, almost exactly.
 */
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
  double sum = 0.;
  int CD = trinomial(2*(d-1),I,J);
  for(int j1 = 0; j1 <= J; ++j1){
    int i1start = std::max(0,I+J-j1-(d-1));
    int i1end = std::min(I,d-1-j1);
    for(int i1 = i1start; i1 <= i1end; ++i1){
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
    int k2start = std::max(0,K-k1-(d-1));
    for (int k2 = k2start; k2 <= K-k1; ++k2){
      for (int j1 = 0; j1 <= J; ++j1){
        int j2start = std::max(0,J-j1-(d-1));
        for (int j2 = j2start; j2 <= J-j1; ++j2){
          int i1end = std::min(I,d-1-j1-k1);
          for (int i1 = 0; i1 <= i1end; ++i1){
            int i2start = std::max(0,I+J+K-i1-j1-k1-j2-k2-(d-1));
            int i2end = std::min(I-i1,d-1-j2-k2);
            for (int i2 = i2start; i2 <= i2end; ++i2){
              int i3 = I-i1-i2;
              int j3 = J-j1-j2;
              int k3 = K-k1-k2;
              sum += quadnomial(d-1,i1,j1,k1)*quadnomial(d-1,i2,j2,k2)
                  *quadnomial(d-1,i3,j3,k3)
                  *getTetPartialJacobianDet(nodes,d,i1,j1,k1,i2,j2,k2,i3,j3,k3);
            }
          }
        }
      }
    }
  }
  return sum*d*d*d/CD;
}

//static double calcMinTriJacDet(int P, apf::NewArray<apf::Vector3>& nodes)
//{
//  double minJ = 1e10;
//  double maxJ = -1e10;
//  for (int I = 0; I <= 2*(P-1); ++I){
//    for (int J = 0; J <= 2*(P-1)-I; ++J){
//      double nijk = Nijk(nodes,P,I,J);
//      minJ = std::min(minJ,nijk);
//      maxJ = std::max(maxJ,nijk);
//    }
//  }
//  return minJ/maxJ;
//}


struct JacobianData
{
  JacobianData(int P) : order(P), I(-1), J(-1), K(-1),
      minJ(1e10), maxJ(-1e10) {};
  int order;
  int I;
  int J;
  int K;
  double minJ;
  double maxJ;
};

static double calcMinTriJacDet(int P, apf::NewArray<double>& nodes,
    JacobianData& data)
{
  for (int I = 0; I <= P; ++I){
    for (int J = 0; J <= P-I; ++J){
      if(nodes[b2[P][I][J]] < data.minJ){
        data.minJ = nodes[b2[P][I][J]];
        data.I = I;
        data.J = J;
      }
      data.maxJ = std::max(data.maxJ,nodes[b2[P][I][J]]);
    }
  }
  return data.minJ;
}


/* nodes is (2(P-1)+1)(2(P-1)+2)/2 = P(2P-1)
 * except for P = 1, which has size 3 due to numbering convention used,
 * such that i=j=k=0 results in index 2
 */
static void getTriJacDetNodes(int P, apf::NewArray<apf::Vector3>& elemNodes,
    apf::NewArray<double>& nodes)
{
  for (int I = 0; I <= 2*(P-1); ++I){
    for (int J = 0; J <= 2*(P-1)-I; ++J){
      nodes[b2[2*(P-1)][I][J]] = Nijk(elemNodes,P,I,J);
    }
  }
}

//static int getTriEdge(int P, int I, int J)
//{
//  unsigned index = b2[P][I][J];
//  assert(index > 2); // should never get a negative jacobian at the corners
//  if (index > 3*P) // interior
//    return -1;
//  assert( (index-3)/(P-1) < 3);
//  return (index-3)/(P-1);
//}

static double getMinTriJacDet(int P, apf::NewArray<double>& nodes, int& iter,
    JacobianData& data)
{

  double minJ = calcMinTriJacDet(P,nodes,data);
  // stop if this is true
  if(iter >= maxAdaptiveIter || data.minJ > minThreshold){
    return data.minJ;
  } else {
    iter++;
    // subdivide, and continue trying
    apf::NewArray<double> subNodes[4];
    int n = (P+1)*(P+2)/2;
    subNodes[0].allocate(n);
    subNodes[1].allocate(n);
    subNodes[2].allocate(n);
    subNodes[3].allocate(n);

    subdivideBezierTriangleJacobianDet(P,nodes,subNodes);
    minJ = getMinTriJacDet(P,subNodes[0],iter,data);
    minJ = std::min(getMinTriJacDet(P,subNodes[1],iter,data),minJ);
    minJ = std::min(getMinTriJacDet(P,subNodes[2],iter,data),minJ);
    minJ = std::min(getMinTriJacDet(P,subNodes[3],iter,data),minJ);
    return minJ;
  }
}

static bool checkTriValidityAtNodeXi(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* me = apf::createMeshElement(m,e);
  apf::FieldShape* fs = m->getShape();
  apf::Matrix3x3 J;
  // First, just check at node xi
  apf::Downward down;
  apf::Vector3 xi, exi;
  for (int d = 0; d <= 2; ++d){
    int nDown = m->getDownward(e,d,down);
    for (int j = 0; j < nDown; ++j){
      int bt = m->getType(down[j]);
      for (int i = 0; i < fs->countNodesOn(bt); ++i){
        fs->getNodeXi(bt,i,xi);
        exi = apf::boundaryToElementXi(m,down[j],e,xi);
        apf::getJacobian(me,exi,J);
        if(J[0][0]*J[1][1]-J[1][0]*J[0][1] < 0.) return false;
      }
    }
  }
  apf::destroyMeshElement(me);
  return true;
}

static bool checkTriValidity(apf::Mesh* m, apf::MeshEntity* e)
{

  int P = m->getShape()->getOrder();
  bool valid = checkTriValidityAtNodeXi(m,e);
  if(!valid) return false;
  if(P == 1) return true;
  // if it is positive, then keep going
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);

  apf::NewArray<apf::Vector3> elemNodes;
  apf::getVectorNodes(elem,elemNodes);
  apf::destroyElement(elem);
  apf::NewArray<double> nodes(P*(2*P-1));
  getTriJacDetNodes(P,elemNodes,nodes);

  int iter = 0;
  JacobianData data(P);
  double minJ = getMinTriJacDet(2*(P-1),nodes,iter,data);
  printf("minj %f %d\n",minJ,iter);
  if(minJ < minThreshold || iter == maxAdaptiveIter)
    return false;
  return true;
}

static bool checkTetValidity(apf::Mesh* /* m*/, apf::MeshEntity* /*e*/)
{
 return false;
}

bool checkValidity(apf::Mesh* m, apf::MeshEntity* e)
{
  if(m->getDimension() == 2)
    return checkTriValidity(m,e);
  else
    return checkTetValidity(m,e);
}

double computeAlternateTriJacobianDet(apf::Mesh* m,
    apf::MeshEntity* e, apf::Vector3& xi)
{
  double detJ = 0.;
  int P = m->getShape()->getOrder();
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);

  for (int I = 0; I <= 2*(P-1); ++I){
    for (int J = 0; J <= 2*(P-1)-I; ++J){
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
  for (int I = 0; I <= 3*(P-1); ++I){
    for (int J = 0; J <= 3*(P-1)-I; ++J){
      for (int K = 0; K <= 3*(P-1)-I-J; ++K){
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
