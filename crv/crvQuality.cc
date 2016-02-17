/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvTables.h"
#include "crvQuality.h"

namespace crv {

static int maxAdaptiveIter = 5;

static int maxElevationLevel = 19;

static double minAcceptable = 0.0;

static double convergenceTolerance = 0.01;

/* This work is based on the approach of Geometric Validity of high-order
 * lagrange finite elements, theory and practical guidance,
 * by George, Borouchaki, and Barral. (2014)
 *
 * The notation follows theirs, almost exactly.
 */
static double getTriPartialJacobianDet(apf::NewArray<apf::Vector3>& nodes,
    int P, int i1, int j1, int i2, int j2)
{
  int p00 = getTriNodeIndex(P,i1+1,j1);
  int p01 = getTriNodeIndex(P,i1,j1+1);
  int p10 = getTriNodeIndex(P,i2+1,j2);
  int p11 = getTriNodeIndex(P,i2,j2);
  return apf::cross(nodes[p01]-nodes[p00],
      nodes[p11]-nodes[p10])[2];
}

static double getTetPartialJacobianDet(apf::NewArray<apf::Vector3>& nodes,
    int P, int i1, int j1, int k1, int i2, int j2, int k2,
    int i3, int j3, int k3)
{
  int p00 = getTetNodeIndex(P,i1+1,j1,k1);
  int p01 = getTetNodeIndex(P,i1,j1+1,k1);
  int p10 = getTetNodeIndex(P,i2+1,j2,k2);
  int p11 = getTetNodeIndex(P,i2,j2,k2+1);
  int p20 = getTetNodeIndex(P,i3+1,j3,k3);
  int p21 = getTetNodeIndex(P,i3,j3,k3);
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

static double calcMinJacDet(int n, apf::NewArray<double>& nodes)
{
  double minJ = 1e10;
  for (int i = 0; i < n; ++i)
    minJ = std::min(minJ,nodes[i]);
  return minJ;
}

static double calcMaxJacDet(int n, apf::NewArray<double>& nodes)
{
  double maxJ = -1e10;
  for (int i = 0; i < n; ++i)
    maxJ = std::max(maxJ,nodes[i]);
  return maxJ;
}
/* nodes is (2(P-1)+1)(2(P-1)+2)/2 = P(2P-1)
 * except for P = 1, which has size 3 due to numbering convention used,
 * such that i=j=k=0 results in index 2
 */
static void getTriJacDetNodes(int P, apf::NewArray<apf::Vector3>& elemNodes,
    apf::NewArray<double>& nodes)
{
  for (int I = 0; I <= 2*(P-1); ++I)
    for (int J = 0; J <= 2*(P-1)-I; ++J)
      nodes[getTriNodeIndex(2*(P-1),I,J)] = Nijk(elemNodes,P,I,J);
}
static void getTetJacDetNodes(int P, apf::NewArray<apf::Vector3>& elemNodes,
    apf::NewArray<double>& nodes)
{
  for (int I = 0; I <= 3*(P-1); ++I)
    for (int J = 0; J <= 3*(P-1)-I; ++J)
      for (int K = 0; K <= 3*(P-1)-I-J; ++K)
         nodes[getTetNodeIndex(3*(P-1),I,J,K)] = Nijkl(elemNodes,P,I,J,K);
}

/*
 * This is the elevation version of this algorithm
 * There is no recursion needed, at least not yet
 *
 */
static void getJacDetByElevation(int type, int P,
    apf::NewArray<double>& nodes, double& minJ, double& maxJ)
{
  /*
   * as a convergence check, use the max dist between points
   */
  double maxDist[2] = {0.,1e10};

  int n = getNumControlPoints(type,P);
  minJ = calcMinJacDet(n,nodes);
  maxJ = calcMaxJacDet(n,nodes);
  maxDist[0] = nodes[1]-nodes[0];
  for (int j = 1; j < n-1; ++j)
    maxDist[0] = std::max(maxDist[0],nodes[j+1]-nodes[j]);

  // declare these two arrays, never need to reallocate
  apf::NewArray<double> elevatedNodes[2];
  elevatedNodes[0].allocate(getNumControlPoints(type,maxElevationLevel));
  elevatedNodes[1].allocate(getNumControlPoints(type,maxElevationLevel));

  // copy for the start
  for(int i = 0; i < n; ++i)
    elevatedNodes[0][i] = nodes[i];

  int i = 0;
  while(P+i < maxElevationLevel && minJ/maxJ < minAcceptable
    && std::fabs(maxDist[(i+1) % 2] - maxDist[i % 2]) > convergenceTolerance){
    // use the modulus to alternate between them,
    // never needing to increase storage
    // for now, only elevate by 1
    elevateBezierJacobianDet(type,P+i,1,
        elevatedNodes[i % 2],
        elevatedNodes[(i+1) % 2]);

    int ni = getNumControlPoints(type,P+i);
    maxDist[(i+1) % 2] = elevatedNodes[(i+1) % 2][1]-elevatedNodes[(i+1) % 2][0];
    for (int j = 1; j < ni-1; ++j)
      maxDist[(i+1) % 2] = std::max(elevatedNodes[(i+1) % 2][j+1]
                         - elevatedNodes[(i+1) % 2][j],maxDist[(i+1) % 2]);

    minJ = calcMinJacDet(ni,elevatedNodes[(i+1) % 2]);
    maxJ = calcMaxJacDet(ni,elevatedNodes[(i+1) % 2]);

    ++i;
  }
}

/*
 * This is the subdivision version, with recursion
 *
 */
static int numSplits[apf::Mesh::TYPES] =
  {0,2,4,0,8,0,0,0};

static void getJacDetBySubdivision(int type, int P,
    int iter, apf::NewArray<double>& nodes,
    double& minJ, double& maxJ, bool& done)
{
  int n = getNumControlPoints(type,P);
  double change = minJ;
  if(!done){
    minJ = calcMinJacDet(n,nodes);
    maxJ = calcMaxJacDet(n,nodes);
    change = minJ - change;
  }

  if(!done && iter < maxAdaptiveIter && minJ/maxJ < minAcceptable
      && std::fabs(change) > convergenceTolerance){
    iter++;
    apf::NewArray<double> subNodes[8];
    for (int i = 0; i < numSplits[type]; ++i)
      subNodes[i].allocate(n);

    subdivideBezierJacobianDet[type](P,nodes,subNodes);

    apf::NewArray<double> newMinJ(numSplits[type]);
    apf::NewArray<double> newMaxJ(numSplits[type]);
    for (int i = 0; i < numSplits[type]; ++i){
      newMinJ[i] = 1e10;
      newMaxJ[i] = -1e10;
    }

    for (int i = 0; i < numSplits[type]; ++i)
      getJacDetBySubdivision(type,P,iter,subNodes[i],
          newMinJ[i],newMaxJ[i],done);

    minJ = newMinJ[0];
    maxJ = newMaxJ[0];
    for (int i = 1; i < numSplits[type]; ++i){
      minJ = std::min(newMinJ[i],minJ);
      maxJ = std::max(newMaxJ[i],maxJ);
    }
  } else if (minJ/maxJ < minAcceptable){
    done = true;
  }
}

static void getJacDetBySubdivisionMatrices(int type, int P,
    int iter, apf::NewArray<double>& c,apf::NewArray<double>& nodes,
    double& minJ, double& maxJ, bool& done, bool& quality)
{
  int n = getNumControlPoints(type,P);
  double change = minJ;
  if(!done){
    minJ = calcMinJacDet(n,nodes);
    maxJ = calcMaxJacDet(n,nodes);
    change = minJ - change;
  }

  if(!done && iter < maxAdaptiveIter && (quality || minJ/maxJ < minAcceptable)
      && std::fabs(change) > convergenceTolerance){

    iter++;
    apf::NewArray<double> subNodes[8];
    for (int i = 0; i < numSplits[type]; ++i)
      subNodes[i].allocate(n);

    subdivideBezierEntityJacobianDet(P,type,c,nodes,subNodes);

    apf::NewArray<double> newMinJ(numSplits[type]);
    apf::NewArray<double> newMaxJ(numSplits[type]);
    for (int i = 0; i < numSplits[type]; ++i){
      newMinJ[i] = 1e10;
      newMaxJ[i] = -1e10;
    }

    for (int i = 0; i < numSplits[type]; ++i)
      getJacDetBySubdivisionMatrices(type,P,iter,c,subNodes[i],
          newMinJ[i],newMaxJ[i],done,quality);

    minJ = newMinJ[0];
    maxJ = newMaxJ[0];
    for (int i = 1; i < numSplits[type]; ++i){
      minJ = std::min(newMinJ[i],minJ);
      maxJ = std::max(newMaxJ[i],maxJ);
    }
  } else if (minJ/maxJ < minAcceptable){
    done = true;
  }
}
/*
 * Distinctly only works for 2D planar meshes in x-y, used as a test for
 * the algorithm used in 3D.
 */
static int checkTriValidityAtNodeXi(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* me = apf::createMeshElement(m,e);
  apf::FieldShape* fs = m->getShape();
  apf::Matrix3x3 J;
  // First, just check at node xi
  apf::Downward down;
  apf::Vector3 xi, exi;
  int qualityTag = 0;
  for (int d = 0; d <= 1; ++d){
    int nDown = m->getDownward(e,d,down);
    for (int j = 0; j < nDown; ++j){
      int bt = m->getType(down[j]);
      for (int i = 0; i < fs->countNodesOn(bt); ++i){
        fs->getNodeXi(bt,i,xi);
        exi = apf::boundaryToElementXi(m,down[j],e,xi);
        apf::getJacobian(me,exi,J);
        if(J[0][0]*J[1][1]-J[1][0]*J[0][1] < minAcceptable){
          qualityTag = 6*d+2+j;
          break;
        }
      }
    }
  }

  for (int i = 0; i < fs->countNodesOn(apf::Mesh::TRIANGLE); ++i){
    fs->getNodeXi(apf::Mesh::TRIANGLE,i,xi);
    apf::getJacobian(me,xi,J);
    if(J[0][0]*J[1][1]-J[1][0]*J[0][1] < minAcceptable){
      qualityTag = 14;
      break;
    }
  }

  apf::destroyMeshElement(me);
  return qualityTag;
}

static int checkTriValidity(apf::Mesh* m, apf::MeshEntity* e,
    int algorithm)
{
  int P = m->getShape()->getOrder();
  if (algorithm >= 2){
    int qualityTag = checkTriValidityAtNodeXi(m,e);
    if(qualityTag > 0) return qualityTag;
  }

  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> elemNodes;
  apf::getVectorNodes(elem,elemNodes);
  // if we are blended, we need to create a full representation
  int blendingOrder = getBlendingOrder(apf::Mesh::TRIANGLE);
  if (blendingOrder > 0 && getNumInternalControlPoints(apf::Mesh::TRIANGLE,P)) {
    getFullRepFromBlended(m,apf::Mesh::TRIANGLE,elemNodes);
  }

  apf::destroyElement(elem);
  apf::NewArray<double> nodes(P*(2*P-1));
  getTriJacDetNodes(P,elemNodes,nodes);

  // if we haven't checked at nodeXi, need to check verts
  if(algorithm < 2){
    apf::Downward verts;
    m->getDownward(e,0,verts);
    for (int i = 0; i < 3; ++i){
      if(nodes[i] < minAcceptable){
        return i+2;
      }
    }
  }

  apf::NewArray<double> subCoefficients;
  if(algorithm == 4){
    getBezierJacobianDetSubdivisionCoefficients(2*(P-1),
        apf::Mesh::EDGE,subCoefficients);
  }

  apf::MeshEntity* edges[3];
  m->getDownward(e,1,edges);
  double minJ = 0, maxJ = 0;
  // Vertices will already be flagged in the first check
  for (int edge = 0; edge < 3; ++edge){
    for (int i = 0; i < 2*(P-1)-1; ++i){
      if (nodes[3+edge*(2*(P-1)-1)+i] < minAcceptable){
        minJ = -1e10;
        apf::NewArray<double> edgeNodes(2*(P-1)+1);
        if(algorithm != 4){
          edgeNodes[0] = nodes[apf::tri_edge_verts[edge][0]];
          edgeNodes[2*(P-1)] = nodes[apf::tri_edge_verts[edge][1]];
          for (int j = 0; j < 2*(P-1)-1; ++j)
            edgeNodes[j+1] = nodes[3+edge*(2*(P-1)-1)+j];
          if(algorithm % 2 == 1){
            getJacDetByElevation(apf::Mesh::EDGE,2*(P-1),edgeNodes,minJ,maxJ);
          } else {
            // allows recursion stop on first "conclusive" invalidity
            bool done = false;
            getJacDetBySubdivision(apf::Mesh::EDGE,2*(P-1),
                0,edgeNodes,minJ,maxJ,done);
          }
        } else {
          edgeNodes[0] = nodes[apf::tri_edge_verts[edge][0]];
          edgeNodes[1] = nodes[apf::tri_edge_verts[edge][1]];
          for (int j = 0; j < 2*(P-1)-1; ++j)
            edgeNodes[j+2] = nodes[3+edge*(2*(P-1)-1)+j];
          bool done = false;
          bool quality = false;
          getJacDetBySubdivisionMatrices(apf::Mesh::EDGE,2*(P-1),
              0,subCoefficients,edgeNodes,minJ,maxJ,done,quality);
        }
        if(minJ < minAcceptable){
          return 8+edge;
        }
      }
    }
  }

  if(algorithm == 4){
    getBezierJacobianDetSubdivisionCoefficients(2*(P-1),
        apf::Mesh::TRIANGLE,subCoefficients);
  }
  bool done = false;
  for (int i = 0; i < (2*P-3)*(2*P-4)/2; ++i){
    if (nodes[6*(P-1)+i] < minAcceptable){
      minJ = -1e10;
      if(algorithm % 2 == 1)
        getJacDetByElevation(apf::Mesh::TRIANGLE,2*(P-1),nodes,minJ,maxJ);
      else if(algorithm == 4){
        bool quality = false;
        getJacDetBySubdivisionMatrices(apf::Mesh::TRIANGLE,2*(P-1),
            0,subCoefficients,nodes,minJ,maxJ,done,quality);
      } else {
        getJacDetBySubdivision(apf::Mesh::TRIANGLE,2*(P-1),
            0,nodes,minJ,maxJ,done);
      }
      if(minJ < minAcceptable){
        return 12;
      }
    }
  }
  return 1;
}

static int checkTetValidityAtNodeXi(apf::Mesh* m, apf::MeshEntity* e)
{
  apf::MeshElement* me = apf::createMeshElement(m,e);
  apf::FieldShape* fs = m->getShape();
  apf::Matrix3x3 J;
  // First, just check at node xi
  apf::Downward down;
  apf::Vector3 xi, exi;
  int qualityTag = 0;
  for (int d = 0; d <= 2; ++d){
    int nDown = m->getDownward(e,d,down);
    for (int j = 0; j < nDown; ++j){
      int bt = m->getType(down[j]);
      for (int i = 0; i < fs->countNodesOn(bt); ++i){
        fs->getNodeXi(bt,i,xi);
        exi = apf::boundaryToElementXi(m,down[j],e,xi);
        apf::getJacobian(me,exi,J);
        if(apf::getJacobianDeterminant(J,3) < minAcceptable){
          qualityTag = 6*d+2+j;
          break;
        }
      }
    }
  }
  for (int i = 0; i < fs->countNodesOn(apf::Mesh::TET); ++i){
    fs->getNodeXi(apf::Mesh::TET,i,xi);
    apf::getJacobian(me,xi,J);
    if(apf::getJacobianDeterminant(J,3) < minAcceptable){
      qualityTag = 20;
      break;
    }
  }

  apf::destroyMeshElement(me);
  return qualityTag;
}

static int checkTetValidity(apf::Mesh* m, apf::MeshEntity* e,
    int algorithm)
{
  int P = m->getShape()->getOrder();
  if (algorithm >= 2){
    int qualityTag = checkTetValidityAtNodeXi(m,e);
    if(qualityTag > 0) return qualityTag;
  }

  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> elemNodes;
  apf::getVectorNodes(elem,elemNodes);
  int blendingOrder = getBlendingOrder(apf::Mesh::TET);
  if (blendingOrder > 0 && getNumInternalControlPoints(apf::Mesh::TET,P)) {
    getFullRepFromBlended(m,apf::Mesh::TET,elemNodes);
  }

  apf::destroyElement(elem);
  // 9*P*P*(P-1)/2+P = (3(P-1)+1)(3(P-1)+2)(3(P-1)+3)/6
  apf::NewArray<double> nodes(9*P*P*(P-1)/2+P);
  getTetJacDetNodes(P,elemNodes,nodes);

  if(algorithm < 2){
    apf::Downward verts;
    m->getDownward(e,0,verts);
    for (int i = 0; i < 4; ++i){
      if(nodes[i] < minAcceptable){
        return 2+i;
      }
    }
  }

  apf::NewArray<double> subCoefficients;
  if(algorithm == 4){
    getBezierJacobianDetSubdivisionCoefficients(3*(P-1),
        apf::Mesh::EDGE,subCoefficients);
  }

  apf::MeshEntity* edges[6];
  m->getDownward(e,1,edges);
  double minJ = 0, maxJ = 0;
  // Vertices will already be flagged in the first check
  for (int edge = 0; edge < 6; ++edge){
    for (int i = 0; i < 3*(P-1)-1; ++i){
      if (nodes[4+edge*(3*(P-1)-1)+i] < minAcceptable){
        minJ = -1e10;
        apf::NewArray<double> edgeNodes(3*(P-1)+1);

        if(algorithm != 4){
          edgeNodes[0] = nodes[apf::tet_edge_verts[edge][0]];
          edgeNodes[3*(P-1)] = nodes[apf::tet_edge_verts[edge][1]];
          for (int j = 0; j < 3*(P-1)-1; ++j)
            edgeNodes[j+1] = nodes[4+edge*(3*(P-1)-1)+j];
          // allows recursion stop on first "conclusive" invalidity
          if(algorithm % 2 == 1)
            getJacDetByElevation(apf::Mesh::EDGE,3*(P-1),edgeNodes,minJ,maxJ);
          else {
            bool done = false;
            getJacDetBySubdivision(apf::Mesh::EDGE,3*(P-1),
                0,edgeNodes,minJ,maxJ,done);
          }
        } else {
          edgeNodes[0] = nodes[apf::tet_edge_verts[edge][0]];
          edgeNodes[1] = nodes[apf::tet_edge_verts[edge][1]];
          for (int j = 0; j < 3*(P-1)-1; ++j)
            edgeNodes[j+2] = nodes[4+edge*(3*(P-1)-1)+j];

          bool done = false;
          bool quality = false;
          getJacDetBySubdivisionMatrices(apf::Mesh::EDGE,3*(P-1),
              0,subCoefficients,edgeNodes,minJ,maxJ,done,quality);
        }
        if(minJ < minAcceptable){
          return 8+edge;
        }
      }
    }
  }
  if(algorithm == 4){
    getBezierJacobianDetSubdivisionCoefficients(3*(P-1),
        apf::Mesh::TRIANGLE,subCoefficients);
  }
  apf::MeshEntity* faces[4];
  m->getDownward(e,2,faces);
  for (int face = 0; face < 4; ++face){
    double minJ = -1e10;
    for (int i = 0; i < (3*P-4)*(3*P-5)/2; ++i){
      if (nodes[18*P-20+face*(3*P-4)*(3*P-5)/2+i] < minAcceptable){
        minJ = -1e10;
        apf::NewArray<double> triNodes((3*P-2)*(3*P-1)/2);
        getTriDetJacNodesFromTetDetJacNodes(face,3*(P-1),nodes,triNodes);
        if(algorithm == 4){
          bool done = false;
          bool quality = false;
          getJacDetBySubdivisionMatrices(apf::Mesh::TRIANGLE,3*(P-1),
              0,subCoefficients,triNodes,minJ,maxJ,done,quality);
        } else if(algorithm % 2 == 1)
          getJacDetByElevation(apf::Mesh::TRIANGLE,3*(P-1),triNodes,minJ,maxJ);
        else {
          bool done = false;
          getJacDetBySubdivision(apf::Mesh::TRIANGLE,3*(P-1),
              0,triNodes,minJ,maxJ,done);
        }
        if(minJ < minAcceptable){
          return 14+face;
        }
      }
    }
  }

  getBezierJacobianDetSubdivisionCoefficients(3*(P-1),
      apf::Mesh::TET,subCoefficients);

  for (int i = 0; i < (3*P-4)*(3*P-5)*(3*P-6)/6; ++i){
    if (nodes[18*P*P-36*P+20+i] < minAcceptable){
      minJ = -1e10;
      if(algorithm % 2 == 0){
        bool done = false;
        bool quality = false;
        getJacDetBySubdivisionMatrices(apf::Mesh::TET,3*(P-1),
            0,subCoefficients,nodes,minJ,maxJ,done,quality);
      } else {
        getJacDetByElevation(apf::Mesh::TET,3*(P-1),nodes,minJ,maxJ);
      }
      if(minJ < minAcceptable){
        return 20;
      }
    }
  }
  return 1;
}

double computeTriJacobianDetFromBezierFormulation(apf::Mesh* m,
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

double computeTetJacobianDetFromBezierFormulation(apf::Mesh* m,
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

double getQuality(int type, int P, apf::NewArray<apf::Vector3>& elemNodes)
{
  int typeDim = apf::Mesh::typeDimension[type];

  apf::NewArray<double> nodes(getNumControlPoints(type,typeDim*(P-1)));
  if(type == apf::Mesh::TRIANGLE){
    getTriJacDetNodes(P,elemNodes,nodes);
  }
  if(type == apf::Mesh::TET){
    getTetJacDetNodes(P,elemNodes,nodes);
  }
  bool done = false;
  double minJ = -1e10, maxJ = -1e10;

  double oldAcceptable = minAcceptable;
  int oldIter = maxAdaptiveIter;
  maxAdaptiveIter = 1;
  minAcceptable = -1e10;
  bool quality = true;

  apf::NewArray<double> subCoefficients;
  getBezierJacobianDetSubdivisionCoefficients(typeDim*(P-1),
      type,subCoefficients);
  getJacDetBySubdivisionMatrices(type,typeDim*(P-1),
      0,subCoefficients,nodes,minJ,maxJ,done,quality);
  done = false;
  minAcceptable = oldAcceptable;
  maxAdaptiveIter = oldIter;
  if(std::fabs(maxJ) > 1e-8)
    return minJ/maxJ;
  else return minJ;
}

double getQuality(apf::Mesh* m, apf::MeshEntity* e)
{
  int P = m->getShape()->getOrder();
  int type = m->getType(e);

  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> elemNodes;
  apf::getVectorNodes(elem,elemNodes);
  int blendingOrder = getBlendingOrder(type);
  if (blendingOrder > 0 && getNumInternalControlPoints(type,P)) {
    getFullRepFromBlended(m,type,elemNodes);
  }
  apf::destroyElement(elem);

  return getQuality(type,P,elemNodes);
}
typedef int (*bezierValidity)(apf::Mesh* m,
    apf::MeshEntity* e, int algorithm);

const bezierValidity checkBezierValidity[apf::Mesh::TYPES] =
{
  NULL,    //vertex
  NULL,    //edge
  checkTriValidity,  //triangle
  NULL,    //quad
  checkTetValidity,  //tet
  NULL,    //hex
  NULL,    //prism
  NULL     //pyramid
};

int checkValidity(apf::Mesh* m, apf::MeshEntity* e,
    int algorithm)
{
  return checkBezierValidity[m->getType(e)](m,e,algorithm);
}

}
