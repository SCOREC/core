#include <apfMatrix.h>
#include <apfMDS.h>
#include <apfNumbering.h>
#include <apfIntegrate.h>
#include <gmi.h>
#include "crvOptimizations.h"
#include "crv.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include "crvDBG.h"
#include <iostream>

/* static int global_counter = 0; */
/* static apf::MeshEntity* tetra[100]; */
/* static int number = 0; */



static void printInvalidities(apf::Mesh2* m, apf::MeshEntity* e[99], apf::MeshEntity* edge, int nat)
{

  return;
  apf::Numbering* n = m->findNumbering("debug_num_edge");
  PCU_ALWAYS_ASSERT(n);
  int num = apf::getNumber(n, edge, 0, 0);
  int tag = m->getModelTag(m->toModel(edge));
  std::cout<<"at edge "<< num <<" tag: "<<tag<< std::endl;

  for (int i = 0; i < nat; i++) {
    std::vector<int> ai = crv::getAllInvalidities(m, e[i]);
    for (std::size_t j = 0; j < ai.size(); j++) {
      printf("%d ", ai[j]);
    }
    printf("\n");
  }
}

static double computeFValNIJKL(apf::Mesh2* m, apf::MeshEntity* e, ma::SizeField* s = 0)
{
  PCU_ALWAYS_ASSERT_VERBOSE(s == 0, "Not implemented for non-zero sizefield!");

  int d = m->getDimension();
  int P = m->getShape()->getOrder();

  apf::NewArray<apf::Vector3> nodes;
  apf::Element* el = apf::createElement(m->getCoordinateField(), e);
  apf::getVectorNodes(el, nodes);
  apf::destroyElement(el);

  double volm = getLinearVolPhys(m, e);

  int weight = 1;
  double sumf = 0;
  if (d == 3) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
  for (int K = 0; K <= d*(P-1); K++) {
    for (int L = 0; L <= d*(P-1); L++) {
      if ((I == J && J == K && I == 0) ||
	(J == K && K == L && J == 0) ||
	(I == K && K == L && I == 0) ||
	(I == J && J == L && I == 0))
        weight = 4;
      else if ((I == J && I == 0) ||
	     (I == K && I == 0) ||
	     (I == L && I == 0) ||
	     (J == K && J == 0) ||
	     (J == L && J == 0) ||
	     (K == L && K == 0))
        weight = 2;
      else
        weight = 1;
      if (I + J + K + L == d*(P-1)) {
        double f = crv::Nijkl(nodes,P,I,J,K)/(6.0*volm) - 1.0;
        sumf = sumf + weight*f*f;
      }
    }
  }
      }
    }
  }

  if (d == 2) {
    for (int I = 0; I <= d*(P-1); I++) {
      for (int J = 0; J <= d*(P-1); J++) {
  for (int K = 0; K <= d*(P-1); K++) {
    if ((I == J && I == 0) ||
        (J == K && J == 0) ||
        (I == K && I == 0))
      weight = 2;
    else
      weight = 1;
    if (I + J + K == d*(P-1)) {
      double f = crv::Nijk(nodes,P,I,J)/(4.0*volm) - 1.0;
      sumf = sumf + weight*f*f;
    }
  }
      }
    }
  }
  return sumf;
}

static double computeFValDetJ(apf::Mesh2* m, apf::MeshEntity* e, ma::SizeField* s)
{
  int order = m->getShape()->getOrder();
  apf::MeshElement* me = apf::createMeshElement(m, e);
  apf::Matrix3x3 J;
  apf::Matrix3x3 T;
  apf::Matrix3x3 Jm;

  double jDet, sum = 0.;
  for (int i = 0; i < apf::countIntPoints(me, order) ; i++) {
    apf::Vector3 qp;
    double w = apf::getIntWeight(me, order, i);
    apf::getIntPoint(me, order, i, qp);

    apf::getJacobian(me, qp, J);
    s->getTransform(me, qp, T);
    Jm = J*T; // Jacobian in metric space
    jDet = apf::getDeterminant(Jm);
    sum += w * (jDet - 1.) * (jDet - 1.);
  }
  return sum;
}


namespace crv{

void CrvInternalEdgeOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvInternalEdgeOptim :: setTol(double tolerance)
{
  tol = tolerance;
}

bool CrvInternalEdgeOptim :: run(int &invaliditySize)
{
  std::vector<int> sizeHolder;
  apf::MeshEntity* adj_array[99];
  apf::Adjacent adj;
  mesh->getAdjacent(edge, 3, adj);
  int thisTetSize = 0;

  for (int i = 0; i < adj.getSize(); i++) {
    adj_array[i] = adj[i];
    std::vector<int> ai = crv::getAllInvalidities(mesh, adj[i]);
    if (adj[i] == tet) thisTetSize = ai.size();
    sizeHolder.push_back(ai.size());

  }

  //std::vector<int> ai = crv::getAllInvalidities(mesh, tet);
  //makeMultipleEntityMesh(mesh, adj_array, edge, "before_cavity_of_edge_", adj.getSize());
  //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "before_cavity_indv_tet_of_edge_", adj.getSize());
  /* printTetNumber(mesh, tet); */
  /* printInvalidities(mesh, adj_array, edge, adj.getSize()); */
  switch (mode) {
    case NIJK:
      objF = new InternalEdgeReshapeObjFunc(adapt, edge, tet, computeFValNIJKL);
      break;
    case DETJ:
      objF = new InternalEdgeReshapeObjFunc(adapt, edge, tet, computeFValDetJ);
      break;
    default:
      break;
  }
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;

  l = new LBFGS(tol, iter, x0, objF);

  apf::MeshEntity* ed[6];
  int thisTETnum = 0;

  for (std::size_t i = 0; i < adj.getSize(); i++) {
    if (adj[i] == tet) thisTETnum = 1;
    else thisTETnum = 0;
    mesh->getDownward(adj[i], 1, ed);
    int edgeIndex = apf::findIn(ed, 6, edge);
    printf("reshape tried on %d edge, TET %d; ", edgeIndex, thisTETnum);
    crv_dbg::printTetNumber(mesh, adj[i], "debug_num_tet");
  }

  bool hasDecreased = false;
  invaliditySize = 0;

  if (l->run() && thisTetSize > 1) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);


    for (int i = 0; i < adj.getSize(); i++) {
      std::vector<int> aiNew = crv::getAllInvalidities(mesh, adj[i]);
      invaliditySize = invaliditySize + aiNew.size();
      hasDecreased = hasDecreased || (aiNew.size() > sizeHolder[i]);
    }

    if (hasDecreased == false ) {
      //invaliditySize = 0;
      //makeMultipleEntityMesh(mesh, adj_array, edge, "after_cavity_of_edge_", adj.getSize());
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "after_cavity_indv_tet_of_edge_", adj.getSize());
      /* printInvalidities(mesh, adj_array, edge, adj.getSize()); */
      /* std::cout<<"--------------------------------------"<<std::endl; */
      return true;
    }
    else {
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "after_cavity_indv_tet_of_edge_", adj.getSize());
      /* objF->restoreInitialNodes(); */
      /* printInvalidities(mesh, adj_array, edge, adj.getSize()); */
      std::cout<<"Size DID NOT decrease"<<std::endl;
      std::cout<<"--------------------------------------"<<std::endl;
      return false;
    }
  }
  else {
    //finalX = l->currentX; 
    //objF->setNodes(finalX); 
    //makeMultipleEntityMesh(mesh, adj_array, edge, "after_cavity_of_edge_", adj.getSize());
    if (thisTetSize == 0) {
      std::cout<<" No Optimization tried"<<std::endl;
      std::cout<<"--------------------------------------"<<std::endl;
      return false;
    }
    std::cout<<"*****Edge Optim FAILURE" <<std::endl;
    std::cout<<"--------------------------------------"<<std::endl;
    return false;
  }
}

void CrvBoundaryEdgeOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvBoundaryEdgeOptim :: setTol(double tolerance)
{
  tol = tolerance;
}

bool CrvBoundaryEdgeOptim :: run(int &invaliditySize)
{
  std::vector<int> sizeHolder;
  apf::MeshEntity* adj_array[99];
  apf::Adjacent adj;
  mesh->getAdjacent(edge, 3, adj);
  int cInvT = 0;
  int thisTetSize = 0;

  for (int i = 0; i < adj.getSize(); i++) {
    adj_array[i] = adj[i];
    std::vector<int> ai = crv::getAllInvalidities(mesh, adj[i]);
    if (adj[i] == tet) thisTetSize = ai.size();
    sizeHolder.push_back(ai.size());
  }

  //std::vector<int> ai = crv::getAllInvalidities(mesh,tet);
  //makeMultipleEntityMesh(mesh, adj_array, edge, "before_cavity_of_edge_", adj.getSize());
  //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "before_cavity_indv_tet_of_edge_", adj.getSize());
  /* printTetNumber(mesh, tet); */
  /* printInvalidities(mesh, adj_array, edge, adj.getSize()); */
  switch (mode) {
    case NIJK:
      objF = new BoundaryEdgeReshapeObjFunc(adapt, edge, tet, computeFValNIJKL);
      break;
    case DETJ:
      objF = new BoundaryEdgeReshapeObjFunc(adapt, edge, tet, computeFValDetJ);
      break;
    default:
      break;
  }
  std::vector<double> x0 = objF->getInitialGuess();



  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  l = new LBFGS(tol, iter, x0, objF);
  apf::MeshEntity* ed[6];
  int thisTETnum = 0;

  for (std::size_t i = 0; i < adj.getSize(); i++) {
    if (adj[i] == tet) thisTETnum = 1;
    else thisTETnum = 0;
    mesh->getDownward(adj[i], 1, ed); 
    int edgeIndex = apf::findIn(ed, 6, edge);
    printf("reshape tried on %d Medge, TET %d ", edgeIndex, thisTETnum);
    /* printTetNumber(mesh, adj[i]); */
  }

  bool hasDecreased = false;
  invaliditySize = 0;

  if (l->run() && thisTetSize > 0) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);

    for (int i = 0; i < adj.getSize(); i++) {
      std::vector<int> aiNew = crv::getAllInvalidities(mesh, adj[i]);
      invaliditySize = invaliditySize + aiNew.size();
      hasDecreased = hasDecreased || (aiNew.size() > sizeHolder[i]);
    }

    if (hasDecreased == false) {
      //makeMultipleEntityMesh(mesh, adj_array, edge, "after_cavity_of_edge_", adj.getSize());
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, edge, "after_cavity_indv_tet_of_edge_", adj.getSize());
      printInvalidities(mesh, adj_array, edge, adj.getSize());
      std::cout<<"--------------------------------------"<<std::endl;
      return true;
    }
    else {
      objF->restoreInitialNodes();
      printInvalidities(mesh, adj_array, edge, adj.getSize());
      std::cout<<"size did not decrease"<<std::endl;
      std::cout<<"--------------------------------------"<<std::endl;
      return false;
    }
  }
  else {
    //finalX = l->currentX;
    //objF->setNodes(finalX);
    if (thisTetSize == 0) {
      std::cout<<"No Optimization tried"<<std::endl;
      std::cout<<"--------------------------------------"<<std::endl;
      return false;
    }
    std::cout<<"*****Optim FAILURE"<<std::endl;
    std::cout<<"--------------------------------------"<<std::endl;
    return false;
  }
}

void CrvFaceOptim :: setMaxIter(int n)
{
  iter = n;
}

void CrvFaceOptim :: setTol(double tolerance)
{
  tol = tolerance;
}

bool CrvFaceOptim :: run(int &invaliditySize)
{
  std::vector<int> sizeHolder;
  apf::MeshEntity* adj_array[99];
  apf::Adjacent adj;
  mesh->getAdjacent(face, 3, adj);
  int thisTetSize = 0;
  for (int i = 0; i < adj.getSize(); i++) {
    adj_array[i] = adj[i];
    //std::vector<int> ai = crv::getAllInvalidities(mesh, adj[i]);
    //if (adj[i] == tet) thisTetSize = ai.size();
    //sizeHolder.push_back(ai.size());   
  }

  std::vector<int> ai = crv::getAllInvalidities(mesh, tet);
  invaliditySize = ai.size();
  //makeMultipleEntityMesh(mesh, adj_array, face, "before_cavity_of_face_", adj.getSize());
  //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, face, "before_cavity_indv_tet_of_face_", adj.getSize());
  /* printTetNumber(mesh, tet); */
  printInvalidities(mesh, adj_array, face, adj.getSize());
  switch (mode) {
    case NIJK:
      objF = new FaceReshapeObjFunc(adapt, face, tet, computeFValNIJKL);
      break;
    case DETJ:
      objF = new FaceReshapeObjFunc(adapt, face, tet, computeFValDetJ);
      break;
    default:
      break;
  }
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  l = new LBFGS(tol, iter, x0, objF);


  apf::MeshEntity* fc[4];
  int thisTETnum = 0;

  for (std::size_t i = 0; i < adj.getSize(); i++) {
    if (adj[i] == tet) thisTETnum = 1;
    else thisTETnum = 0;
    mesh->getDownward(adj[i], 2, fc);
    int faceIndex = apf::findIn(fc, 4, face);
    printf("reshape tried on %d face, TET %d ",faceIndex, thisTETnum);
    /* printTetNumber(mesh, adj[i]); */
  }

  //bool hasDecreased = false;
  //invaliditySize = 0;

  //if (l->run() && thisTetSize > 0) {
  if (l->run() && invaliditySize > 0) {
    finalX = l->currentX;
    fval = l->fValAfter;
    objF->setNodes(finalX);

    std::vector<int> aiNew = crv::getAllInvalidities(mesh, tet);
    /*
    for (int i = 0; i < adj.getSize(); i++) {
      std::vector<int> aiNew = crv::getAllInvalidities(mesh, adj[i]);
      invaliditySize = invaliditySize + aiNew.size();
      hasDecreased = hasDecreased || (aiNew.size() > sizeHolder[i]);
    }
*/
   // if (hasDecreased == false ) {
    if (aiNew.size() <= invaliditySize) {
      invaliditySize = aiNew.size();
      //makeMultipleEntityMesh(mesh, adj_array, face, "after_cavity_of_face_", adj.getSize());
      //makeIndividualTetsFromFacesOrEdges(mesh, adj_array, face, "after_cavity_indv_tet_of_face_", adj.getSize());
      printInvalidities(mesh, adj_array, face, adj.getSize());
      std::cout<<"----------------------------------------------------"<<std::endl;

      return true;
    }
    else {
      objF->restoreInitialNodes();
      printInvalidities(mesh, adj_array, face, adj.getSize());
      std::cout<<"Size didNOT decrease"<<std::endl;
      std::cout<<"----------------------------------------------------"<<std::endl;
      return false;
    }
  }
  else {
    if (thisTetSize == 0) {
      std::cout<<" No Optimization tried"<<std::endl;
      std::cout<<"-------------------------------------------"<<std::endl;
      return false;
    }
    std::cout<<"*****FaceOptim FAILURE"<<std::endl;
    return false;
  }
}


}
