/* #include "crvEdgeOptim.h" */
/* #include "LBFGS.h" */
#include "crvOptimizations.h"
#include "crv.h"
#include "gmi.h"
#include "apfMDS.h"
#include "apfNumbering.h"
#include "crvQuality.h"
#include "crvBezier.h"
#include "crvMath.h"
#include "crvDBG.h"
#include <iostream>
#include "apfMatrix.h"
#include <Eigen/Core>
#include <LBFGS.h>

/* static int global_counter = 0; */
/* static apf::MeshEntity* tetra[100]; */
/* static int number = 0; */


using Eigen::VectorXd;
using namespace LBFGSpp;


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

class LBFGSFunctor
{
  private:
    int n;
    crv::ObjFunction* obj;
  public:
    LBFGSFunctor(int _n, crv::ObjFunction* _obj) : n(_n), obj(_obj) {}
    double operator () (const VectorXd &x, VectorXd &grad)
    {
      std::vector<double> x_vec;
      x_vec.clear();
      for (int i = 0; i < n; i++)
        x_vec.push_back(x[i]);
      std::vector<double> g_vec;
      g_vec = obj->getGrad(x_vec);
      for (int i = 0; i < n; i++)
        grad[i] = g_vec[i];

      return obj->getValue(x_vec);
    }
};

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
  InternalEdgeReshapeObjFunc *objF = new InternalEdgeReshapeObjFunc(mesh, edge, tet);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  LBFGSFunctor fun(x0.size(), objF);
  LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 100;
  LBFGSSolver<double> solver(param);
  VectorXd x_init = VectorXd::Zero(x0.size());
  for (int i = 0; i < x0.size(); i++) {
    x_init[i] = x0[i];
    finalX.push_back(0);
  }

  double fx;
  int niter = solver.minimize(fun, x_init, fx);

  printf("niter and fx are %d and %f\n", niter, fx);

  /* LBFGS *l = new LBFGS(tol, iter, x0, objF); */

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

  if (niter <= param.max_iterations && thisTetSize > 0) {
    for (int i = 0; i < x0.size(); i++)
      finalX[i] = x_init[i];
    fval = fx;
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

  BoundaryEdgeReshapeObjFunc *objF = new BoundaryEdgeReshapeObjFunc(mesh, edge, tet);
  std::vector<double> x0 = objF->getInitialGuess();

  LBFGSFunctor fun(x0.size(), objF);
  LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 100;
  LBFGSSolver<double> solver(param);
  VectorXd x_init = VectorXd::Zero(x0.size());
  for (int i = 0; i < x0.size(); i++) {
    x_init[i] = x0[i];
    finalX.push_back(0);
  }

  double fx;
  int niter = solver.minimize(fun, x_init, fx);

  printf("niter and fx are %d and %f\n", niter, fx);



  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  /* LBFGS *l = new LBFGS(tol, iter, x0, objF); */
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

  if (niter < param.max_iterations && thisTetSize > 0) {
    for (int i = 0; i < x0.size(); i++)
      finalX[i] = x_init[i];
    fval = fx;
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
  FaceReshapeObjFunc *objF = new FaceReshapeObjFunc(mesh, face, tet);
  std::vector<double> x0 = objF->getInitialGuess();
  //double f0 = objF->getValue(x0);
  //std::cout<< "fval at x0 " << f0<<std::endl;
  /* LBFGS *l = new LBFGS(tol, iter, x0, objF); */

  LBFGSFunctor fun(x0.size(), objF);
  LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 100;
  LBFGSSolver<double> solver(param);
  VectorXd x_init = VectorXd::Zero(x0.size());
  for (int i = 0; i < x0.size(); i++) {
    x_init[i] = x0[i];
    finalX.push_back(0);
  }

  double fx;
  int niter = solver.minimize(fun, x_init, fx);

  printf("niter and fx are %d and %f\n", niter, fx);



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
  if (niter < param.max_iterations && invaliditySize > 0) {
    for (int i = 0; i < x0.size(); i++)
      finalX[i] = x_init[i];
    fval = fx;
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
