#ifndef CRVEDGEOPTIM_H
#define CRVEDGEOPTIM_H

#include <iostream>
#include <pcu_util.h>
#include <vector>
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "LBFGS.h"
#include "apfShape.h"

namespace crv{

class CrvEdgeReshapeObjFunc : public ObjFunction
{
  public:
    CrvEdgeReshapeObjFunc(apf::Mesh2* m, apf::MeshEntity* e) : 
  	mesh(m), edge(e)
  {
    P = mesh->getShape()->getOrder();
    d = mesh->getDimension();
    getSpaceDim();
    getVolume();
    getInitEdgeN();
    getInitFaceN();
    getInitTetN();
  }
  public:
    int getSpaceDim();
    double getValue(std::vector<double> &x);
    std::vector<double> getGrad(std::vector<double> &x);
    std::vector<double> getInitialGuess();
    void setNodes(std::vector<double> &x);
    void restoreInitialNodes();
    int P; //order
    int d; //dimension
    
  private:
    void getInitEdgeN();
    void getInitFaceN();
    void getInitTetN();
    std::vector<double> convertNodeVectorToX(std::vector<apf::Vector3> en, std::vector<apf::Vector3> fn, std::vector<apf::Vector3> tn);
    std::vector<apf::Vector3> convertXtoNodeVector(const std::vector<double> &x);
    void blendTris(const std::vector<apf::Vector3> &edgeNodes, std::vector<apf::Vector3> &faceNodes);
    //void blendTets(const std::vector<apf::Vector3> &edgeNodes, const std::vector<apf::Vector3> faceNodes, std::vector<apf::Vector3> &tetNodes);
    void updateNodes(std::vector<apf::Vector3> ed, std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te);
    std::vector<double> getVolume();
    double computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    std::vector<double> vol;
    std::vector<apf::Vector3> ien;
    std::vector<apf::Vector3> ifn;
    std::vector<apf::Vector3> itn;
};  

class CrvEdgeOptim
{
  public:
    CrvEdgeOptim(apf::Mesh2* m, apf::MeshEntity* e) :
    	mesh(m), edge(e) {}
  public:
    void setMaxIter(int n);
    void setTol(double t);
    bool run();
  public:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};
}

#endif
