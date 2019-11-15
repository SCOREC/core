#ifndef CRVFACEOPTIM_H
#define CRVFACEOPTIM_H

#include <iostream>
#include <pcu_util.h>
#include <vector>
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "LBFGS.h"
#include "apfShape.h"

namespace crv{

class CrvFaceReshapeObjFunc : public ObjFunction
{
  public:
    CrvFaceReshapeObjFunc(apf::Mesh2* m, apf::MeshEntity* f, apf::MeshEntity* t) : 
  	mesh(m), face(f), tet(t)
  {
    P = mesh->getShape()->getOrder();
    d = mesh->getDimension();
    getSpaceDim();
    getVolume();
    getInitFaceN();
    getInitTetN();
  }
  ~CrvFaceReshapeObjFunc() {}

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
    void getInitFaceN();
    void getInitTetN();
    std::vector<double> convertNodeVectorToX(std::vector<apf::Vector3> fn);
    std::vector<apf::Vector3> convertXtoNodeVector(const std::vector<double> &x);
    //void blendTets(const std::vector<apf::Vector3> &edgeNodes, const std::vector<apf::Vector3> faceNodes, std::vector<apf::Vector3> &tetNodes);
    void updateNodes(std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te);
    std::vector<double> getVolume();
    double computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* face;
    apf::MeshEntity* tet;
    std::vector<double> vol;
    std::vector<apf::Vector3> ifn;
    std::vector<apf::Vector3> itn;
};  

class CrvFaceOptim
{
  public:
    CrvFaceOptim(apf::Mesh2* m, apf::MeshEntity* f, apf::MeshEntity* t) :
    	mesh(m), face(f), tet(t) {}
    ~CrvFaceOptim(){}

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    apf::Mesh2* mesh;
    apf::MeshEntity* face;
    apf::MeshEntity* tet;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};
}

#endif
