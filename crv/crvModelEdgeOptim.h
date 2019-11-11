#ifndef CRVMODELEDGEOPTIM_H
#define CRVMODELEDGEOPTIM_H

#include <iostream>
#include <pcu_util.h>
#include <vector>
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "LBFGS.h"
#include "apfShape.h"

namespace crv{

class CrvModelEdgeReshapeObjFunc : public ObjFunction
{
  public:
    CrvModelEdgeReshapeObjFunc(apf::Mesh2* m, apf::MeshEntity* e, apf::MeshEntity* t) : 
  	mesh(m), edge(e), tet(t)
  {
    P = mesh->getShape()->getOrder();
    d = mesh->getDimension();
    getSpaceDim();
    getVolume();
    getInitEdgeN();
    getInitFaceN();
    getInitTetN();
  }
  ~CrvModelEdgeReshapeObjFunc(){}

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
    std::vector<double> getParamCoords();
    std::vector<apf::Vector3> convertParamCoordsToNodeVector(const std::vector<double> &x);
    std::vector<apf::Vector3> convertXtoNodeVector(const std::vector<double> &x);
    void blendTris(const std::vector<apf::Vector3> &edgeNodes, std::vector<apf::Vector3> &faceNodes);
    //void blendTets(const std::vector<apf::Vector3> &edgeNodes, const std::vector<apf::Vector3> faceNodes, std::vector<apf::Vector3> &tetNodes);
    std::vector<apf::Vector3> getFaceControlPointsFromInterpolatingPoints(apf::MeshEntity* face, const std::vector<apf::Vector3> &faceInterpolatingP);
    void updateNodes(std::vector<apf::Vector3> ed, std::vector<apf::Vector3> fa, std::vector<apf::Vector3> te, bool isInitialX);
    std::vector<double> getVolume();
    double computeFValOfElement(apf::NewArray<apf::Vector3> &nodes, double volm);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    std::vector<double> vol;
    std::vector<apf::Vector3> ien;
    std::vector<apf::Vector3> ifn;
    std::vector<apf::Vector3> itpfn;
    std::vector<apf::Vector3> itn;

};  

class CrvModelEdgeOptim
{
  public:
    CrvModelEdgeOptim(apf::Mesh2* m, apf::MeshEntity* e, apf::MeshEntity* t) :
    	mesh(m), edge(e), tet(t) {}
    ~CrvModelEdgeOptim(){}

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run();
  public:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};
}

#endif
