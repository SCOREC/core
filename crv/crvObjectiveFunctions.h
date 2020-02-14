#ifndef CRVOBJECTIVEFUNCTIONS_H
#define CRVOBJECTIVEFUNCTIONS_H

#include <iostream>
#include <vector>

#include <apf.h>
#include <apfMesh.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfShape.h>
#include <apfNumbering.h>

using namespace std;

namespace crv
{

class ObjFunction
{
  public:
    ObjFunction(){};
    virtual int getSpaceDim() = 0;
    virtual double getValue(const vector<double> &x) = 0;
    // TODO :: can we do this once for all the objective functions?
    virtual vector<double> getGrad(const vector<double> &_x) = 0;
};


// TODO: can these be reused
/* void getInitEdgeN(); */
/* void getInitFaceN(); */
/* void getInitTetN(); */

// Internal Edge Objective Functions
class InternalEdgeReshapeObjFunc : public ObjFunction
{
  public:
    InternalEdgeReshapeObjFunc(apf::Mesh2* m, apf::MeshEntity* e, apf::MeshEntity* t) :
    mesh(m), edge(e), tet(t)
    {
      P = mesh->getShape()->getOrder();
      d = mesh->getDimension();
      getSpaceDim();
      getInitEdgeN();
      getInitFaceN();
      getInitTetN();
    }
    ~InternalEdgeReshapeObjFunc(){}
    int getSpaceDim();
    double getValue(const vector<double> &x);
    vector<double> getGrad(const vector<double> &x);
    vector<double> getInitialGuess();
    void setNodes(const vector<double> &x);
    void restoreInitialNodes();
    int P; //order
    int d; //dimension
  private:
    void getInitEdgeN();
    void getInitFaceN();
    void getInitTetN();
    vector<double> convertNodeVectorToX(
    	const vector<apf::Vector3> &en,
    	const vector<apf::Vector3> &fn,
    	const vector<apf::Vector3> &tn);
    vector<apf::Vector3> convertXtoNodeVector(const vector<double> &x);
    void blendTris(
    	const vector<apf::Vector3> &edgeNodes,
    	vector<apf::Vector3> &faceNodes);
    /* void blendTets( */
    /* 	const vector<apf::Vector3> &edgeNodes, */
    /* 	const vector<apf::Vector3> faceNodes, */
    /* 	vector<apf::Vector3> &tetNodes); */
    void updateNodes(
    	const vector<apf::Vector3> &ed,
    	const vector<apf::Vector3> &fa,
    	const vector<apf::Vector3> &te);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    vector<apf::Vector3> ien;
    vector<apf::Vector3> ifn;
    vector<apf::Vector3> itn;
};


// Boundary Edge Objective Function
class BoundaryEdgeReshapeObjFunc : public ObjFunction
{
  public:
    BoundaryEdgeReshapeObjFunc(
    	apf::Mesh2* m,
    	apf::MeshEntity* e,
    	apf::MeshEntity* t) :
    mesh(m), edge(e), tet(t)
    {
      P = mesh->getShape()->getOrder();
      d = mesh->getDimension();
      getSpaceDim();
      getInitEdgeN();
      getInitFaceN();
      getInitTetN();
    }
    ~BoundaryEdgeReshapeObjFunc(){}
    int getSpaceDim();
    double getValue(const vector<double> &x);
    vector<double> getGrad(const vector<double> &x);
    vector<double> getInitialGuess();
    void setNodes(const vector<double> &x);
    void restoreInitialNodes();
    int P; //order
    int d; //dimension
  private:
    void getInitEdgeN();
    void getInitFaceN();
    void getInitTetN();
    vector<double> getParamCoords();
    vector<apf::Vector3> convertParamCoordsToNodeVector(
    	const vector<double> &x);
    vector<apf::Vector3> convertXtoNodeVector(const vector<double> &x);
    void blendTris(
    	const vector<apf::Vector3> &edgeNodes,
    	vector<apf::Vector3> &faceNodes);
    /* void blendTets( */
    /* 	const vector<apf::Vector3> &edgeNodes, */
    /* 	const vector<apf::Vector3> faceNodes, */
    /* 	vector<apf::Vector3> &tetNodes); */
    vector<apf::Vector3> getFaceControlPointsFromInterpolatingPoints(
    	apf::MeshEntity* face,
    	const vector<apf::Vector3> &faceInterpolatingP);
    void updateNodes(
    	const vector<apf::Vector3> &ed,
    	const vector<apf::Vector3> &fa,
    	const vector<apf::Vector3> &te,
    	bool isInitialX);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    vector<apf::Vector3> ien;
    vector<apf::Vector3> ifn;
    vector<apf::Vector3> itpfn;
    vector<apf::Vector3> itn;
};

class FaceReshapeObjFunc : public ObjFunction
{
  public:
    FaceReshapeObjFunc(apf::Mesh2* m, apf::MeshEntity* f, apf::MeshEntity* t) :
    mesh(m), face(f), tet(t)
    {
      P = mesh->getShape()->getOrder();
      d = mesh->getDimension();
      getSpaceDim();
      getInitFaceN();
      getInitTetN();
    }
    ~FaceReshapeObjFunc() {}
    int getSpaceDim();
    double getValue(const vector<double> &x);
    vector<double> getGrad(const vector<double> &_x);
    vector<double> getInitialGuess();
    void setNodes(const vector<double> &x);
    void restoreInitialNodes();
    int P; //order
    int d; //dimension
  private:
    void getInitFaceN();
    void getInitTetN();
    vector<double> convertNodeVectorToX(const vector<apf::Vector3> &fn);
    vector<apf::Vector3> convertXtoNodeVector(const vector<double> &x);
    //void blendTets(const std::vector<apf::Vector3> &edgeNodes, const std::vector<apf::Vector3> faceNodes, std::vector<apf::Vector3> &tetNodes);
    void updateNodes(
    	const vector<apf::Vector3> &fa,
    	const vector<apf::Vector3> &te);
  protected:
    apf::Mesh2* mesh;
    apf::MeshEntity* face;
    apf::MeshEntity* tet;
    vector<apf::Vector3> ifn;
    vector<apf::Vector3> itn;
};

}

#endif
