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

static double getLinearVolPhys(apf::Mesh2* m, apf::MeshEntity* e)
{
  apf::MeshEntity* vs[12];
  int n = m->getDownward(e, 0, vs);
  apf::Vector3 coords[12];
  for (int i = 0; i < n; i++) {
    m->getPoint(vs[i], 0, coords[i]);
  }

  if (m->getType(e) == apf::Mesh::TRIANGLE)
  {
    return apf::cross(coords[1]-coords[0], coords[2]-coords[0]).getLength()/2.;
  }
  else if (m->getType(e) == apf::Mesh::TET)
  {
    apf::Matrix3x3 J;
    J[0] = coords[1] - coords[0];
    J[1] = coords[2] - coords[0];
    J[2] = coords[3] - coords[0];
    return apf::getDeterminant(J) / 6.;
  }
  else
    PCU_ALWAYS_ASSERT_VERBOSE(0,
      "Not implemented for entities of type other than tri or tet!");
  return 0.;
}

namespace crv
{

class ObjFunction
{
  public:
    ObjFunction(){};
    virtual int getSpaceDim() = 0;
    virtual double getTol() = 0;
    virtual double getValue(const vector<double> &x) = 0;
    // TODO :: can we do this once for all the objective functions?
    vector<double> getGrad(const vector<double> &_x)
    {
      double h;
      vector<double> x = _x;
      double eps = this->getTol();
      vector<double> g;
      for (size_t i = 0; i < x.size(); i++) {
	h = abs(x[i]) > eps ? eps * abs(x[i]) : eps;

	// forward diff
	x[i] += h;
	double ff = this->getValue(x);
	x[i] -= h;

	// backward diff
	x[i] -= h;
	double fb = this->getValue(x);
	x[i] += h;

	g.push_back( (ff - fb) / 2./ h );
      }
      return(g);
    }
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
      tol = cbrt(getLinearVolPhys(m, t)) * 1.e-3;
    }
    ~InternalEdgeReshapeObjFunc(){}
    int getSpaceDim();
    double getTol() {return tol;}
    double getValue(const vector<double> &x);
    /* vector<double> getGrad(const vector<double> &x); */
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
    double tol;
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
      tol = cbrt(getLinearVolPhys(m, t)) * 1.e-3;
    }
    ~BoundaryEdgeReshapeObjFunc(){}
    int getSpaceDim();
    double getTol() {return tol;}
    double getValue(const vector<double> &x);
    /* vector<double> getGrad(const vector<double> &x); */
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
    double tol;
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
      tol = cbrt(getLinearVolPhys(m, t)) * 1.e-3;
    }
    ~FaceReshapeObjFunc() {}
    int getSpaceDim();
    double getTol() {return tol;}
    double getValue(const vector<double> &x);
    /* vector<double> getGrad(const vector<double> &_x); */
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
    double tol;
};

}

#endif
