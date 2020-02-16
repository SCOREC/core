#ifndef CRVOPTIMIZATIONS_H
#define CRVOPTIMIZATIONS_H

#include <iostream>
#include <pcu_util.h>
#include <vector>
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"

#include "LBFGS.h"
#include "crvObjectiveFunctions.h"


namespace crv {

class CrvEntityOptim
{
  public:
    CrvEntityOptim() {}
    virtual void setMaxIter(int n) = 0;
    virtual void setTol(double tolerance) = 0;
    virtual bool run(int &invaliditySize) = 0;
  protected:
    LBFGS* l;
};

class CrvInternalEdgeOptim : public CrvEntityOptim
{
  public:
    CrvInternalEdgeOptim(
    	crv::Adapt* a,
    	apf::MeshEntity* e,
    	apf::MeshEntity* t,
    	int m) :
    	adapt(a), edge(e), tet(t), mode(m)
    {
      mesh = adapt->mesh;
    }
    ~CrvInternalEdgeOptim()
    {
      delete objF;
      delete l;
    }

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    crv::Adapt* adapt;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    int mode;
    apf::Mesh2* mesh;
    InternalEdgeReshapeObjFunc* objF;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};

class CrvBoundaryEdgeOptim : public CrvEntityOptim
{
  public:
    CrvBoundaryEdgeOptim(
    	crv::Adapt* a,
    	apf::MeshEntity* e,
    	apf::MeshEntity* t,
    	int m) :
      adapt(a), edge(e), tet(t), mode(m)
    {
      mesh = adapt->mesh;
    }
    ~CrvBoundaryEdgeOptim()
    {
      delete objF;
      delete l;
    }

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    crv::Adapt* adapt;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    int mode;
    apf::Mesh2* mesh;
    BoundaryEdgeReshapeObjFunc* objF;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};

class CrvFaceOptim : public CrvEntityOptim
{
  public:
    CrvFaceOptim(
    	crv::Adapt* a,
    	apf::MeshEntity* f,
    	apf::MeshEntity* t,
    	int m) :
      adapt(a), face(f), tet(t), mode(m)
    {
      mesh = adapt->mesh;
    }
    ~CrvFaceOptim()
    {
      delete objF;
      delete l;
    }

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    crv::Adapt* adapt;
    apf::MeshEntity* face;
    apf::MeshEntity* tet;
    int mode;
    apf::Mesh2* mesh;
    FaceReshapeObjFunc* objF;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};

}

#endif
