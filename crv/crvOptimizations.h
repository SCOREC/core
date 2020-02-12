#ifndef CRVOPTIMIZATIONS_H
#define CRVOPTIMIZATIONS_H

#include <iostream>
#include <pcu_util.h>
#include <vector>
#include "apf.h"
#include "apfMesh.h"
#include "apfMesh2.h"
#include "apfShape.h"

#include "crvObjectiveFunctions.h"


namespace crv {

class CrvEntityOptim
{
  public:
    CrvEntityOptim() {}
    virtual void setMaxIter(int n) = 0;
    virtual void setTol(double tolerance) = 0;
    virtual bool run(int &invaliditySize) = 0;
};

class CrvInternalEdgeOptim : public CrvEntityOptim
{
  public:
    CrvInternalEdgeOptim(apf::Mesh2* m, apf::MeshEntity* e, apf::MeshEntity* t) :
    	mesh(m), edge(e), tet(t) {}
    ~CrvInternalEdgeOptim(){}

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};

class CrvBoundaryEdgeOptim : public CrvEntityOptim
{
  public:
    CrvBoundaryEdgeOptim(apf::Mesh2* m, apf::MeshEntity* e, apf::MeshEntity* t) :
      mesh(m), edge(e), tet(t) {}
    ~CrvBoundaryEdgeOptim(){}

  public:
    void setMaxIter(int n);
    void setTol(double tolerance);
    bool run(int &invaliditySize);
  public:
    apf::Mesh2* mesh;
    apf::MeshEntity* edge;
    apf::MeshEntity* tet;
    int iter;
    double tol;
    std::vector<double> finalX;
    double fval;
};

class CrvFaceOptim : public CrvEntityOptim
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
