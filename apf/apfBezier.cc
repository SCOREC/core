/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <cassert>
#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>

#include "apf.h"
#include "apfMesh.h"
#include "apfNew.h"
#include "apfShape.h"

namespace apf {

// negative -> flipped relative to canonical
// relies on e0 being always ordered correctly
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

template <int P>
class BezierShape : public FieldShape
{
public:
  class Vertex : public EntityShape
  {
  public:
    void getValues(Mesh*, MeshEntity*,
        Vector3 const&, NewArray<double>& values) const
    {
    	values.allocate(1);
    	values[0] = 1.0;
    }
    void getLocalGradients(Mesh*, MeshEntity*,
        Vector3 const&, NewArray<Vector3>&) const
    {
    }
    int countNodes() const {return 1;}
  };
  class Edge : public EntityShape
  {
  public:
    void getValues(Mesh* /*m*/, MeshEntity* /*e*/,
        Vector3 const& xi, NewArray<double>& values) const
    {
      double t = 0.5*(xi[0]+1.);
      values.allocate(P+1);
      for(int i = 1; i < P; ++i)
        values[i+1] = binomial(P,i)
        * pow(1.0-t,P-i)*pow(t, i);
      values[0] = pow(1-t, P);
      values[1] = pow(t, P);

    }
    void getLocalGradients(Mesh* /*m*/, MeshEntity* /*e*/,
        Vector3 const& xi, NewArray<Vector3>& grads) const
    {
      double t = 0.5*(xi[0]+1.);
      grads.allocate(P+1);
      for(int i = 1; i < P; ++i)
        grads[i+1] = Vector3(binomial(P,i) * (i-P*t)
            * pow(1.0-t,P-1-i)*pow(t, i-1)/2.,0,0);
      grads[0] = Vector3(-P*pow(1-t, P-1)/2.,0,0);
      grads[1] = Vector3(P*pow(t, P-1)/2.,0,0);
    }
    int countNodes() const {return P+1;}
    void alignSharedNodes(Mesh*,
        MeshEntity*, MeshEntity*, int order[])
    {
      (void)order;
    }
  };
  class Triangle : public EntityShape
  {
  public:
    Triangle()
    {
      int m1[] = {1,2,0};
      int m2[] = {1,4,2,3,5,0};
      int m3[] = {1,5,6,2,4,9,7,3,8,0};
      int m4[] = {1,6,7,8,2,5,13,14,9,4,12,10,3,11,0};
      int m5[] = {1,7,8,9,10,2,6,17,18,19,11,5,16,20,12,4,15,13,3,14,0};
      int m6[] = {1,8,9,10,11,12,2,7,21,22,23,24,13,6,20,27,25,14,5,19,
          26,15,4,18,16,3,17,0};
      int* maps[6] = {m1,m2,m3,m4,m5,m6};
      for(int i = 0; i < (P+1)*(P+2)/2; ++i)
        map[i] = maps[P-1][i];
    }
    void getValues(Mesh* m, MeshEntity* e, Vector3 const& xi,
        NewArray<double>& values) const
    {
      values.allocate((P+1)*(P+2)/2);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){

        for(int i = 0; i < P+1; ++i)
          for(int j = 0; j < P+1-i; ++j)
            values[map[j*(P+1)+i-j*(j-1)/2]] =
                factorial(P)/factorial(i)/factorial(j)/factorial(P-i-j)
                *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],P-i-j);

      } else {
        // Triangular Blending

        for(int i = 0; i < 3; ++i)
          values[i] = -xii[i];
        // zero the rest, the face node weight is always zero
        for(int i = 3; i < (P+1)*(P+2)/2; ++i)
          values[i] = 0.0;

        double x;
        Vector3 xv;
        NewArray<double> v;

        int const (*tev)[2] = tri_edge_verts;

        for(int i = 0; i < 3; ++i){
          x = xii[tev[i][0]]+xii[tev[i][1]];
          if(x < 1e-13) continue;
          xv[0] = 2.0*(xii[tev[i][1]]/x)-1.0;

          getBezier(3,P)->getEntityShape(Mesh::EDGE)
                ->getValues(0,0,xv,v);

          for(int j = 0; j < 2; ++j)
            values[tev[i][j]]   += v[j]*x;
          for(int j = 0; j < (P-1); ++j)
            values[3+i*(P-1)+j] += v[2+j]*x;
        }
      }
    }
    void getLocalGradients(Mesh* m, MeshEntity* e, Vector3 const& xi,
        NewArray<Vector3>& grads) const
    {
      grads.allocate((P+1)*(P+2)/2);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      Vector3 gxii[3] = {Vector3(-1,-1,0),Vector3(1,0,0),Vector3(0,1,0)};

      ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        for(int i = 0; i < P+1; ++i)
          for(int j = 0; j < P+1-i; ++j)
            grads[map[j*(P+1)+i-j*(j-1)/2]] =
              gxii[0]*factorial(P)/factorial(i)/factorial(j)/factorial(P-i-j)
              *pow(xii[0],i-1)*pow(xii[1],j)*pow(xii[2],P-i-j-1)
              *(i*(1.-xii[1])-(P-j)*xii[0]) +
              gxii[1]*factorial(P)/factorial(i)/factorial(j)/factorial(P-i-j)
              *pow(xii[0],i)*pow(xii[1],j-1)*pow(xii[2],P-i-j-1)
              *(j*(1.-xii[0])-(P-i)*xii[1]);

      } else {
        for(int i = 0; i < 3; ++i)
          grads[i] = gxii[i]*-1;

        for(int i = 3; i < (P+1)*(P+2)/2; ++i)
          grads[i] = Vector3(0,0,0);

        double x;
        Vector3 xv, gx;

        NewArray<double> v;
        NewArray<Vector3> gv;
        MeshEntity * edges[3];

        int const (*tev)[2] = tri_edge_verts;

        m->getDownward(e,1,edges);
        for(int i = 0; i < 3; ++i){
          x = xii[tev[i][0]]+xii[tev[i][1]];
          if(x < 1e-13) continue;
          gx = gxii[tev[i][0]]+gxii[tev[i][1]];

          xv[0] = 2.0*(xii[tev[i][1]]/x)-1.0;

          getBezier(3,P)->getEntityShape(Mesh::EDGE)
                    ->getValues(m,edges[i],xv,v);

          getBezier(3,P)->getEntityShape(Mesh::EDGE)
                    ->getLocalGradients(m,edges[i],xv,gv);

          for(int j = 0; j < 2; ++j)
            grads[tev[i][j]]   += gx*v[j]
              + (gxii[tev[i][1]]-gx*xii[tev[i][1]]/x)*gv[j][0]*2.;
          for(int j = 0; j < (P-1); ++j)
            grads[3+i*(P-1)+j] += gx*v[j+2]
              + (gxii[tev[i][1]]-gx*xii[tev[i][1]]/x)*gv[j+2][0]*2.;
        }
      }
    }
    int countNodes() const {return (P+1)*(P+2)/2;}
    void alignSharedNodes(Mesh* m,
        MeshEntity* elem, MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(flip){
        for(int i = 0; i < P-1; ++i)
          order[i] = P-2-i;
      } else {
        for(int i = 0; i < P-1; ++i)
          order[i] = i;
      }
    }
  private:
    int map[(P+1)*(P+2)/2];
  };
  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh* m, MeshEntity* e,
        Vector3 const& xi, NewArray<double>& values) const
    {
    	values.allocate(2*P*P+2);

      double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};

      for(int i = 0; i < 4; ++i)
        values[i] = xii[i];
      for(int i = 4; i < 2*P*P+2; ++i)
        values[i] = 0.0;

      double x;
      Vector3 xv;
      NewArray<double> v;

      MeshEntity* faces[4];

      int nE = P-1;
      int nF = (P-1)*(P-2)/2;

      int const (*tev)[2] = tet_edge_verts;
      int const (*ttv)[3] = tet_tri_verts;

      for(int i = 0; i < 6; ++i){
        x = xii[tev[i][0]]+xii[tev[i][1]];
        if(x < 1e-13) continue;

        xv[0] = 2.0*(xii[tev[i][1]]/x)-1.0;

        getBezier(3,P)->getEntityShape(Mesh::EDGE)
        	->getValues(0,0,xv,v);

        for(int j = 0; j < 2; ++j) // vertices
          values[tev[i][j]] += -v[j]*x;
        for(int j = 0; j < nE; ++j)// edge nodes
          values[4+i*nE+j]  += -v[2+j]*x;

      }

      m->getDownward(e,2,faces);
      for(int i = 0; i < 4; ++i){

        x = 0.;
        for(int j = 0; j < 3; ++j)
          x += xii[ttv[i][j]];
        if(x < 1e-13) continue;

        int shift = (m->getModelType(m->toModel(faces[i]))
            == m->getDimension());
        for(int j = 0; j < 3; ++j)
          xv[j] = xii[ttv[i][(j+shift) % 3]];

        xv = xv/x;

        getBezier(3,P)->getEntityShape(Mesh::TRIANGLE)
        	->getValues(m,faces[i],xv,v);

        for(int j = 0; j < 3; ++j) // vertices
          values[ttv[i][j]] += v[j]*x;

        // Edge contributions from faces
        // Edges are the first 3*nE entries per face
        // also not necessarily aligned, canonically
        for(int k = 0; k < 3; ++k){
        	int l = tet_tri_edges[i][k];
        	if(flip_tet_tri_edges[i][k] == false){
        		for(int j = 0; j < nE; ++j)
        			values[4+l*nE+j] += v[3+k*nE+j]*x;
        	} else { // we need to flip
        		for(int j = 0; j < nE; ++j)
        			values[4+l*nE+j] += v[3+k*nE+nE-1-j]*x;
        	}
        }
        for(int j = 0; j < nF; ++j) // face nodes
          values[4+6*nE+i*nF+j] +=  v[3+3*nE+j]*x;
      } // done faces
    }
    void getLocalGradients(Mesh* m, MeshEntity* e,
        Vector3 const& xi,
        NewArray<Vector3>& grads) const
    {
      grads.allocate(2*P*P+2);

      double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
      Vector3 gxii[4] = {Vector3(-1,-1,-1),Vector3(1,0,0),
        Vector3(0,1,0),Vector3(0,0,1)};

      for(int i = 0; i < 4; ++i)
        grads[i] = gxii[i];
      for(int i = 4; i < 2*P*P+2; ++i)
        grads[i] = Vector3(0,0,0);

      double x;
      Vector3 xv, gx;
      Vector3 gxv[2];
      NewArray<double> v;
      NewArray<Vector3> gv;

      MeshEntity* faces[4];

      int nE = P-1;
      int nF = (P-1)*(P-2)/2;

      int const (*tev)[2] = tet_edge_verts;
      int const (*ttv)[3] = tet_tri_verts;

      for(int i = 0; i < 6; ++i){
        x = xii[tev[i][0]]+xii[tev[i][1]];
        if(x < 1e-13) continue;

        gx = gxii[tev[i][0]]+gxii[tev[i][1]];

        xv[0] = 2.0*(xii[tev[i][1]]/x)-1.0;

        getBezier(3,P)->getEntityShape(Mesh::EDGE)
          ->getValues(0,0,xv,v);

        getBezier(3,P)->getEntityShape(Mesh::EDGE)
          ->getLocalGradients(0,0,xv,gv);

        for(int j = 0; j < 2; ++j) // vertices
          grads[tev[i][j]] += gx*(-v[j])
            - (gxii[tev[i][1]]-gx*xii[tev[i][1]]/x)*gv[j][0]*2.;

        for(int j = 0; j < nE; ++j)// edge nodes
          grads[4+i*nE+j]  += gx*(-v[j+2])
            - (gxii[tev[i][1]]-gx*xii[tev[i][1]]/x)*gv[j+2][0]*2.;
      }

      m->getDownward(e,2,faces);
      for(int i = 0; i < 4; ++i){
        x = 0.;
        for(int j = 0; j < 3; ++j)
          x += xii[ttv[i][j]];

        if(x < 1e-13) continue;

        int shift =
            (m->getModelType(m->toModel(faces[i])) == m->getDimension());
        for(int j = 0; j < 3; ++j)
          xv[j] = xii[ttv[i][(j+shift) % 3]];

        gx = gxii[ttv[i][0]] + gxii[ttv[i][1]] + gxii[ttv[i][2]];

        xv = xv/x;
        // actually x*gxv
        gxv[0] = gxii[ttv[i][(0+shift) % 3]]-gx*xv[0];
        gxv[1] = gxii[ttv[i][(1+shift) % 3]]-gx*xv[1];

        getBezier(3,P)->getEntityShape(Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);

        getBezier(3,P)->getEntityShape(Mesh::TRIANGLE)
          ->getLocalGradients(m,faces[i],xv,gv);

        for(int j = 0; j < 3; ++j) // vertices
          grads[ttv[i][j]] += gx*v[j]
            + gxv[0]*gv[j][0] + gxv[1]*gv[j][1];

        for(int k = 0; k < 3; ++k){
          int l = tet_tri_edges[i][k];
          if(flip_tet_tri_edges[i][k] == false){
            for(int j = 0; j < nE; ++j)
              grads[4+l*nE+j] += gx*v[3+k*nE+j]
                + gxv[0]*gv[3+k*nE+j][0]
                + gxv[1]*gv[3+k*nE+j][1];

          } else { // we need to flip
            for(int j = 0; j < nE; ++j)
              grads[4+l*nE+j] += gx*v[3+k*nE+nE-1-j]
                + gxv[0]*gv[3+k*nE+nE-1-j][0]
                + gxv[1]*gv[3+k*nE+nE-1-j][1];
          }
        }
        for(int j = 0; j < nF; ++j) // face nodes
          grads[4+6*nE+i*nF+j] += gx*v[3+3*nE+j]
            + gxv[0]*gv[3+3*nE+j][0]
            + gxv[1]*gv[3+3*nE+j][1];
      } // done faces
    }
    int countNodes() const {return 2*P*P+2;}
    void alignSharedNodes(Mesh* m,
        MeshEntity* elem, MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
    	getAlignment(m,elem,shared,which,flip,rotate);
      if(m->getType(shared) == Mesh::EDGE){
        if(!flip)
          for(int i = 0; i < P-1; ++i)
            order[i] = i;
        else
          for(int i = 0; i < P-1 ; ++i)
            order[i] = P-2-i;
        return;
      }
      // must be a triangle
      int n = (P-1)*(P-2)/2;
      int l = n/3; //loops
      if(!flip)
        for(int i = 0; i < n; ++i)
          order[i] = (i+l*(3-rotate)) % (3*l);
      else {
        int shift[4] = {0,0,1,4};
        for(int i = 0; i < n; ++i)
          order[i] = (n-1-i+(n-shift[l])-l*rotate) % (3*l);
      }
      if(n % l) order[3*l] = 3*l;
    }
  };
  EntityShape* getEntityShape(int type)
  {
    static Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    static Tetrahedron tet;
    static EntityShape* shapes[Mesh::TYPES] =
    {&vertex,   //vertex
     &edge,     //edge
     &triangle, //triangle
     NULL,      //quad
     &tet,      //tet
     NULL,      //hex
     NULL,      //prism
     NULL};     //pyramid
    return shapes[type];
  }
  bool hasNodesIn(int dimension)
  {
    if ((dimension == 0)||
        ((dimension == 1) && P > 1)||
        ((dimension == 2) && P > 2))
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    static int nodes[Mesh::TYPES] =
    {1,                 //vertex
     P-1,               //edge
     (P-1)*(P-2)/2,     //triangle
     0,                 //quad
     0,                 //tet
     0,                 //hex
     0,                 //prism
     0};                //pyramid
    return nodes[type];
  }
  int getOrder() {return P;}
};

template <int P>
class BezierCurve : public BezierShape<P>
{
public:
  const char* getName() const {return name.c_str();}
  BezierCurve<P>() {
    std::stringstream ss;
    ss << "BezierCurve_" << P;
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  void getNodeXi(int type, int node, Vector3& xi)
  {
    static double eP2[1] = {0.0};
    static double eP3[2] = {-0.4306648,0.4306648};
    static double eP4[3] = {-0.6363260,0.0,0.6363260};
    static double eP5[4] = {-0.7485748,-0.2765187,0.2765187,0.7485748};
    static double eP6[5] = {-0.8161268,-0.4568660,0.0,
        0.4568660,0.8161268};

    static double* edgePoints[6] =
    {eP2, eP2, eP3, eP4, eP5, eP6 };

    if(type == Mesh::EDGE && P > 1){
      xi = Vector3(edgePoints[P-1][node],0,0);
    } else {
      xi = Vector3(0,0,0);
    }
  }
protected:
  std::string name;
};
template <int P>
class BezierSurface : public BezierShape<P>
{
public:
  const char* getName() const {return name.c_str();}
  BezierSurface() {
    std::stringstream ss;
    ss << "BezierSurface_" << P;
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  void getNodeXi(int type, int node, Vector3& xi)
  {
    static double eP2[1] = {0.0};
    static double eP3[2] = {-0.4503914,0.4503914};
    static double eP4[3] = {-0.6612048,0.0,0.6612048};
    static double eP5[4] = {-0.7732854,-0.2863522,0.2863522,0.7732854};
    static double eP6[5] = {-0.8388042,-0.469821,0.0,
      0.469821,0.8388042};
    static double* edgePoints[6] =
    {eP2, eP2, eP3, eP4, eP5, eP6 };
    if (type == Mesh::EDGE) {
      xi = Vector3(edgePoints[P-1][node],0,0);
    } else if (type == Mesh::TRIANGLE) {
      xi = Vector3(1./3.,1./3.,1./3.);
      if(node == (P-1)*(P-2)/2-1 && P % 3 == 0){
        return;
      } else { // technically only two of these numbers are needed
        switch (P) {
          case 1:
          case 2:
          case 3:
            fail("expected P >= 4");
          case 4:
            xi[(node+2) % 3] = 0.5582239;
            xi[(node+0) % 3] = 0.220880;
            xi[(node+1) % 3] = 0.220880;
            break;
          case 5:
            if(node % 2 == 0) {
              xi[(node/2+2) % 3] = 0.6949657;
              xi[(node/2+0) % 3] = 0.1525171;
              xi[(node/2+1) % 3] = 0.1525171;
            } else {
              xi[((node-1)/2+2) % 3] = 0.4168658;
              xi[((node-1)/2+0) % 3] = 0.4168658;
              xi[((node-1)/2+1) % 3] = 0.1662683;
            }
            break;
          case 6:
            if (node % 3 == 0) {
              xi[(node/3+2) % 3] = 0.7805723;
              xi[(node/3+0) % 3] = 0.1097139;
              xi[(node/3+1) % 3] = 0.1097139;
            } else if ((node-1) % 3 == 0) {
              xi[((node-1)/3+2) % 3] = 0.5586077;
              xi[((node-1)/3+0) % 3] = 0.3157892;
              xi[((node-1)/3+1) % 3] = 0.1256031;
            } else if ((node-2) % 3 == 0) {
              xi[((node-2)/3+2) % 3] = 0.3157892;
              xi[((node-2)/3+0) % 3] = 0.5586077;
              xi[((node-2)/3+1) % 3] = 0.1256031;
            }
            break;
        }
      }
    } else {
      xi = Vector3(0,0,0);
    }
  }
protected:
  std::string name;
};

template <int P>
class GregorySurface : public FieldShape
{
public:
  const char* getName() const {return name.c_str();}
  GregorySurface() {
    std::stringstream ss;
    ss << "GregorySurface_" << P;
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  class Triangle : public EntityShape
  {
  public:
    Triangle()
    {
      int m1[] = {1,2,0};
      int m2[] = {1,4,2,3,5,0};
      int m3[] = {1,5,6,2,4,9,7,3,8,0};
      int m4[] = {1,6,7,8,2,5,13,14,9,4,12,10,3,11,0};
      int m5[] = {1,7,8,9,10,2,6,17,18,19,11,5,16,20,12,4,15,13,3,14,0};
      int m6[] = {1,8,9,10,11,12,2,7,21,22,23,24,13,6,20,27,25,14,5,19,
          26,15,4,18,16,3,17,0};
      int* maps[6] = {m1,m2,m3,m4,m5,m6};
      for(int i = 0; i < (P+1)*(P+2)/2; ++i)
        map[i] = maps[P-1][i];
    }
    void getValues(Mesh* /*m*/, MeshEntity* /*e*/, Vector3 const& /*xi*/,
        NewArray<double>& values) const
    {
      values.allocate((P+1)*(P+2)/2);

    }
    void getLocalGradients(Mesh* /*m*/, MeshEntity* /*e*/, Vector3 const& /*xi*/,
        NewArray<Vector3>& grads) const
    {
      grads.allocate((P+1)*(P+2)/2);

    }
    int countNodes() const {return (P+1)*(P+2)/2;}
    void alignSharedNodes(Mesh* m,
        MeshEntity* elem, MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(flip){
        for(int i = 0; i < P-1; ++i)
          order[i] = P-2-i;
      } else {
        for(int i = 0; i < P-1; ++i)
          order[i] = i;
      }
    }
  private:
    int map[(P+1)*(P+2)/2];
  };
  class Tetrahedron : public EntityShape
  {
  public:

    void getValues(Mesh* /*m*/, MeshEntity* /*e*/,
        Vector3 const& /*xi*/, NewArray<double>& values) const
    {
      values.allocate(2*P*P+2);
    }
    void getLocalGradients(Mesh* /*m*/, MeshEntity* /*e*/,
        Vector3 const& /*xi*/,
        NewArray<Vector3>& grads) const
    {
      grads.allocate(2*P*P+2);

    }
    int countNodes() const {return 2*P*P+2;}
    void alignSharedNodes(Mesh* m,
        MeshEntity* elem, MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(m->getType(shared) == Mesh::EDGE){
        if(!flip)
          for(int i = 0; i < P-1; ++i)
            order[i] = i;
        else
          for(int i = 0; i < P-1 ; ++i)
            order[i] = P-2-i;
        return;
      }
      // must be a triangle
      int n = (P-1)*(P-2)/2;
      int l = n/3; //loops
      if(!flip)
        for(int i = 0; i < n; ++i)
          order[i] = (i+l*(3-rotate)) % (3*l);
      else {
        int shift[4] = {0,0,1,4};
        for(int i = 0; i < n; ++i)
          order[i] = (n-1-i+(n-shift[l])-l*rotate) % (3*l);
      }
      if(n % l) order[3*l] = 3*l;
    }
  };
  EntityShape* getEntityShape(int type)
  {
    static BezierShape<1>::Vertex vertex;
    static BezierShape<4>::Edge edge4;
    static BezierShape<5>::Edge edge5;
    static BezierShape<6>::Edge edge6;

    static Triangle triangle;
    static Tetrahedron tet;
    static EntityShape* shapes[Mesh::TYPES] =
    {&vertex,   //vertex
     NULL,      //edge
     &triangle, //triangle
     NULL,      //quad
     &tet,      //tet
     NULL,      //hex
     NULL,      //prism
     NULL};     //pyramid
    static EntityShape* edges[3] = {&edge4,&edge5,&edge6};

    if(type == Mesh::EDGE)
      return edges[P-4];
    return shapes[type];
  }
  bool hasNodesIn(int dimension)
  {
    if ((dimension == 0)||
        ((dimension == 1) && P > 1)||
        ((dimension == 2) && P > 2))
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    static int nodes[Mesh::TYPES] =
    {1,                 //vertex
     P-1,               //edge
     (P-1)*(P-2)/2,     //triangle
     0,                 //quad
     0,                 //tet
     0,                 //hex
     0,                 //prism
     0};                //pyramid
    return nodes[type];
  }
  int getOrder() {return P;}
protected:
  std::string name;
};

static void getBezierCurveInterPtsToCtrlPts(int order,
    NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -0.970273514,0.333333333,2.71895067,-1.08201049,
      0.333333333,-0.970273514,-1.08201049,2.71895067};
  double e4[15] = {
      -1.43042029,-0.25,3.3954584,-1.46967987,0.754641763,
      0.953613524,0.953613524,-2.76673344,4.62623983,-2.76673344,
      -0.25,-1.43042029,0.754641763,-1.46967987,3.3954584};
  double e5[24] = {
      -1.88592269,0.2,4.05614416,-1.81653638,1.02954238,-0.58322747,
      1.85476912,-0.942961345,-5.01939998,6.96205914,-4.562341,2.70787405,
      -0.942961345,1.85476912,2.70787405,-4.562341,6.96205914,-5.01939998,
      0.2,-1.88592269,-0.58322747,1.02954238,-1.81653638,4.05614416};
  double e6[35] = {
      -2.338908,-0.166666667,4.70907763,-2.14695479,
      1.2670886,-0.800405899,0.476769119,
      3.03457283,0.935563201,-7.82909978,9.74813268,
      -6.60581336,4.37362215,-2.65697772,
      -2.27592963,-2.27592963,6.3088041,-9.70710792,
      12.3484669,-9.70710792,6.3088041,
      0.935563201,3.03457283,-2.65697772,4.37362215,
      -6.60581336,9.74813268,-7.82909978,
      -0.166666667,-2.338908,0.476769119,-0.800405899,
      1.2670886,-2.14695479,4.70907763};
  double* table[5] = {
      e2,e3,e4,e5,e6};
  int nb = order-1;
  int ni = order+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[order-2][i*ni+j];

}
static void getBezierShapeInterPtsToCtrlPts(int order, int type,
    NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -1.00596379,0.333333333,2.69317846,-1.020548,
      0.333333333,-1.00596379,-1.020548,2.69317846};
  double e4[15] = {
      -1.52680421,-0.25,3.37567603,-1.28732567,0.688453848,
      1.01786947,1.01786947,-2.70941992,4.3831009,-2.70941992,
      -0.25,-1.52680421,0.688453848,-1.28732567,3.37567603};
  double e5[24] = {
      -2.06136018,0.2,4.05936189,-1.52513018,0.846118039,-0.518989558,
      2.07392858,-1.03068009,-5.04403528,6.48508084,-4.12567866,2.64138461,
      -1.03068009,2.07392858,2.64138461,-4.12567866,6.48508084,-5.04403528,
      0.2,-2.06136018,-0.518989558,0.846118039,-1.52513018,4.05936189};
  double e6[35] = {
      -2.60465922,-0.166666667,4.7431726,-1.74847212,
      0.99151406,-0.630691219,0.415802564,
      3.50945493,1.04186369,-8.00777337,8.99109434,
      -5.80152935,3.85156704,-2.58467729,
      -2.6320912,-2.6320912,6.39664545,-8.91824704,
      11.3073856,-8.91824704,6.39664545,
      1.04186369,3.50945493,-2.58467729,3.85156704,
      -5.80152935,8.99109434,-8.00777337,
      -0.166666667,-2.60465922,0.415802564,-0.630691219,
      0.99151406,-1.74847212,4.7431726};
  double f3[10] = {
      0.505963791,0.505963791,0.505963791,-0.836315229,-0.836315229,
      -0.836315229,-0.836315229,-0.836315229,-0.836315229,4.5};
  double f4[45] = {
      1.47384491,-0.487291885,-0.487379988,-2.15710583,-0.789454397,
      1.00218533,0.456974235,0.41607983,0.457065362,1.00239351,
      -0.78966413,-2.15716739,7.06638632,-2.00331581,-2.00355008,
      -0.487291885,1.47384491,-0.487379988,1.00218533,-0.789454397,
      -2.15710583,-2.15716739,-0.78966413,1.00239351,0.457065362,
      0.41607983,0.456974235,-2.00331581,7.06638632,-2.00355008,
      -0.487302875,-0.487302875,1.47397924,0.456969833,0.415988284,
      0.456969833,1.00222185,-0.789441831,-2.15739071,-2.15739071,
      -0.789441831,1.00222185,-2.00331581,-2.00331581,7.06655155};
  double f5[126] = {
      2.95550925,0.485063144,0.485064069,-4.01515326,-0.633917331,
      1.14097233,-1.04127809,-0.319296612,-0.217804382,-0.217804388,
      -0.319297318,-1.04128021,1.14097477,-0.633919719,-4.01515431,
      10.1689603,-3.08692404,1.17750746,0.897196394,1.17750932,
      -3.08692737,-1.87982915,-1.87982915,0.569548438,3.55892311,-2.43777981,
      -2.43777981,3.55892311,1.77258998,0.975113153,0.114906986,
      -0.833742192,-0.833742192,0.114906986,0.975113153,1.77258998,
      -6.13761798,12.4785213,-6.13761798,-2.28036489,2.24753181,
      -2.28036489,0.485063144,2.95550925,0.485064069,-1.04127809,1.14097233,
      -0.633917331,-4.01515326,-4.01515431,-0.633919719,1.14097477,
      -1.04128021,-0.319297318,-0.217804388,-0.217804382,-0.319296612,
      1.17750746,-3.08692404,10.1689603,-3.08692737,1.17750932,
      0.897196394,0.569547683,-1.87982942,-1.87983114,-0.833740807,0.114906984,
      0.975110986,1.77258965,3.55892419,-2.43777941,-2.43778357,
      3.55892769,1.77259212,0.975112397,0.114906311,-0.833740942,
      2.24752905,-2.28036294,-6.13761791,12.4785246,-6.13762225,
      -2.28036323,0.485063281,0.485063281,2.95551097,-0.319296574,-0.217803972,
      -0.217803972,-0.319296574,-1.04127851,1.14097258,-0.633917018,
      -4.015157,-4.015157,-0.633917018,1.14097258,-1.04127851,
      1.17750756,0.897195548,1.17750756,-3.08692488,10.1689625,
      -3.08692488,-1.87982942,0.569547683,-1.87983114,1.77258965,0.975110986,
      0.114906984,-0.833740807,-0.833740942,0.114906311,0.975112397,
      1.77259212,3.55892769,-2.43778357,-2.43777941,3.55892419,
      -6.13761791,-2.28036294,2.24752905,-2.28036323,-6.13762225,
      12.4785246};
  double f6[280] = {
      4.99079517,-0.479884384,-0.47988312,-6.45824664,-0.363965285,
      1.22336735,-1.23728258,1.04863932,0.239902581,0.147640596,
      0.0937167895,0.147640457,0.239901872,1.04863643,-1.23727948,
      1.22336374,-0.363960892,-6.45824438,13.8928375,-4.27628545,
      1.80891892,-0.815239911,-0.490056818,-0.490055805,-0.815237976,
      1.80891479,-4.27627655,1.32762381,
      -4.68945914,2.30341,-0.668682176,8.40126116,-5.32763856,
      -2.38073429,4.59608707,-4.67688886,-1.6138599,-0.831300245,
      -0.353270382,0.0169225284,0.78380706,1.11543265,-0.639627184,
      -0.406338791,1.63350587,4.42017843,-13.101009,20.0696901,
      -10.7012814,5.23778052,2.40699204,0.724701049,-2.11625937,
      3.07107055,-2.09328984,-4.18119991,
      2.30341,-4.68945914,-0.668682176,-4.67688886,4.59608707,
      -2.38073429,-5.32763856,8.40126116,4.42017843,1.63350587,
      -0.406338791,-0.639627184,1.11543265,0.78380706,0.0169225284,
      -0.353270382,-0.831300245,-1.6138599,5.23778052,-10.7012814,
      20.0696901,-13.101009,-2.09328984,3.07107055,-2.11625937,
      0.724701049,2.40699204,-4.18119991,
      -0.479884384,4.99079517,-0.47988312,1.04863932,-1.23728258,
      1.22336735,-0.363965285,-6.45824664,-6.45824438,-0.363960892,
      1.22336374,-1.23727948,1.04863643,0.239901872,0.147640457,
      0.0937167895,0.147640596,0.239902581,-0.815239911,1.80891892,
      -4.27628545,13.8928375,-4.27627655,1.80891479,-0.815237976,
      -0.490055805,-0.490056818,1.32762381,
      -0.668681688,-4.68945906,2.30340803,1.11543116,-0.639623988,
      -0.406342307,1.6335101,4.42018047,8.40125889,-5.32764256,
      -2.38073106,4.59608452,-4.67688443,-1.61385789,-0.83130176,
      -0.353268755,0.0169208258,0.783807682,-2.11625804,3.07106668,
      -2.09328154,-13.1010089,20.0696815,-10.7012774,5.23777847,
      2.40698919,0.724701462,-4.18119971,
      -0.668681902,2.30340988,-4.6894572,0.783808233,0.0169213392,
      -0.353269103,-0.831301836,-1.61386073,-4.6768877,4.59608832,
      -2.38073512,-5.32763833,8.4012556,4.42017598,1.63350856,
      -0.406341222,-0.639625077,1.11543118,-2.11625834,0.724702974,
      2.40698863,5.23778035,-10.7012776,20.0696879,-13.1010089,
      -2.09328137,3.07106572,-4.1812002,
      -0.479884213,-0.479884213,4.99079329,0.239902612,0.147640874,
      0.09371655,0.147640874,0.239902612,1.04863883,-1.2372824,
      1.22336699,-0.363964411,-6.45824151,-6.45824151,-0.363964411,
      1.22336699,-1.2372824,1.04863883,-0.815239673,-0.490056323,
      -0.490056323,-0.815239673,1.80891797,-4.27628462,13.8928377,
      -4.27628462,1.80891797,1.3276242,
      2.30340988,-0.668681902,-4.6894572,-1.61386073,-0.831301836,
      -0.353269103,0.0169213392,0.783808233,1.11543118,-0.639625077,
      -0.406341222,1.63350856,4.42017598,8.4012556,-5.32763833,
      -2.38073512,4.59608832,-4.6768877,5.23778035,2.40698863,
      0.724702974,-2.11625834,3.07106572,-2.09328137,-13.1010089,
      20.0696879,-10.7012776,-4.1812002,
      -4.68945906,-0.668681688,2.30340803,4.42018047,1.6335101,
      -0.406342307,-0.639623988,1.11543116,0.783807682,0.0169208258,
      -0.353268755,-0.83130176,-1.61385789,-4.67688443,4.59608452,
      -2.38073106,-5.32764256,8.40125889,-13.1010089,-2.09328154,
      3.07106668,-2.11625804,0.724701462,2.40698919,5.23777847,
      -10.7012774,20.0696815,-4.18119971,
      2.7404102,2.7404102,2.74041055,-3.89671938,0.852566248,
      2.6291994,0.852566248,-3.89671938,-3.89671807,0.852569479,
      2.62919661,0.852568921,-3.89671824,-3.89671824,0.852568921,
      2.62919661,0.852569479,-3.89671807,9.21852962,-7.99944935,
      -7.99944935,9.21852962,-7.99944374,-7.99945144,9.21853253,
      -7.99945144,-7.99944374,23.4971758};
  double* table[10] = {
      e2,e3,e4,e5,e6,NULL,f3,f4,f5,f6};
  int nb = (type == Mesh::TRIANGLE) ? (order-1)*(order-2)/2 : order-1;
  int ni = (type == Mesh::TRIANGLE) ? (order+1)*(order+2)/2 : order+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[5*(type-1)+(order-2)][i*ni+j];
}

void getTransformationCoefficients(int order, int dim, int type,
    NewArray<double>& c){
  if(dim == 2)
    getBezierCurveInterPtsToCtrlPts(order,c);
  else
    getBezierShapeInterPtsToCtrlPts(order,type,c);
}

static FieldShape* getBezierCurve(int order)
{
  static BezierCurve<1> bezierCurve1;
  static BezierCurve<2> bezierCurve2;
  static BezierCurve<3> bezierCurve3;
  static BezierCurve<4> bezierCurve4;
  static BezierCurve<5> bezierCurve5;
  static BezierCurve<6> bezierCurve6;

  FieldShape* curveTable[6] = {&bezierCurve1,&bezierCurve2,
      &bezierCurve3,&bezierCurve4,&bezierCurve5,&bezierCurve6};
  return curveTable[order-1];
}

static FieldShape* getBezierSurface(int order)
{
  static BezierSurface<1> bezierSurface1;
  static BezierSurface<2> bezierSurface2;
  static BezierSurface<3> bezierSurface3;
  static BezierSurface<4> bezierSurface4;
  static BezierSurface<5> bezierSurface5;
  static BezierSurface<6> bezierSurface6;

  FieldShape* surfaceTable[6] = {&bezierSurface1,&bezierSurface2,
      &bezierSurface3,&bezierSurface4,&bezierSurface5,&bezierSurface6};
  return surfaceTable[order-1];
}

FieldShape* getBezier(int dimension, int order)
{
  return (dimension == 2) ? getBezierCurve(order) : getBezierSurface(order);
}

FieldShape* getGregory(int order)
{
  static GregorySurface<4> gregorySurface4;
  static GregorySurface<5> gregorySurface5;
  static GregorySurface<6> gregorySurface6;

  FieldShape* surfaceTable[3] = {&gregorySurface4, &gregorySurface5,
      &gregorySurface6};
  return surfaceTable[order-4];
}
}//namespace apf
