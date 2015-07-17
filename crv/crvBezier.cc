/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <math.h>
#include <stdio.h>
#include "crvBezier.h"
/* see bezier.tex */

namespace crv {

static int P = 1;

static bool useBlend()
{
  return (getBlendingOrder() != 0);
}

/* Elevates a bezier curve from order n to order n+r */
void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r)
{

  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),edge);

  apf::Vector3 pt;
  apf::NewArray<apf::Vector3> p;
  apf::getVectorNodes(elem,p);

  assert(m->getType(edge) == apf::Mesh::EDGE);

  // reorder p into geometric ordering
  apf::NewArray<int> map(n+1);
  map[0] = 0; map[n] = 1;

  for(int i = 1; i < n; ++i)
    map[i] = i+1;
  for(int i = 1; i < n+r; ++i){
    pt.zero();
    for(int j = std::max(0,i-r); j <= std::min(i,n); ++j)
      pt += p[map[j]]*binomial(n,j)*
      binomial(r,i-j)/binomial(n+r,i);
    m->setPoint(edge,i-1,pt);
  }

  apf::destroyElement(elem);
}

class BezierShape : public apf::FieldShape
{
public:
  class Vertex : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<double>& values) const
    {
      values.allocate(1);
      values[0] = 1.0;
    }
    void getLocalGradients(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
    }
    int countNodes() const {return 1;}
    void alignSharedNodes(apf::Mesh*,
        apf::MeshEntity*, apf::MeshEntity*, int order[])
    {
      (void)order;
    }
  };
  class Edge : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      double t = 0.5*(xi[0]+1.);
      values.allocate(P+1);
      for(int i = 1; i < P; ++i)
        values[i+1] = binomial(P,i)
        * pow(1.0-t,P-i)*pow(t, i);
      values[0] = pow(1-t, P);
      values[1] = pow(t, P);

    }
    void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      double t = 0.5*(xi[0]+1.);
      grads.allocate(P+1);
      for(int i = 1; i < P; ++i)
        grads[i+1] = apf::Vector3(binomial(P,i) * (i-P*t)
            * pow(1.0-t,P-1-i)*pow(t, i-1)/2.,0,0);
      grads[0] = apf::Vector3(-P*pow(1-t, P-1)/2.,0,0);
      grads[1] = apf::Vector3(P*pow(t, P-1)/2.,0,0);
    }
    int countNodes() const {return P+1;}
    void alignSharedNodes(apf::Mesh*,
        apf::MeshEntity*, apf::MeshEntity*, int order[])
    {
      (void)order;
    }
  };
  class Triangle : public apf::EntityShape
  {
  public:
    Triangle()
    {
      int m1[] = {2,0,1};
      int m2[] = {2,5,0,4,3,1};
      int m3[] = {2,7,8,0,6,9,3,5,4,1};
      int m4[] = {2,9,10,11,0,8,14,12,3,7,13,4,6,5,1};
      int m5[] = {2,11,12,13,14,0,10,19,20,15,3,9,18,16,4,8,17,5,7,6,1};
      int m6[] = {2,13,14,15,16,17,0,12,24,25,26,18,3,11,23,27,19,4,10,
          22,20,5,9,21,6,8,7,1};

      int* maps[6] = {m1,m2,m3,m4,m5,m6};
      for(int j = 1; j <= 6; ++j){
        map[j-1].allocate(curved_face_total[BEZIER][j-1]);
        for(int i = 0; i < curved_face_total[BEZIER][j-1]; ++i)
          map[j-1][i] = maps[j-1][i];
      }
    }
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(curved_face_total[BEZIER][P-1]);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      if(!useBlend() || m->getModelType(m->toModel(e)) != m->getDimension()){

        for(int i = 0; i < P+1; ++i)
          for(int j = 0; j < P+1-i; ++j)
            values[map[P-1][j*(P+1)+i-j*(j-1)/2]] =
                binomial(P,i)*binomial(P-i,j)
                *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],P-i-j);

      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(curved_face_total[BEZIER][P-1]);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

      if(!useBlend() || m->getModelType(m->toModel(e)) != m->getDimension()){
        for(int i = 1; i < P+1; ++i)
          for(int j = 1; j < P-i; ++j)
            grads[map[P-1][j*(P+1)+i-j*(j-1)/2]] =
                gxii[0]*binomial(P,i)*binomial(P-i,j)
                *pow(xii[0],i-1)*pow(xii[1],j)*pow(xii[2],P-i-j-1)
                *(i*(1.-xii[1])-(P-j)*xii[0]) +
                gxii[1]*binomial(P,i)*binomial(P-i,j)
                *pow(xii[0],i)*pow(xii[1],j-1)*pow(xii[2],P-i-j-1)
                *(j*(1.-xii[0])-(P-i)*xii[1]);

        // i = 0
        for(int j = 1; j < P; ++j)
          grads[map[P-1][j*(P+1)-j*(j-1)/2]] =
            (gxii[0]*(j-P)*xii[1] +
            gxii[1]*(j*(1.-xii[0])-P*xii[1]))
            *binomial(P,j)
            *pow(xii[1],j-1)*pow(xii[2],P-j-1);

        // j = 0
        for(int i = 1; i < P; ++i)
          grads[map[P-1][i]] =
              (gxii[0]*(i*(1.-xii[1])-P*xii[0]) +
               gxii[1]*(i-P)*xii[0])
              *binomial(P,i)
              *pow(xii[0],i-1)*pow(xii[2],P-i-1);

        // k = 0
        for(int i = 1, j = P-1; i < P; ++i, --j)
          grads[map[P-1][j*(P+1)+i-j*(j-1)/2]] =
              (gxii[0]*i*xii[1] + gxii[1]*j*xii[0])
              *binomial(P,i)*binomial(P-i,j)
              *pow(xii[0],i-1)*pow(xii[1],j-1);

        // i = j = 0
        grads[map[P-1][0]] = gxii[2]*P*pow(xii[2],P-1);
        // i = k = 0
        grads[map[P-1][((P+1)*(P+2))/2-1]] = gxii[1]*P*pow(xii[1],P-1);
        // j = k = 0
        grads[map[P-1][P]] = gxii[0]*P*pow(xii[0],P-1);

      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

    }
    int countNodes() const {return curved_face_total[BEZIER][P-1];}
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
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
    apf::NewArray<int> map[6];
  };
  class Tetrahedron : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      if(!useBlend() && P <= 4){
        values.allocate(curved_tet_total[BEZIER][P-1]);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        for(int i = 0; i < 4; ++i)
          values[i] = pow(xii[i],P);

        int nE = P-1;
        int nF = curved_face_internal[BEZIER][P-1];
        int nT = curved_tet_internal[BEZIER][P-1];
        int const (*tev)[2] = apf::tet_edge_verts;
        int const (*ttv)[3] = apf::tet_tri_verts;

        for(int i = 0; i < 6; ++i){
          for(int j = 0; j < nE; ++j){// edge nodes
            int ijkl[4] = {0,0,0,0};
            ijkl[tev[i][0]] = P-j-1;
            ijkl[tev[i][1]] = j+1;
            values[4+i*nE+j] = binomial(P,ijkl[0])*binomial(P-ijkl[0],ijkl[1])
                *binomial(P-ijkl[0]-ijkl[1],ijkl[2])
                *pow(xii[0],ijkl[0])*pow(xii[1],ijkl[1])
                *pow(xii[2],ijkl[2])*pow(xii[3],ijkl[3]);
          }
        }

        for(int i = 0; i < 4; ++i){
          for(int j = 0; j < nF; ++j){ // face nodes
            int ijkl[4] = {0,0,0,0};
            ijkl[ttv[i][0]] = 1;
            ijkl[ttv[i][1]] = 1;
            ijkl[ttv[i][2]] = 1;
            if(P == 4)
              ijkl[ttv[i][j]] += 1;
            values[4+6*nE+i*nF+j] = binomial(P,ijkl[0])
                *binomial(P-ijkl[0],ijkl[1])
                *binomial(P-ijkl[0]-ijkl[1],ijkl[2])
                *pow(xii[0],ijkl[0])*pow(xii[1],ijkl[1])
                *pow(xii[2],ijkl[2])*pow(xii[3],ijkl[3]);
          }
        } // done faces

        for(int j = 0; j < nT; ++j){
          int ijkl[4] = {1,1,1,1};
          values[4+6*nE+4*nF+j] = binomial(P,ijkl[0])
              *binomial(P-ijkl[0],ijkl[1])
              *binomial(P-ijkl[0]-ijkl[1],ijkl[2])
              *pow(xii[0],ijkl[0])*pow(xii[1],ijkl[1])
              *pow(xii[2],ijkl[2])*pow(xii[3],ijkl[3]);
        }
      } else {
        values.allocate(blended_tet_total[BEZIER][P-1]);
        BlendedTetGetValues(m,e,xi,values);
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      if(!useBlend() && P <= 4){
        grads.allocate(curved_tet_total[BEZIER][P-1]);


        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
            apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

        for(int i = 0; i < 4; ++i)
          grads[i] = gxii[i]*P*pow(xii[i],P);

        for(int i = 4; i < curved_tet_total[BEZIER][P-1]; ++i)
          grads[i].zero();

      } else {
        grads.allocate(blended_tet_total[BEZIER][P-1]);
        BlendedTetGetLocalGradients(m,e,xi,grads);
      }
    }
    int countNodes() const {
      if(!useBlend())
        return curved_tet_total[BEZIER][P-1];
      else
        return blended_tet_total[BEZIER][P-1];
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(m->getType(shared) == apf::Mesh::EDGE){
        if(!flip)
          for(int i = 0; i < P-1; ++i)
            order[i] = i;
        else
          for(int i = 0; i < P-1 ; ++i)
            order[i] = P-2-i;
        return;
      }
      // must be a triangle
      int n = curved_face_internal[BEZIER][P-1];
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
  apf::EntityShape* getEntityShape(int type)
  {
    static Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    static Tetrahedron tet;
    static apf::EntityShape* shapes[apf::Mesh::TYPES] =
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
        ((dimension == 2) && P > 2)||
        ((dimension == 3) && P > 3
            && !useBlend()))
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    switch (type) {
      case apf::Mesh::VERTEX:
        return 1;
      case apf::Mesh::EDGE:
        return P-1;
      case apf::Mesh::TRIANGLE:
        return curved_face_internal[BEZIER][P-1];
      case apf::Mesh::TET:
        if(!useBlend()){
          return curved_tet_internal[BEZIER][P-1];
        } else
          return 0;
      default:
        return 0;
    }
  }
  int getOrder() {return std::max(P,getBlendingOrder());}
};

class BezierCurve : public BezierShape
{
public:
  const char* getName() const {return name.c_str();}
  BezierCurve() {
    std::stringstream ss;
    ss << "BezierCurve";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    static double eP2[1] = {0.0};
    static double eP3[2] = {-0.4306648,0.4306648};
    static double eP4[3] = {-0.6363260,0.0,0.6363260};
    static double eP5[4] = {-0.7485748,-0.2765187,0.2765187,0.7485748};
    static double eP6[5] = {-0.8161268,-0.4568660,0.0,
        0.4568660,0.8161268};

    static double* edgePoints[6] =
    {eP2, eP2, eP3, eP4, eP5, eP6 };

    if(type == apf::Mesh::EDGE && P > 1){
      xi[0] = edgePoints[P-1][node];
    } else {
      getBezierNodeXi(type,P,node,xi);
    }
  }
protected:
  std::string name;
};

class BezierSurface : public BezierShape
{
public:
  const char* getName() const {return name.c_str();}
  BezierSurface() {
    std::stringstream ss;
    ss << "BezierSurface";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    getBezierNodeXi(type,P,node,xi);
  }
protected:
  std::string name;
};

class GregorySurface4 : public apf::FieldShape
{
public:
  const char* getName() const {return name.c_str();}
  GregorySurface4() {
    std::stringstream ss;
    ss << "GregorySurface";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  class Triangle : public apf::EntityShape
  {
  public:
    Triangle()
    {
      int m[] = {2,9,10,11,0,8,14,12,3,7,13,4,6,5,1};

      for(int i = 0; i < 15; ++i)
        map[i] = m[i];

      index[0][0] = 2; index[1][0] = 0; index[2][0] = 1;
      index[0][1] = 1; index[1][1] = 2; index[2][1] = 0;

      pairs[0][0] = 0; pairs[1][0] = 1; pairs[2][0] = 2;
      pairs[0][1] = 5; pairs[1][1] = 3; pairs[2][1] = 4;

    }
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(curved_face_total[GREGORY][3]);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        for(int i = 0; i < 5; ++i)
          for(int j = 0; j < 5-i; ++j)
            values[map[j*5+i-j*(j-1)/2]] =
                binomial(4,i)*binomial(4-i,j)
                *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],4-i-j);

        for(int i = 0; i < 3; ++i){
          double xiix,x = xii[index[i][0]] + xii[index[i][1]];
          double bernstein = values[12+pairs[i][0]];
          values[12+pairs[i][1]] = 0.;

          if(x < 1e-12)
            xiix = 0.5;
          else
            xiix = xii[index[i][0]]/x;
          values[12+pairs[i][0]] = bernstein*xiix;
          values[12+pairs[i][1]] = bernstein*(1.-xiix);
        }
      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(curved_face_total[GREGORY][3]);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};
      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        for(int i = 1; i < 5; ++i)
          for(int j = 1; j < 5-i; ++j)
            grads[map[j*5+i-j*(j-1)/2]] =
                gxii[0]*binomial(4,i)*binomial(4-i,j)
                *pow(xii[0],i-1)*pow(xii[1],j)*pow(xii[2],4-i-j-1)
                *(i*(1.-xii[1])-(4-j)*xii[0]) +
                gxii[1]*binomial(4,i)*binomial(4-i,j)
                *pow(xii[0],i)*pow(xii[1],j-1)*pow(xii[2],4-i-j-1)
                *(j*(1.-xii[0])-(4-i)*xii[1]);

        // i = 0
        for(int j = 1; j < 4; ++j)
          grads[map[j*(4+1)-j*(j-1)/2]] =
            (gxii[0]*(j-4.)*xii[1] +
            gxii[1]*(j*(1.-xii[0])-4.*xii[1]))
            *binomial(4,j)
            *pow(xii[1],j-1)*pow(xii[2],4-j-1);

        // j = 0
        for(int i = 1; i < 4; ++i)
          grads[map[i]] =
              (gxii[0]*(i*(1.-xii[1])-4.*xii[0]) +
               gxii[1]*(i-4.)*xii[0])
              *binomial(4,i)
              *pow(xii[0],i-1)*pow(xii[2],4-i-1);

        // k = 0
        for(int i = 1, j = 3; i < 4; ++i, --j)
          grads[map[j*5+i-j*(j-1)/2]] =
              (gxii[0]*i*xii[1] + gxii[1]*j*xii[0])
              *binomial(4,i)*binomial(4-i,j)
              *pow(xii[0],i-1)*pow(xii[1],j-1);

        // i = j = 0
        grads[map[0]] = (gxii[0]+gxii[1])*(-4.)*pow(xii[2],3);
        // i = k = 0
        grads[map[((4+1)*(4+2))/2-1]] = gxii[1]*4.*pow(xii[1],3);
        // j = k = 0
        grads[map[4]] = gxii[0]*4.*pow(xii[0],3);

        double x, xiix;
        apf::Vector3 xv, gx;
        apf::NewArray<double> v;

        getValues(m,e,xi,v);

        for(int i = 0; i < 3; ++i){
          x  = xii[index[i][0]] + xii[index[i][1]];
          gx = gxii[index[i][0]] + gxii[index[i][1]];

          apf::Vector3 bernstein = grads[12+pairs[i][0]];
          grads[12+pairs[i][1]].zero();

          if(x < 1e-12)
              xiix = 0.5;
            else
              xiix = xii[index[i][0]]/x;

          grads[12+pairs[i][0]] = bernstein*xiix;
          grads[12+pairs[i][1]] = bernstein*(1.-xiix);

        }
      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

    }
    int countNodes() const
    {
      return curved_face_total[GREGORY][3];
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(flip){
        for(int i = 0; i < 3; ++i)
          order[i] = 2-i;
      } else {
        for(int i = 0; i < 3; ++i)
          order[i] = i;
      }
    }
  private:
    int map[15];
    int index[3][2];
    int pairs[3][2];
  };
  class Tetrahedron : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      values.allocate(blended_tet_total[GREGORY][3]);
      BlendedTetGetValues(m,e,xi,values);
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(blended_tet_total[GREGORY][3]);
      BlendedTetGetLocalGradients(m,e,xi,grads);
    }
    int countNodes() const {return blended_tet_total[GREGORY][3];}
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(m->getType(shared) == apf::Mesh::EDGE){
        if(!flip)
          for(int i = 0; i < 3; ++i)
            order[i] = i;
        else
          for(int i = 0; i < 3; ++i)
            order[i] = 2-i;
        return;
      }
      int orients[6][6] =
      {{0,1,2,3,4,5},{2,0,1,5,3,4},{1,2,0,4,5,3},
       {4,3,5,1,0,2},{3,5,4,0,2,1},{5,4,3,2,1,0}};
      for(int i = 0; i < 6; ++i)
        order[i] = orients[flip*3+rotate][i];
    }
  };
  apf::EntityShape* getEntityShape(int type)
  {
    static BezierShape::Vertex vertex;
    static BezierShape::Edge edge;
    static Triangle triangle;
    static Tetrahedron tet;
    static apf::EntityShape* shapes[apf::Mesh::TYPES] =
    {&vertex,   //vertex
     &edge,    //edge
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
    if (dimension == 3)
      return false;
    else
      return true;
  }
  int countNodesOn(int type)
  {
    static int nodes[apf::Mesh::TYPES] =
    {1,                 //vertex
     3,                 //edge
     6,                 //triangle
     0,                 //quad
     0,                 //tet
     0,                 //hex
     0,                 //prism
     0};                //pyramid
    return nodes[type];
  }
  /* These don't make sense for gregory patches */
  void getNodeXi(int type, int node, apf::Vector3& xi)
   {
    static double edgePoints[3] = {-0.6612048,0.0,0.6612048};
    if (type == apf::Mesh::EDGE) {
      xi[0] = edgePoints[node];
    } else if (type == apf::Mesh::TRIANGLE) {
      getBezierNodeXi(type,4,node,xi);
      if (node < 3){
        xi[(node+2) % 3] = 0.5582239;
        xi[(node+0) % 3] = 0.22088805;
        xi[(node+1) % 3] = 0.22088805;
      } else {
        xi[(node+1) % 3] = 0.5582239;
        xi[(node+2) % 3] = 0.22088805;
        xi[(node+0) % 3] = 0.22088805;
      }
    } else
      xi.zero();
   }
  int getOrder() {return std::max(4,getBlendingOrder());}
protected:
  std::string name;
};

class Nurbs : public apf::FieldShape
{
public:
  const char* getName() const {return name.c_str();}
  Nurbs() {
    std::stringstream ss;
    ss << "Nurbs";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  class Edge : public apf::EntityShape
  {
  public:
    Edge() {
      weights.allocate(P+1);
      for(int i = 0; i < P+1; ++i)
        weights[i] = 1.;
    }
    void getValues(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      getBezier(3,P)->getEntityShape(apf::Mesh::EDGE)->getValues(0,0,xi,values);

      double sum = 0.;
      for(int i = 0; i < P+1; ++i){
        values[i] *= weights[i];
        sum += values[i];
      }
      for(int i = 0; i < P+1; ++i){
        values[i] /= sum;
      }
    }
    void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      apf::NewArray<double> values;
      getBezier(3,P)->
          getEntityShape(apf::Mesh::EDGE)->getLocalGradients(0,0,xi,grads);
      getBezier(3,P)->getEntityShape(apf::Mesh::EDGE)->getValues(0,0,xi,values);

      double sum = 0.;
      for(int i = 0; i < P+1; ++i){
        values[i] *= weights[i];
        sum += values[i];
      }

      apf::Vector3 gsum = apf::Vector3(0,0,0);
      for(int i = 0; i < P+1; ++i){
        grads[i] = grads[i]*weights[i];
        gsum += grads[i];
      }
      for(int i = 0; i < P+1; ++i){
        grads[i] = (grads[i]*sum-gsum*values[i])/sum/sum;
      }
    }
    int countNodes() const {return P+1;}
    void alignSharedNodes(apf::Mesh*,
        apf::MeshEntity*, apf::MeshEntity*, int order[])
    {
      (void)order;
    }
    void setWeights(apf::NewArray<double>& w)
    {
      weights.allocate(P+1);
      for(int i = 0; i < P+1; ++i)
        weights[i] = w[i];
    }
  private:
    apf::NewArray<double> weights;
  };
  class Triangle : public apf::EntityShape
  {
  public:
    Triangle()
    {
      int m1[] = {2,0,1};
      int m2[] = {2,5,0,4,3,1};
      int m3[] = {2,7,8,0,6,9,3,5,4,1};
      int m4[] = {2,9,10,11,0,8,14,12,3,7,13,4,6,5,1};
      int m5[] = {2,11,12,13,14,0,10,19,20,15,3,9,18,16,4,8,17,5,7,6,1};
      int m6[] = {2,13,14,15,16,17,0,12,24,25,26,18,3,11,23,27,19,4,10,
          22,20,5,9,21,6,8,7,1};

      int* maps[6] = {m1,m2,m3,m4,m5,m6};
      for(int j = 1; j <= 6; ++j){
        map[j-1].allocate(curved_face_total[BEZIER][j-1]);
        for(int i = 0; i < curved_face_total[BEZIER][j-1]; ++i)
          map[j-1][i] = maps[j-1][i];
      weights.allocate(curved_face_total[BEZIER][P-1]);
      for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i)
        weights[i] = 1.;
      }
    }
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      getBezier(3,P)->
          getEntityShape(apf::Mesh::TRIANGLE)->getValues(m,e,xi,values);

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){

        double sum = 0.;
        for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i){
          values[i] *= weights[i];
          sum += values[i];
        }
        for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i){
          values[i]/=sum;
        }
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      getBezier(3,P)->
          getEntityShape(apf::Mesh::TRIANGLE)->getLocalGradients(m,e,xi,grads);

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        apf::NewArray<double> values;
        getBezier(3,P)->
            getEntityShape(apf::Mesh::TRIANGLE)->getValues(m,e,xi,values);

        for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i){
          values[i] *= weights[i];
          grads[i] = grads[i]*weights[i];
        }
        double sum = 0.;
        apf::Vector3 gsum = apf::Vector3(0,0,0);
        for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i){
          gsum += grads[i];
          sum += values[i];
        }
        for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i){
          grads[i] = (grads[i]*sum-gsum*values[i])/sum/sum;
        }
      }
    }
    int countNodes() const
    {
      return curved_face_total[BEZIER][P-1];
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
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
    void setWeights(apf::NewArray<double>& w)
    {
      weights.allocate(curved_face_total[BEZIER][P-1]);
      for(int i = 0; i < curved_face_total[BEZIER][P-1]; ++i)
        weights[i] = w[i];
    }
  private:
    apf::NewArray<double> weights;
    apf::NewArray<int> map[6];
  };
  apf::EntityShape* getEntityShape(int type)
  {
    static BezierShape::Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    static BezierShape::Tetrahedron tet;
    static apf::EntityShape* shapes[apf::Mesh::TYPES] =
    {&vertex,   //vertex
     &edge,    //edge
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
    if (dimension == 3)
      return false;
    else
      return true;
  }
  int countNodesOn(int type)
  {
    return getBezier(3,P)->countNodesOn(type);
  }
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    getBezierNodeXi(type,P,node,xi);
  }
  int getOrder() {return std::max(P,getBlendingOrder());}
  void setEdgeWeights(apf::NewArray<double>& weights)
  {
    Edge* edge = static_cast<Edge*>(getEntityShape(apf::Mesh::EDGE));
    edge->setWeights(weights);
  }
  void setTriangleWeights(apf::NewArray<double>& weights)
  {
    Triangle* triangle =
        static_cast<Triangle*>(getEntityShape(apf::Mesh::TRIANGLE));
    triangle->setWeights(weights);
  }
  protected:
    std::string name;
};

static void setOrder(const int order)
{
  P = order;
}

apf::FieldShape* getBezier(int dimension, int order)
{
  crv::setOrder(order);

  static crv::BezierCurve bezierCurve;
  static crv::BezierSurface bezierSurface;

  if(dimension == 2)
    return &bezierCurve;
  else
    return &bezierSurface;
}

// Third order is a modified bezier.
// Technically a Bezier but points set using
// maCurveMesh G1 methodology
// Exists for debugging purposes

apf::FieldShape* getGregory(int order)
{
  static GregorySurface4 gregorySurface;

  if(order == 4){
    setOrder(order);
    return &gregorySurface;
  } if(order == 3)
    return getBezier(3,order);

  return NULL;
}

apf::FieldShape* getNurbs(int order)
{
  setOrder(order);
  static Nurbs nurbs;
  return &nurbs;
}

void setNurbsEdgeWeights(apf::NewArray<double>& weights)
{
  Nurbs* nurbs = dynamic_cast<Nurbs*>(getNurbs(P));
  nurbs->setEdgeWeights(weights);
}

void setNurbsTriangleWeights(apf::NewArray<double>& weights)
{
  Nurbs* nurbs = static_cast<Nurbs*>(getNurbs(P));
  nurbs->setTriangleWeights(weights);
}

} // namespace crv
