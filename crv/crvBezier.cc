/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <math.h>
#include "crvBezier.h"
#include "crvBezierShapes.h"
#include "crvMath.h"
#include "crvShape.h"
#include "crvTables.h"

namespace crv {

/* The order used is static and set here. Rather than
   define a large number of classes for each order,
   since shape functions are static, we choose not to template either,
   but rather force the order of the mesh to be constant,
   when variable order shape functions are permitted, this
   will need to be significantly changed

   For functionality such as elevating the mesh order,
   this needs to be changed frequently, as permitting multiple
   instances different orders of bezier shape functions is not possible
 */
static int P = 1;

static bool useBlending(int type)
{
  return (getBlendingOrder(type) != 0);
}

void getFullRepFromBlended(int type,
    apf::NewArray<double>& transformCoefficients,
    apf::NewArray<apf::Vector3>& elemNodes)
{
  int n = getNumControlPoints(type,P);
  int ne = getNumInternalControlPoints(type,P);
  apf::NewArray<apf::Vector3> newNodes(ne);
  apf::NewArray<apf::Vector3> nodes(n-ne);
  convertInterpolationPoints(n-ne,ne,elemNodes,transformCoefficients,newNodes);
  for (int i = 0; i < n-ne; ++i)
    nodes[i] = elemNodes[i];

  elemNodes.allocate(n);
  for (int i = 0; i < n-ne; ++i)
    elemNodes[i] = nodes[i];
  for (int i = 0; i < ne; ++i)
    elemNodes[n-ne+i] = newNodes[i];
}


static void alignEdgeWithTri(apf::Mesh* m, apf::MeshEntity* elem,
    apf::MeshEntity* shared, int order[])
{
  int which,rotate;
  bool flip;
  getAlignment(m,elem,shared,which,flip,rotate);
  if(!flip)
    for(int i = 0; i < P-1; ++i)
      order[i] = i;
  else
    for(int i = 0; i < P-1; ++i)
      order[i] = P-2-i;
  return;
}

class Bezier : public apf::FieldShape
{
public:
  const char* getName() const {return name.c_str();}
  Bezier() {
    std::stringstream ss;
    ss << "Bezier";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
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
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Bezier shapes");
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
      values.allocate(P+1);
      bezier[apf::Mesh::EDGE](P,xi,values);
    }
    void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(P+1);
      bezierGrads[apf::Mesh::EDGE](P,xi,grads);
    }
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Bezier shapes");
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
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate((P+1)*(P+2)/2);

      if(!useBlending(apf::Mesh::TRIANGLE)
          || isBoundaryEntity(m,e)){
        bezier[apf::Mesh::TRIANGLE](P,xi,values);
      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate((P+1)*(P+2)/2);

      if(!useBlending(apf::Mesh::TRIANGLE)
          || isBoundaryEntity(m,e)){
        bezierGrads[apf::Mesh::TRIANGLE](P,xi,grads);
      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

    }
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Bezier shapes");
    }
    int countNodes() const {return getNumControlPoints(apf::Mesh::TRIANGLE,P);}
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      alignEdgeWithTri(m,elem,shared,order);
    }
  };
  class Tetrahedron : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      if(!useBlending(apf::Mesh::TET)){
        values.allocate((P+1)*(P+2)*(P+3)/6);
        bezier[apf::Mesh::TET](P,xi,values);
      } else {
        values.allocate(2*P*P+2);
        BlendedTetGetValues(m,e,xi,values);
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      if(!useBlending(apf::Mesh::TET)){
        grads.allocate((P+1)*(P+2)*(P+3)/6);
        bezierGrads[apf::Mesh::TET](P,xi,grads);
      } else {
        grads.allocate(2*P*P+2);
        BlendedTetGetLocalGradients(m,e,xi,grads);
      }
    }
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Bezier shapes");
    }
    int countNodes() const {
      if(!useBlending(apf::Mesh::TET))
        return (P+1)*(P+2)*(P+3)/6;
      else
        return 2*P*P+2;
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
      int n = (P-1)*(P-2)/2;
      if(P <= 6)
        for(int i = 0; i < n; ++i)
          order[i] = tet_tri[P][flip][rotate][i];
      else {
        // since we don't have a table, compute it all
        int index0, index1;
        if(!flip){
          index0 = (3-rotate) % 3;
          index1 = (4-rotate) % 3;
        } else {
          index0 = (rotate+2) % 3;
          index1 = (rotate+1) % 3;
        }
        int index = 0;
        for(int i = 0; i <= P-3; ++i)
          for(int j = 0; j <= P-3-i; ++j){
            int ijk[3] = {i,j,P-3-i-j};
            order[index] = ijk[index0]*(P-2)-ijk[index0]*(ijk[index0]-1)/2
              +ijk[index1];
            index++;
          }
      }
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
    if ((dimension < P && dimension < 3)
        || (P > 3 && !useBlending(apf::Mesh::TET)))
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
        return (P-1)*(P-2)/2;
      case apf::Mesh::TET:
        if(!useBlending(apf::Mesh::TET)){
          return (P-1)*(P-2)*(P-3)/6;
        } else
          return 0;
      default:
        return 0;
    }
  }
  int getOrder() {return P;}
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    getBezierNodeXi(type,P,node,xi);
  }
  protected:
    std::string name;
};

static int g_index[3][2] = {{2,1},{0,2},{1,0}};
static int g_pairs[3][2] = {{0,5},{1,3},{2,4}};

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
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(18);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      if (!useBlending(apf::Mesh::TRIANGLE)
          || isBoundaryEntity(m,e)){
        apf::NewArray<double> bvalues;
        getBezier(4)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getValues(m,e,xi,bvalues);

        for(int i = 0; i < 15; ++i)
          values[i] = bvalues[i];

        for(int i = 0; i < 3; ++i){
          double xiix,x = xii[g_index[i][0]] + xii[g_index[i][1]];
          double bernstein = values[12+g_pairs[i][0]];
          values[12+g_pairs[i][1]] = 0.;

          if(x < 1e-12)
            xiix = 0.5;
          else
            xiix = xii[g_index[i][0]]/x;
          values[12+g_pairs[i][0]] = bernstein*xiix;
          values[12+g_pairs[i][1]] = bernstein*(1.-xiix);
        }
      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(18);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

      if (!useBlending(apf::Mesh::TRIANGLE)
          || isBoundaryEntity(m,e)){
        apf::NewArray<apf::Vector3> bgrads;
        apf::NewArray<double> bvalues;

        getBezier(4)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getLocalGradients(m,e,xi,bgrads);
        getBezier(4)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getValues(m,e,xi,bvalues);

        for(int i = 0; i < 15; ++i)
          grads[i] = bgrads[i];

        double x, xiix;
        apf::Vector3 xv, gx;

        for(int i = 0; i < 3; ++i){
          x  = xii[g_index[i][0]] + xii[g_index[i][1]];
          gx = gxii[g_index[i][0]] + gxii[g_index[i][1]];

          apf::Vector3 bernsteinGrad = grads[12+g_pairs[i][0]];
          double bernstein = bvalues[12+g_pairs[i][0]];
          if(x < 1e-12){
            xiix = 0.5;
            xv.zero();
          } else {
            xiix = xii[g_index[i][0]]/x;
            xv = (gxii[g_index[i][0]]*xii[g_index[i][1]]
                 -gxii[g_index[i][1]]*xii[g_index[i][0]])
                 *bernstein/x/x;
          }

          grads[12+g_pairs[i][0]] = bernsteinGrad*xiix + xv;
          grads[12+g_pairs[i][1]] = bernsteinGrad*(1.-xiix) - xv;
        }
      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

    }
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Gregory shapes");
    }
    int countNodes() const
    {
      return 18;
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      alignEdgeWithTri(m,elem,shared,order);
    }
  };
  class Tetrahedron : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<double>& values) const
    {
      if(!useBlending(apf::Mesh::TET)){
        values.allocate(47);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        for(int i = 0; i < 4; ++i)
          values[i] = intpow(xii[i],P);

        int nE = 3;
        int nF = 6;
        int nT = 1;

        int const (*tev)[2] = apf::tet_edge_verts;
        int const (*ttv)[3] = apf::tet_tri_verts;

        for(int a = 0; a < 6; ++a)
          for(int b = 0; b < nE; ++b) // edge nodes
            values[4+a*nE+b] = binomial(P,b+1)
              *Bij(P-b-1,b+1,xii[tev[a][0]],xii[tev[a][1]]);

        for(int a = 0; a < 4; ++a){
          for(int b = 0; b < 3; ++b){ // face nodes
            double bernstein = trinomial(P,P-2,1)
              *Bijk(P-2,1,1,xii[ttv[a][b]],xii[ttv[a][(b+1) % 3]],
                  xii[ttv[a][(b+2) % 3]]);

            double xiix;
            double x = xii[ttv[a][g_index[b][0]]]
                     + xii[ttv[a][g_index[b][1]]];

            if(x < 1e-12)
              xiix = 0.5;
            else
              xiix = xii[ttv[a][g_index[b][0]]]/x;

            values[4+6*nE+a*nF+g_pairs[b][0]] = bernstein*xiix;
            values[4+6*nE+a*nF+g_pairs[b][1]] = bernstein*(1.-xiix);
          }
        }

        for(int a = 0; a < nT; ++a){ // internal nodes
          int ijkl[4] = {1,1,1,1};
          ijkl[a] += P-4;
          values[4+6*nE+4*nF+a] = quadnomial(P,P-3,1,1)
              *Bijkl(ijkl,xii);
        }
      } else {
        values.allocate(46);
        BlendedTetGetValues(m,e,xi,values);
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      if(!useBlending(apf::Mesh::TET)){
        grads.allocate(47);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
            apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

        for(int i = 0; i < 4; ++i)
          grads[i] = gxii[i]*P*intpow(xii[i],P-1);

        int nE = 3;
        int nF = 6;
        int nT = 1;
        apf::Vector3 xv;

        int const (*tev)[2] = apf::tet_edge_verts;
        int const (*ttv)[3] = apf::tet_tri_verts;

        for(int a = 0; a < 6; ++a)
          for(int b = 0; b < nE; ++b) // edge nodes
            grads[4+a*nE+b] = gxii[tev[a][0]]*binomial(P,b+1)*(P-b-1)
                              *Bij(P-b-2,b+1,xii[tev[a][0]],xii[tev[a][1]])
                            + gxii[tev[a][1]]*binomial(P,b+1)*(b+1)
                              *Bij(P-b-1,b,xii[tev[a][0]],xii[tev[a][1]]);

        for(int a = 0; a < 4; ++a){
          for(int b = 0; b < 3; ++b){ // face nodes
            double xi_t[3] = {xii[ttv[a][0]],xii[ttv[a][1]],xii[ttv[a][2]]};
            grads[4+6*nE+a*nF+g_pairs[b][0]].zero();
            grads[4+6*nE+a*nF+g_pairs[b][1]].zero();
            double xiix;
            double x = xii[ttv[a][g_index[b][0]]]
                     + xii[ttv[a][g_index[b][1]]];

            double bernstein = trinomial(P,P-2,1)
              *Bijk(P-2,1,1,xii[ttv[a][b]],xii[ttv[a][(b+1) % 3]],
                  xii[ttv[a][(b+2) % 3]]);

            if(x < 1e-12){
              xiix = 0.5;
              xv.zero();
            } else {
              xiix = xii[ttv[a][g_index[b][0]]]/x;
              xv = (gxii[ttv[a][g_index[b][0]]]*xii[ttv[a][g_index[b][1]]]
                   -gxii[ttv[a][g_index[b][1]]]*xii[ttv[a][g_index[b][0]]])
                   *bernstein/x/x;
            }

            for(int c = 0; c < 3; ++c){
              int ijk[3] = {1,1,1};
              ijk[b] += P-3; ijk[c] -= 1;
              grads[4+6*nE+a*nF+g_pairs[b][0]] += gxii[ttv[a][c]]
                  *trinomial(P,P-2,1)
                  *(ijk[c]+1)*Bijk(ijk,xi_t)*xiix + xv;
              grads[4+6*nE+a*nF+g_pairs[b][1]] += gxii[ttv[a][c]]
                  *trinomial(P,P-2,1)
                  *(ijk[c]+1)*Bijk(ijk,xi_t)*(1.-xiix) - xv;
            }
          }
        } // done faces

        for(int a = 0; a < nT; ++a){
          grads[4+6*nE+4*nF+a].zero();
          for(int b = 0; b < 4; ++b)
            grads[4+6*nE+4*nF+a] += gxii[b]*quadnomial(P,P-3,1,1)*(P-3)
                *Bijkl(P-4,1,1,1,xii[b],xii[(b+1) % 4],
                xii[(b+2) % 4],xii[(b+3) % 4]);
        }
      } else {
        grads.allocate(46);
        BlendedTetGetLocalGradients(m,e,xi,grads);
      }
    }
    void getVectorValues(apf::Mesh*, apf::MeshEntity*,
        apf::Vector3 const&, apf::NewArray<apf::Vector3>&) const
    {
      fail("getVectorValues is not implemented for Gregory shapes");
    }
    int countNodes() const {
      if(!useBlending(apf::Mesh::TET))
        return 47;
      else
        return 46;
    }
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
    static Bezier::Vertex vertex;
    static Bezier::Edge edge;
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
    if (useBlending(apf::Mesh::TET) && dimension == 3)
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
     !useBlending(4),       //tet
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
        xi[(node+0) % 3] = 0.5582239;
        xi[(node+1) % 3] = 0.22088805;
        xi[(node+2) % 3] = 0.22088805;
      }
    } else if (type == apf::Mesh::TET)
      xi = apf::Vector3(0.25,0.25,0.25);
    else
      xi.zero();
  }
  int getOrder() {return 4;}
protected:
  std::string name;
};

void setOrder(const int order)
{
  P = order;
}
int getOrder()
{
  return P;
}

apf::FieldShape* getBezier(int order)
{
  if(order < 1 ||order > 19)
    fail("order must be in [1,19]\n");
  setOrder(order);
  static crv::Bezier bezier;
  return &bezier;
}

apf::FieldShape* getGregory()
{
  setOrder(4);
  static GregorySurface4 gregorySurface4;
  return &gregorySurface4;
}

} // namespace crv
