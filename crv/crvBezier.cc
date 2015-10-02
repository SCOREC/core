/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <math.h>
#include "crvBezier.h"
#include "crvTables.h"
/* see bezier.tex */

namespace crv {

static int P = 1;

static bool useBlend()
{
  return (getBlendingOrder() != 0);
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
    ss << "Bezier" ;
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
        values[i+1] = binomial(P,i)*Bij(P-i,i,1.-t,t);
      values[0] = pow(1-t, P);
      values[1] = pow(t, P);

    }
    void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      double t = 0.5*(xi[0]+1.);
      grads.allocate(P+1);
      for(int i = 1; i < P; ++i)
        grads[i+1] = apf::Vector3(binomial(P,i)*(i-P*t)
            *Bij(P-1-i,i-1,1.-t,t)/2.,0,0);
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
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate((P+1)*(P+2)/2);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      if(!useBlend() || m->getModelType(m->toModel(e)) != m->getDimension()){

        for(int i = 0; i < P+1; ++i)
          for(int j = 0; j < P+1-i; ++j)
            values[getTriNodeIndex(P,i,j)] =
                trinomial(P,i,j)*Bijk(i,j,P-i-j,xii[0],xii[1],xii[2]);

      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate((P+1)*(P+2)/2);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

      if(!useBlend() || m->getModelType(m->toModel(e)) != m->getDimension()){
        for(int i = 0; i < 3; ++i)
          grads[i] = gxii[i]*P*pow(xii[i],P-1);

        for(int i = 1; i < P+1; ++i)
          for(int j = 1; j < P-i; ++j)
            grads[getTriNodeIndex(P,i,j)] =
                gxii[0]*trinomial(P,i,j)*(i*(1.-xii[1])-(P-j)*xii[0])
                *Bijk(i-1,j,P-i-j-1,xii[0],xii[1],xii[2]) +
                gxii[1]*trinomial(P,i,j)*(j*(1.-xii[0])-(P-i)*xii[1])
                *Bijk(i,j-1,P-i-j-1,xii[0],xii[1],xii[2]);

        // i = 0
        for(int j = 1; j < P; ++j)
          grads[3+2*(P-1)-j] =
            (gxii[0]*(j-P)*xii[1] +
            gxii[1]*(j*(1.-xii[0])-P*xii[1]))
            *binomial(P,j)*Bij(j-1,P-j-1,xii[1],xii[2]);

        // j = 0
        for(int i = 1; i < P; ++i)
          grads[3+2*(P-1)-1+i] =
              (gxii[0]*(i*(1.-xii[1])-P*xii[0]) +
               gxii[1]*(i-P)*xii[0])
              *binomial(P,i)*Bij(i-1,P-i-1,xii[0],xii[2]);

        // k = 0
        for(int i = 1, j = P-1; i < P; ++i, --j)
          grads[3+(P-1)-i] =
              (gxii[0]*i*xii[1] + gxii[1]*j*xii[0])
              *trinomial(P,i,j)*Bij(i-1,j-1,xii[0],xii[1]);

      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

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
      if(!useBlend()){
        values.allocate((P+1)*(P+2)*(P+3)/6);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        for(int i = 0; i < 4; ++i)
          values[i] = pow(xii[i],P);

        int nE = P-1;

        int const (*tev)[2] = apf::tet_edge_verts;

        for(int a = 0; a < 6; ++a)
          for(int b = 0; b < nE; ++b) // edge nodes
            values[4+a*nE+b] = binomial(P,b+1)
              *Bij(P-b-1,b+1,xii[tev[a][0]],xii[tev[a][1]]);

        // face 0, l = 0
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j)
              values[computeTetNodeIndex(P,i,j,P-i-j)] = trinomial(P,i,j)
                  *Bijk(i,j,P-i-j,xii[0],xii[1],xii[2]);
        // face 1, k = 0
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j)
              values[computeTetNodeIndex(P,i,j,0)] = trinomial(P,i,j)
                  *Bijk(i,j,P-i-j,xii[0],xii[1],xii[3]);
        // face 2, i = 0
        for(int j = 1; j <= P-1; ++j)
          for(int k = 1; k <= P-1-j; ++k)
              values[computeTetNodeIndex(P,0,j,k)] = trinomial(P,j,k)
                  *Bijk(j,k,P-j-k,xii[1],xii[2],xii[3]);
        // face 3, j = 0
        for(int i = 1; i <= P-1; ++i)
            for(int k = 1; k <= P-1-i; ++k)
                values[computeTetNodeIndex(P,i,0,k)] = trinomial(P,i,k)
                    *Bijk(i,k,P-i-k,xii[0],xii[2],xii[3]);

        // internal nodes
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j)
            for(int k = 1; k <= P-1-i-j; ++k)
              values[computeTetNodeIndex(P,i,j,k)] = quadnomial(P,i,j,k)
                  *Bijkl(i,j,k,P-i-j-k,xii[0],xii[1],xii[2],xii[3]);


      } else {
        values.allocate(2*P*P+2);
        BlendedTetGetValues(m,e,xi,values);
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      if(!useBlend()){
        grads.allocate((P+1)*(P+2)*(P+3)/6);

        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
            apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

        for(int i = 0; i < 4; ++i)
          grads[i] = gxii[i]*P*pow(xii[i],P-1);

        int nE = P-1;

        int const (*tev)[2] = apf::tet_edge_verts;
        int const (*ttv)[3] = apf::tet_tri_verts;

        for(int a = 0; a < 6; ++a)
          for(int b = 0; b < nE; ++b) // edge nodes
            grads[4+a*nE+b] = gxii[tev[a][0]]*binomial(P,b+1)*(P-b-1)
                              *Bij(P-b-2,b+1,xii[tev[a][0]],xii[tev[a][1]])
                            + gxii[tev[a][1]]*binomial(P,b+1)*(b+1)
                              *Bij(P-b-1,b,xii[tev[a][0]],xii[tev[a][1]]);

        // face 0, l = 0
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j){
            int index = computeTetNodeIndex(P,i,j,P-i-j);
            grads[index].zero();
            int ijk[3] = {i,j,P-i-j};
            for(int b = 0; b < 3; ++b){
              grads[index] += gxii[ttv[0][b]]*ijk[b]
                             *Bijk(ijk[b % 3]-1,ijk[(b+1) % 3],ijk[(b+2) % 3],
                              xii[ttv[0][b % 3]],xii[ttv[0][(b+1) % 3]],
                              xii[ttv[0][(b+2) % 3]]);

            }
            grads[index] = grads[index]*trinomial(P,i,j);
          }
        // face 1, k = 0
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j){
            int index = computeTetNodeIndex(P,i,j,0);
            grads[index].zero();
            int ijk[3] = {i,j,P-i-j};
            for(int b = 0; b < 3; ++b){
              grads[index] += gxii[ttv[1][b]]*ijk[b]
                             *Bijk(ijk[b % 3]-1,ijk[(b+1) % 3],ijk[(b+2) % 3],
                              xii[ttv[1][b % 3]],xii[ttv[1][(b+1) % 3]],
                              xii[ttv[1][(b+2) % 3]]);

            }
            grads[index] = grads[index]*trinomial(P,i,j);
          }
        // face 2, i = 0
        for(int j = 1; j <= P-1; ++j)
          for(int k = 1; k <= P-1-j; ++k){
            int index = computeTetNodeIndex(P,0,j,k);
            grads[index].zero();
            int jkl[3] = {j,k,P-j-k};
            for(int b = 0; b < 3; ++b){
              grads[index] += gxii[ttv[2][b]]*jkl[b]
                             *Bijk(jkl[b % 3]-1,jkl[(b+1) % 3],jkl[(b+2) % 3],
                              xii[ttv[2][b % 3]],xii[ttv[2][(b+1) % 3]],
                              xii[ttv[2][(b+2) % 3]]);

            }
            grads[index] = grads[index]*trinomial(P,j,k);
          }
        // face 3, j = 0
        for(int i = 1; i <= P-1; ++i)
            for(int k = 1; k <= P-1-i; ++k){
              int index = computeTetNodeIndex(P,i,0,k);
              grads[index].zero();
              int ikl[3] = {i,k,P-i-k};
              for(int b = 0; b < 3; ++b){
                grads[index] += gxii[ttv[3][b]]*ikl[b]
                               *Bijk(ikl[b % 3]-1,ikl[(b+1) % 3],ikl[(b+2) % 3],
                                xii[ttv[3][b % 3]],xii[ttv[3][(b+1) % 3]],
                                xii[ttv[3][(b+2) % 3]]);

              }
              grads[index] = grads[index]*trinomial(P,i,k);
            }
        // internal nodes
        for(int i = 1; i <= P-1; ++i)
          for(int j = 1; j <= P-1-i; ++j)
            for(int k = 1; k <= P-1-i-j; ++k){
              int index = computeTetNodeIndex(P,i,j,k);
              grads[index].zero();
              int ijkl[4] = {i,j,k,P-i-j-k};
              for(int b = 0; b < 4; ++b){
                grads[index] += gxii[b]*ijkl[b]
                   *Bijkl(ijkl[b % 4]-1,ijkl[(b+1) % 4],ijkl[(b+2) % 4],ijkl[(b+3) % 4],
                       xii[b],xii[(b+1) % 4],xii[(b+2) % 4],xii[(b+3) % 4]);

              }
              grads[index] = grads[index]*quadnomial(P,i,j,k);
            }
      } else {
        grads.allocate(2*P*P+2);
        BlendedTetGetLocalGradients(m,e,xi,grads);
      }
    }
    int countNodes() const {
      if(!useBlend())
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
    if ((dimension < P && dimension < 3) || (P > 3 && !useBlend()))
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
        if(!useBlend()){
          return (P-1)*(P-2)*(P-3)/6;
        } else
          return 0;
      default:
        return 0;
    }
  }
  int getOrder() {return std::max(P,getBlendingOrder());}
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    getBezierNodeXi(type,P,node,xi);
  }
  protected:
    std::string name;
};

static int g_index[3][2] = {{2,1},{0,2},{1,0}};
static int g_pairs[3][2] = {{0,5},{1,3},{2,4}};

class GregorySurface3 : public apf::FieldShape
{
public:
  const char* getName() const {return name.c_str();}
  GregorySurface3() {
    std::stringstream ss;
    ss << "GregorySurface3";
    name = ss.str();
    this->registerSelf(name.c_str());
  }
  class Triangle : public apf::EntityShape
  {
  public:
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(15);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      apf::ModelEntity* g = m->toModel(e);
      if (!useBlend() || m->getModelType(g) != m->getDimension()){
        apf::NewArray<double> bvalues;
        getBezier(3)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getValues(m,e,xi,bvalues);

        for(int i = 0; i < 9; ++i)
          values[i] = bvalues[i];
        for(int i = 9; i < 15; ++i)
          values[i] = 0.;
        double bernstein = bvalues[9];
        for(int i = 0; i < 3; ++i){
          double xiix,x = xii[g_index[i][0]] + xii[g_index[i][1]];

          if(x < 1e-12)
            xiix = 0.5;
          else
            xiix = xii[g_index[i][0]]/x;

          values[9+g_pairs[i][0]] = bernstein*xiix*xii[i];
          values[9+g_pairs[i][1]] = bernstein*(1.-xiix)*xii[i];
        }

      } else
        BlendedTriangleGetValues(m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(15);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};
      apf::ModelEntity* g = m->toModel(e);
      if (!useBlend() || m->getModelType(g) != m->getDimension()){
        apf::NewArray<apf::Vector3> bgrads;
        apf::NewArray<double> values;
        getBezier(3)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getLocalGradients(m,e,xi,bgrads);
        getBezier(3)->getEntityShape(apf::Mesh::TRIANGLE)
            ->getValues(m,e,xi,values);

        for(int i = 0; i < 9; ++i)
          grads[i] = bgrads[i];
        for(int i = 9; i < 15; ++i)
          grads[i].zero();
        double x, xiix;

        apf::Vector3 xv, gx;
        apf::Vector3 bernsteinGrad = bgrads[9];
        double bernstein = values[9];

        for(int i = 0; i < 3; ++i){
          x  = xii[g_index[i][0]] + xii[g_index[i][1]];
          gx = gxii[g_index[i][0]] + gxii[g_index[i][1]];

          if(x < 1e-12){
            xiix = 0.5;
            xv.zero();
          } else {
            xiix = xii[g_index[i][0]]/x;
            xv = (gxii[g_index[i][0]]*xii[g_index[i][1]]
                 -gxii[g_index[i][1]]*xii[g_index[i][0]])
                 *bernstein/x/x*xii[i];
          }

          grads[9+g_pairs[i][0]] = bernsteinGrad*xiix*xii[i] + xv + gxii[i]*bernstein*xiix;
          grads[9+g_pairs[i][1]] = bernsteinGrad*(1.-xiix)*xii[i] - xv + gxii[i]*bernstein*(1.-xiix);

        }
      } else
        BlendedTriangleGetLocalGradients(m,e,xi,grads);

    }
    int countNodes() const
    {
      return 15;
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
      if(!useBlend()){
        values.allocate(40);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        for(int i = 0; i < 4; ++i)
          values[i] = pow(xii[i],P);
        for(int i = 4; i < 40; ++i)
          values[i] = 0.;
        int nE = 2;
        int nF = 15;

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

            values[4+6*nE+a*nF+g_pairs[b][0]] = bernstein*xiix*xii[ttv[a][b]];
            values[4+6*nE+a*nF+g_pairs[b][1]] = bernstein*(1.-xiix)*xii[ttv[a][b]];
          }
        }
      } else {
        values.allocate(40);
        BlendedTetGetValues(m,e,xi,values);
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      if(!useBlend()){
        grads.allocate(40);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
            apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

        for(int i = 0; i < 4; ++i)
          grads[i] = gxii[i]*P*pow(xii[i],P-1);
        for(int i = 4; i < 40; ++i)
          grads[i].zero();
        int nE = 2;
        int nF = 15;

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

            if(x < 1e-12){
              xiix = 0.5;
            } else {
              xiix = xii[ttv[a][g_index[b][0]]]/x;
            }

            for(int c = 0; c < 3; ++c){
              int ijk[3] = {1,1,1};
              ijk[b] += P-3; ijk[c] -= 1;
              grads[4+6*nE+a*nF+g_pairs[b][0]] += gxii[ttv[a][c]]
                  *trinomial(P,P-2,1)
                  *(ijk[c]+1)*Bijk(ijk,xi_t)*xiix*xii[ttv[a][b]];
              grads[4+6*nE+a*nF+g_pairs[b][1]] += gxii[ttv[a][c]]
                  *trinomial(P,P-2,1)
                  *(ijk[c]+1)*Bijk(ijk,xi_t)*(1.-xiix)*xii[ttv[a][b]];
            }
          }
        } // done faces
      } else {
        grads.allocate(40);
        BlendedTetGetLocalGradients(m,e,xi,grads);
      }
    }
    int countNodes() const {
      return 40;
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      int which,rotate;
      bool flip;
      getAlignment(m,elem,shared,which,flip,rotate);
      if(m->getType(shared) == apf::Mesh::EDGE){
        if(!flip)
          for(int i = 0; i < 2; ++i)
            order[i] = i;
        else
          for(int i = 0; i < 2; ++i)
            order[i] = 1-i;
        return;
      }
      static int orients[6][6] =
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
    if (dimension == 3)
      return false;
    else
      return true;
  }
  int countNodesOn(int type)
  {
    static int nodes[apf::Mesh::TYPES] =
    {1,                 //vertex
     2,                 //edge
     6,                 //triangle
     0,                 //quad
     0,       //tet
     0,                 //hex
     0,                 //prism
     0};                //pyramid
    return nodes[type];
  }
  /* These don't make sense for gregory patches */
  void getNodeXi(int type, int node, apf::Vector3& xi)
  {
    static double edgePoints[2] = {-0.4503914,0.4503914};
    if (type == apf::Mesh::EDGE) {
      xi[0] = edgePoints[node];
    } else if (type == apf::Mesh::TRIANGLE) {
      xi[0] = 1./3.;
      xi[1] = 1./3.;
    }
    else
      xi.zero();
  }
  int getOrder() {return std::max(3,getBlendingOrder());}
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
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(18);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      apf::ModelEntity* g = m->toModel(e);
      if (!useBlend() || m->getModelType(g) != m->getDimension()){
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
      apf::ModelEntity* g = m->toModel(e);
      if (!useBlend() || m->getModelType(g) != m->getDimension()){
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
      if(!useBlend()){
        values.allocate(47);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        for(int i = 0; i < 4; ++i)
          values[i] = pow(xii[i],P);

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
      if(!useBlend()){
        grads.allocate(47);
        double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
        apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
            apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

        for(int i = 0; i < 4; ++i)
          grads[i] = gxii[i]*P*pow(xii[i],P-1);

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
    int countNodes() const {
      if(!useBlend())
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
    if (useBlend() && dimension == 3)
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
     !useBlend(),       //tet
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
      getBezier(P)->getEntityShape(apf::Mesh::EDGE)->getValues(0,0,xi,values);

      double sum = 0.;
      for(int i = 0; i < P+1; ++i){
        values[i] *= weights[i];
        sum += values[i];
      }
      for(int i = 0; i < P+1; ++i)
        values[i] /= sum;
    }
    void getLocalGradients(apf::Mesh* /*m*/, apf::MeshEntity* /*e*/,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      apf::NewArray<double> values;
      getBezier(P)->
          getEntityShape(apf::Mesh::EDGE)->getLocalGradients(0,0,xi,grads);
      getBezier(P)->getEntityShape(apf::Mesh::EDGE)->getValues(0,0,xi,values);

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
      for(int i = 0; i < P+1; ++i)
        grads[i] = (grads[i]*sum-gsum*values[i])/sum/sum;
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
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      getBezier(P)->
          getEntityShape(apf::Mesh::TRIANGLE)->getValues(m,e,xi,values);

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){

        double sum = 0.;
        for(int i = 0; i < (P+1)*(P+2)/2; ++i){
          values[i] *= weights[i];
          sum += values[i];
        }
        for(int i = 0; i < (P+1)*(P+2)/2; ++i){
          values[i]/=sum;
        }
      }
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      getBezier(P)->
          getEntityShape(apf::Mesh::TRIANGLE)->getLocalGradients(m,e,xi,grads);

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        apf::NewArray<double> values;
        getBezier(P)->
            getEntityShape(apf::Mesh::TRIANGLE)->getValues(m,e,xi,values);
        int n = (P+1)*(P+2)/2;
        for(int i = 0; i < n; ++i){
          values[i] *= weights[i];
          grads[i] = grads[i]*weights[i];
        }
        double sum = 0.;
        apf::Vector3 gsum = apf::Vector3(0,0,0);
        for(int i = 0; i < n; ++i){
          gsum += grads[i];
          sum += values[i];
        }
        for(int i = 0; i < n; ++i){
          grads[i] = (grads[i]*sum-gsum*values[i])/sum/sum;
        }
      }
    }
    int countNodes() const
    {
      return (P+1)*(P+2)/2;
    }
    void alignSharedNodes(apf::Mesh* m,
        apf::MeshEntity* elem, apf::MeshEntity* shared, int order[])
    {
      alignEdgeWithTri(m,elem,shared,order);
    }
    void setWeights(apf::NewArray<double>& w)
    {
      weights.allocate((P+1)*(P+2)/2);
      for(int i = 0; i < (P+1)*(P+2)/2; ++i)
        weights[i] = w[i];
    }
  private:
    apf::NewArray<double> weights;
  };
  apf::EntityShape* getEntityShape(int type)
  {
    static Bezier::Vertex vertex;
    static Edge edge;
    static Triangle triangle;
    static Bezier::Tetrahedron tet;
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
    return getBezier(P)->countNodesOn(type);
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

apf::FieldShape* getBezier(int order)
{
  if(order < 1 ||order > 19)
    fail("order must be in [1,19]\n");
  crv::setOrder(order);
  static crv::Bezier bezier;
  return &bezier;
}

apf::FieldShape* getGregory(int order)
{
  static GregorySurface3 gregorySurface3;
  static GregorySurface4 gregorySurface4;
  if(order == 3){
     setOrder(order);
     return &gregorySurface3;
  }
  if(order == 4){
    setOrder(order);
    return &gregorySurface4;
  }
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
