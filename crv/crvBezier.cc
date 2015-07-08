/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include <math.h>

#include "crv.h"

/* see bezier.tex */

namespace crv {
// negative -> flipped relative to canonical
// relies on e0 being always ordered correctly
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

enum {
  CURVED_BEZIER,
  CURVED_GREGORY,
  CURVED_TYPES
};

// numbers of nodes on
static int const curved_face_internal[2][6] =
{{0,0,1,3,6,10},{0,0,0,6,0,0}};

// total numbers of nodes
static int const curved_face_total[2][6] =
{{3,6,10,15,21,28},{0,0,0,18,0,0}};

static int const blended_tet_total[2][6] =
{{4,10,20,34,52,74},{0,0,0,46,0,0}};

static double const blendingTol = 1.e-12;

static double B = 2.;

static int P = 0;

static void BlendedTriangleGetValues(const int type,
    apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values)
{
  // Triangular Blending
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  for(int i = 0; i < 3; ++i)
    values[i] = -pow(xii[i],B);
  // zero the rest, the face node weight is always zero
  for(int i = 3; i < curved_face_total[type][P-1]; ++i)
    values[i] = 0.0;

  double x, xiix;
  apf::Vector3 xv;
  apf::NewArray<double> v;

  int const (*tev)[2] = apf::tri_edge_verts;

  apf::MeshEntity* edges[3];
  m->getDownward(e,1,edges);
  for(int i = 0; i < 3; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];

    if(x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix-1.0;

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
          ->getValues(m,edges[i],xv,v);

    for(int j = 0; j < 2; ++j)
      values[tev[i][j]]   += v[j]*pow(x,B);
    for(int j = 0; j < (P-1); ++j)
      values[3+i*(P-1)+j] = v[2+j]*pow(x,B);
  }
}

static void BlendedTriangleGetLocalGradients(const int type,
    apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads)
{

  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
  apf::Vector3 gxii[3] = {apf::Vector3(-1,-1,0),
      apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

  for(int i = 0; i < 3; ++i)
    grads[i] = gxii[i]*-pow(xii[i],B-1.)*B;

  for(int i = 3; i < curved_face_total[type][P-1]; ++i)
    grads[i].zero();

  double x, xiix;
  apf::Vector3 xv, gx;

  apf::NewArray<double> v;
  apf::NewArray<apf::Vector3> gv;
  apf::MeshEntity* edges[3];

  int const (*tev)[2] = apf::tri_edge_verts;

  m->getDownward(e,1,edges);
  for(int i = 0; i < 3; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];
    if(x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix-1.0;

    gx = gxii[tev[i][0]]+gxii[tev[i][1]];

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
      ->getValues(m,edges[i],xv,v);

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
      ->getLocalGradients(m,edges[i],xv,gv);

    for(int j = 0; j < 2; ++j)
      grads[tev[i][j]]   += gx*B*pow(x,B-1.)*v[j]
        + (gxii[tev[i][1]]-gx*xiix)*gv[j][0]*2.*pow(x,B-1.);

    for(int j = 0; j < (P-1); ++j)
      grads[3+i*(P-1)+j] = gx*B*pow(x,B-1.)*v[j+2]
        + (gxii[tev[i][1]]-gx*xiix)*gv[j+2][0]*2.*pow(x,B-1.);

  }

}

static void BlendedTetrahedronGetValues(const int type,
    apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values)
{

  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};

  for(int i = 0; i < 4; ++i)
    values[i] = pow(xii[i],B);

  for(int i = 4; i < blended_tet_total[type][P-1]; ++i)
    values[i] = 0.0;

  double x;
  apf::Vector3 xv, xiix;
  apf::NewArray<double> v;

  apf::MeshEntity* faces[4];

  int nE = P-1;
  int nF = curved_face_internal[type][P-1];

  int const (*tev)[2] = apf::tet_edge_verts;
  int const (*ttv)[3] = apf::tet_tri_verts;

  apf::MeshEntity* edges[6];
  m->getDownward(e,1,edges);
  for(int i = 0; i < 6; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];
    if(x < blendingTol)
      xiix[0] = 0.5;
    else
      xiix[0] = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix[0]-1.0;

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
      ->getValues(m,edges[i],xv,v);

    for(int j = 0; j < 2; ++j) // vertices
      values[tev[i][j]] += -v[j]*pow(x,B);
    for(int j = 0; j < nE; ++j)// edge nodes
      values[4+i*nE+j]  += -v[2+j]*pow(x,B);

  }

  m->getDownward(e,2,faces);
  for(int i = 0; i < 4; ++i){

    x = 0.;
    for(int j = 0; j < 3; ++j)
      x += xii[ttv[i][j]];

    if(x < blendingTol)
      xv = apf::Vector3(1./3.,1./3.,1./3.);
    else {
      for(int j = 0; j < 3; ++j)
        xv[j] = xii[ttv[i][(j+1) % 3]]/x;
    }

    if(type == CURVED_BEZIER)
      getBezier(3,P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);
    else if(type == CURVED_GREGORY)
      getGregory(P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);

    for(int j = 0; j < 3; ++j) // vertices
      values[ttv[i][j]] += v[j]*pow(x,B);

    // Edge contributions from faces
    // Edges are the first 3*nE entries per face
    // also not necessarily aligned, canonically
    for(int k = 0; k < 3; ++k){
      int l = tet_tri_edges[i][k];
      if(flip_tet_tri_edges[i][k] == false){
        for(int j = 0; j < nE; ++j)
          values[4+l*nE+j] += v[3+k*nE+j]*pow(x,B);
      } else { // we need to flip
        for(int j = 0; j < nE; ++j)
          values[4+l*nE+j] += v[3+k*nE+nE-1-j]*pow(x,B);
      }
    }
    for(int j = 0; j < nF; ++j) // face nodes
      values[4+6*nE+i*nF+j] +=  v[3+3*nE+j]*pow(x,B);
  } // done faces

}

static void BlendedTetrahedronGetLocalGradients(const int type,
    apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads)
{
  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
  apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

  for(int i = 0; i < 4; ++i)
    grads[i] = gxii[i]*pow(xii[i],B-1.)*B;

  for(int i = 4; i < blended_tet_total[type][P-1]; ++i)
    grads[i].zero();

  double x;
  apf::Vector3 xv, xiix, gx;
  apf::Vector3 gxv[2];
  apf::NewArray<double> v;
  apf::NewArray<apf::Vector3> gv;

  apf::MeshEntity* faces[4];

  int nE = P-1;
  int nF = curved_face_internal[type][P-1];

  int const (*tev)[2] = apf::tet_edge_verts;
  int const (*ttv)[3] = apf::tet_tri_verts;

  apf::MeshEntity* edges[6];
  m->getDownward(e,1,edges);
  for(int i = 0; i < 6; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];
    if(x < blendingTol)
      xiix[0] = 0.5;
    else
      xiix[0] = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix[0]-1.0;
    gx = gxii[tev[i][0]]+gxii[tev[i][1]];

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
        ->getValues(m,edges[i],xv,v);

    getBezier(3,P,B)->getEntityShape(apf::Mesh::EDGE)
        ->getLocalGradients(m,edges[i],xv,gv);

    for(int j = 0; j < 2; ++j) // vertices
      grads[tev[i][j]] += gx*B*pow(x,B-1.)*(-v[j])
      - (gxii[tev[i][1]]-gx*xiix[0])*gv[j][0]*2.*pow(x,B-1.);

    for(int j = 0; j < nE; ++j)// edge nodes
      grads[4+i*nE+j]  += gx*B*pow(x,B-1.)*(-v[j+2])
      - (gxii[tev[i][1]]-gx*xiix[0])*gv[j+2][0]*2.*pow(x,B-1.);
  }

  m->getDownward(e,2,faces);
  for(int i = 0; i < 4; ++i){
    x = 0.;
    for(int j = 0; j < 3; ++j)
      x += xii[ttv[i][j]];

    if(x < blendingTol)
      xv = apf::Vector3(1./3.,1./3.,1./3.);
    else {
      for(int j = 0; j < 3; ++j)
        xv[j] = xii[ttv[i][(j+1) % 3]]/x;
    }

    gx = gxii[ttv[i][0]] + gxii[ttv[i][1]] + gxii[ttv[i][2]];

    gxv[0] = (gxii[ttv[i][1]]-gx*xv[0])*pow(x,B-1.);
    gxv[1] = (gxii[ttv[i][2]]-gx*xv[1])*pow(x,B-1.);

    if(type == CURVED_BEZIER){
      getBezier(3,P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);
      getBezier(3,P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getLocalGradients(m,faces[i],xv,gv);
    }
    else if(type == CURVED_GREGORY){
      getGregory(P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);
      getGregory(P,B)->getEntityShape(apf::Mesh::TRIANGLE)
          ->getLocalGradients(m,faces[i],xv,gv);
    }

    for(int j = 0; j < 3; ++j) // vertices
      grads[ttv[i][j]] += gx*B*pow(x,B-1.)*v[j]
        + gxv[0]*gv[j][0] + gxv[1]*gv[j][1];

    for(int k = 0; k < 3; ++k){
      int l = tet_tri_edges[i][k];
      if(flip_tet_tri_edges[i][k] == false){
        for(int j = 0; j < nE; ++j)
          grads[4+l*nE+j] += gx*B*pow(x,B-1.)*v[3+k*nE+j]
            + gxv[0]*gv[3+k*nE+j][0]
            + gxv[1]*gv[3+k*nE+j][1];

      } else { // we need to flip
        for(int j = 0; j < nE; ++j)
          grads[4+l*nE+j] += gx*B*pow(x,B-1.)*v[3+k*nE+nE-1-j]
            + gxv[0]*gv[3+k*nE+nE-1-j][0]
            + gxv[1]*gv[3+k*nE+nE-1-j][1];
      }
    }
    for(int j = 0; j < nF; ++j) // face nodes
      grads[4+6*nE+i*nF+j] += gx*B*pow(x,B-1.)*v[3+3*nE+j]
        + gxv[0]*gv[3+3*nE+j][0]
        + gxv[1]*gv[3+3*nE+j][1];
  } // done faces

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
        map[j-1].allocate(curved_face_total[CURVED_BEZIER][j-1]);             
        for(int i = 0; i < curved_face_total[CURVED_BEZIER][j-1]; ++i)
          map[j-1][i] = maps[j-1][i];
      }
    }
    void getValues(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<double>& values) const
    {
      values.allocate(curved_face_total[CURVED_BEZIER][P-1]);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){

        for(int i = 0; i < P+1; ++i)
          for(int j = 0; j < P+1-i; ++j)
            values[map[P-1][j*(P+1)+i-j*(j-1)/2]] =
                binomial(P,i)*binomial(P-i,j)
                *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],P-i-j);

      } else
        BlendedTriangleGetValues(CURVED_BEZIER,m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e, apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(curved_face_total[CURVED_BEZIER][P-1]);

      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
      apf::Vector3 gxii[3] =
        {apf::Vector3(-1,-1,0),apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

      apf::ModelEntity* g = m->toModel(e);

      if (m->getModelType(g) != m->getDimension()){
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
        grads[map[P-1][0]] = (gxii[0]+gxii[1])*(-P)*pow(xii[2],P-1);
        // i = k = 0
        grads[map[P-1][((P+1)*(P+2))/2-1]] = gxii[1]*P*pow(xii[1],P-1);
        // j = k = 0
        grads[map[P-1][P]] = gxii[0]*P*pow(xii[0],P-1);

      } else
        BlendedTriangleGetLocalGradients(CURVED_BEZIER,m,e,xi,grads);

    }
    int countNodes() const {return curved_face_total[CURVED_BEZIER][P-1];}
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
      values.allocate(blended_tet_total[CURVED_BEZIER][P-1]);
      BlendedTetrahedronGetValues(CURVED_BEZIER,m,e,xi,values);
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(blended_tet_total[CURVED_BEZIER][P-1]);
      BlendedTetrahedronGetLocalGradients(CURVED_BEZIER,m,e,xi,grads);
    }
    int countNodes() const {return blended_tet_total[CURVED_BEZIER][P-1];}
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
      int n = curved_face_internal[CURVED_BEZIER][P-1];
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
        ((dimension == 2) && P > 2))
      return true;
    else
      return false;
  }
  int countNodesOn(int type)
  {
    switch (type) {
      case 0:
        return 1;
      case 1:
        return P-1;
      case 2:
        return curved_face_internal[CURVED_BEZIER][P-1];
      default:
        return 0;
    }
  }
  int getOrder() {return std::max(P,(int)B);}
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
      xi.zero();
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
    static double eP2[1] = {0.0};
    static double eP3[2] = {-0.4503914,0.4503914};
    static double eP4[3] = {-0.6612048,0.0,0.6612048};
    static double eP5[4] = {-0.7732854,-0.2863522,0.2863522,0.7732854};
    static double eP6[5] = {-0.8388042,-0.469821,0.0,
      0.469821,0.8388042};
    static double* edgePoints[6] =
    {eP2, eP2, eP3, eP4, eP5, eP6 };
    if (type == apf::Mesh::EDGE) {
      xi[0] = edgePoints[P-1][node];
    } else if (type == apf::Mesh::TRIANGLE) {
      xi = apf::Vector3(1./3.,1./3.,1./3.);
      if(node == curved_face_internal[CURVED_BEZIER][P-1]-1 && P % 3 == 0){
        return;
      } else { // technically only two of these numbers are needed
        switch (P) {
          case 1:
          case 2:
          case 3:
            fail("expected P >= 4");
          case 4:
            xi[(node+2) % 3] = 0.5582239;
            xi[(node+0) % 3] = 0.22088805;
            xi[(node+1) % 3] = 0.22088805;
            break;
          case 5:
            if(node % 2 == 0) {
              xi[(node/2+2) % 3] = 0.6949657;
              xi[(node/2+0) % 3] = 0.15251715;
              xi[(node/2+1) % 3] = 0.15251715;
            } else {
              xi[((node-1)/2+2) % 3] = 0.4168658;
              xi[((node-1)/2+0) % 3] = 0.4168658;
              xi[((node-1)/2+1) % 3] = 0.1662684;
            }
            break;
          case 6:
            if (node % 3 == 0) {
              xi[(node/3+2) % 3] = 0.7805723;
              xi[(node/3+0) % 3] = 0.10971385;
              xi[(node/3+1) % 3] = 0.10971385;
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
      xi.zero();
    }
  }
protected:
  std::string name;
};

static void getBezierCurveInterPtsToCtrlPts(apf::NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -0.970273514083553,0.333333333333333,2.71895067382449,-1.08201049307427,
      0.333333333333333,-0.970273514083553,-1.08201049307427,2.71895067382449};
  double e4[15] = {
      -1.4304202857228,-0.25,3.39545839723405,-1.46967987431139,
      0.754641762800137,0.953613523815196,0.953613523815197,-2.76673344002279,4.62623983241519,
      -2.76673344002279,-0.25,-1.4304202857228,0.754641762800137,-1.46967987431139,
      3.39545839723405};
  double e5[24] = {
      -1.88592269024942,0.2,4.05614415979432,-1.81653638123435,
      1.0295423816296,-0.583227469940158,1.85476912333284,-0.942961345124708,-5.01939997635205,6.96205913930752,
      -4.56234099538677,2.70787405422317,-0.942961345124708,1.85476912333285,2.70787405422317,-4.56234099538677,
      6.96205913930752,-5.01939997635206,0.2,-1.88592269024942,-0.583227469940158,1.0295423816296,
      -1.81653638123435,4.05614415979432};
  double e6[35] = {
      -2.33890800235808,-0.166666666666667,4.70907763497668,-2.14695478588352,
      1.2670886004356,-0.80040589915343,0.476769118649422,3.03457283388393,0.935563200943235,-7.82909978199834,9.74813267975089,
      -6.60581336123903,4.37362214799981,-2.65697771934049,-2.27592962541295,-2.27592962541295,6.3088040999163,-9.70710791530195,
      12.3484668815972,-9.70710791530194,6.3088040999163,0.935563200943235,3.03457283388393,-2.65697771934049,4.37362214799981,
      -6.60581336123903,9.74813267975088,-7.82909978199834,-0.166666666666667,-2.33890800235809,0.476769118649422,-0.80040589915343,
      1.2670886004356,-2.14695478588352,4.70907763497668};
  double* table[5] = {
      e2,e3,e4,e5,e6};
  int nb = P-1;
  int ni = P+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[P-2][i*ni+j];

}
static void getBezierShapeInterPtsToCtrlPts(int type,
    apf::NewArray<double> & c)
{
  double e2[3] = {-0.5,-0.5,2};
  double e3[8] = {
      -1.00596379148431,0.333333333333333,2.69317845753742,-1.02054799938644,
      0.333333333333333,-1.00596379148431,-1.02054799938644,2.69317845753742};
  double e4[15] = {
      -1.52680420766155,-0.25,3.37567603341243,-1.28732567375034,
      0.688453847999458,1.01786947177436,1.01786947177436,-2.70941992094126,4.38310089833379,
      -2.70941992094126,-0.25,-1.52680420766155,0.688453847999459,-1.28732567375034,
      3.37567603341243};
  double e5[24] = {
      -2.06136018481524,0.2,4.05936188849008,-1.52513018373641,
      0.846118038541145,-0.518989558479574,2.07392858296555,-1.03068009240762,-5.04403528329688,6.4850808350761,
      -4.1256786562572,2.64138461392004,-1.03068009240762,2.07392858296555,2.64138461392004,-4.1256786562572,
      6.48508083507611,-5.04403528329688,0.2,-2.06136018481524,-0.518989558479574,0.846118038541145,
      -1.52513018373641,4.05936188849008};
  double e6[35] = {
      -2.60465921875445,-0.166666666666667,4.74317259573776,-1.74847211573074,
      0.99151406014263,-0.630691218757935,0.415802564029398,3.50945493261743,1.04186368750178,-8.00777336834505,8.99109434126308,
      -5.80152934540333,3.85156704410589,-2.5846772917398,-2.63209119946307,-2.63209119946307,6.39664544713349,-8.91824703868012,
      11.3073855820194,-8.91824703868012,6.39664544713349,1.04186368750178,3.50945493261743,-2.5846772917398,3.85156704410589,
      -5.80152934540333,8.99109434126308,-8.00777336834505,-0.166666666666667,-2.60465921875445,0.415802564029398,-0.630691218757935,
      0.991514060142631,-1.74847211573074,4.74317259573777};
  double f3[10] = {
      0.5059637914843115,0.5059637914843116,0.5059637914843119,-0.8363152290754889,
      -0.8363152290754891,-0.8363152290754891,-0.8363152290754897,-0.8363152290754891,
      -0.8363152290754898,4.5};
  double f4[45] = {
      1.473866405971784,-0.4873245458152482,-0.4873245458152483,-2.157170326583087,
      -0.7895371825786641,1.002268609355065,0.4569922360188257,0.4160673296503333,
      0.4569922360188257,1.002268609355065,-0.7895371825786643,-2.157170326583087,
      7.066481928037422,-2.00343662222666,-2.00343662222666,
      -0.4873245458152483,1.473866405971784,-0.4873245458152482,1.002268609355065,
      -0.7895371825786646,-2.157170326583087,-2.157170326583087,-0.7895371825786643,
      1.002268609355065,0.4569922360188255,0.4160673296503333,0.4569922360188255,
      -2.00343662222666,7.066481928037421,-2.00343662222666,
      -0.4873245458152483,-0.4873245458152481,1.473866405971784,0.4569922360188258,
      0.4160673296503332,0.4569922360188255,1.002268609355065,-0.7895371825786646,
      -2.157170326583087,-2.157170326583088,-0.7895371825786643,1.002268609355065,
      -2.00343662222666,-2.00343662222666,7.066481928037422};

  double f5[126] = {
      2.955509375394112,0.4850633218816851,0.4850633218816851,-4.015153971419451,
      -0.6339185304220427,1.140973028173454,-1.041278651590119,-0.3192965376118003,
      -0.2178043751582093,-0.2178043751582093,-0.3192965376118008,-1.041278651590118,
      1.140973028173453,-0.633918530422042,-4.015153971419451,10.16896081482358,
      -3.086925777444234,1.177507607441262,0.8971975820812148,1.177507607441263,
      -3.086925777444236,
      -1.87982956469617,-1.879829564696171,0.5695483481139566,3.558924203285685,
      -2.437779912677446,-2.437779912677446,3.558924203285686,1.772590081128005,
      0.9751122788728792,0.1149068856949735,-0.833741865064019,-0.8337418650640189,
      0.1149068856949739,0.9751122788728801,1.772590081128005,-6.137619043989807,
      12.47852594090006,-6.137619043989807,-2.280366067778754,2.247531721435293,
      -2.280366067778758,
      0.4850633218816852,2.955509375394111,0.4850633218816852,-1.041278651590119,
      1.140973028173453,-0.6339185304220422,-4.01515397141945,-4.015153971419449,
      -0.6339185304220425,1.140973028173453,-1.041278651590118,-0.3192965376118005,
      -0.2178043751582091,-0.2178043751582096,-0.3192965376118004,1.177507607441262,
      -3.086925777444234,10.16896081482358,-3.086925777444233,1.177507607441262,
      0.8971975820812153,
      0.569548348113957,-1.87982956469617,-1.879829564696171,-0.8337418650640194,
      0.114906885694974,0.97511227887288,1.772590081128004,3.558924203285685,
      -2.437779912677446,-2.437779912677443,3.558924203285685,1.772590081128006,
      0.9751122788728793,0.1149068856949746,-0.8337418650640199,2.247531721435295,
      -2.280366067778758,-6.137619043989805,12.47852594090005,-6.137619043989806,
      -2.280366067778759,
      0.4850633218816858,0.4850633218816853,2.955509375394112,-0.3192965376118007,
      -0.2178043751582097,-0.2178043751582098,-0.3192965376118001,-1.041278651590119,
      1.140973028173454,-0.6339185304220429,-4.015153971419452,-4.01515397141945,
      -0.6339185304220436,1.140973028173455,-1.04127865159012,1.177507607441263,
      0.8971975820812174,1.177507607441262,-3.086925777444234,10.16896081482359,
      -3.086925777444237,
      -1.879829564696172,0.5695483481139572,-1.879829564696172,1.772590081128005,
      0.9751122788728803,0.114906885694975,-0.8337418650640201,-0.8337418650640197,
      0.114906885694974,0.9751122788728801,1.772590081128007,3.558924203285685,
      -2.437779912677443,-2.437779912677449,3.558924203285688,-6.137619043989808,
      -2.280366067778761,2.247531721435295,-2.280366067778758,-6.137619043989812,
      12.4785259409000};
  double f6[280] = {
      4.990795106388393,-0.4798837984147291,-0.4798837984147287,-6.458245423578343,
      -0.3639631629214212,1.22336563258444,-1.237281132017878,1.048637951819922,
      0.2399023612066681,0.1476405704559505,0.09371675625975406,0.1476405704559514,
      0.2399023612066667,1.048637951819923,-1.23728113201788,1.223365632584442,
      -0.3639631629214235,-6.458245423578343,13.89283735202935,-4.276280903149708,
      1.808916611395355,-0.8152389231864423,-0.4900563063939803,-0.4900563063939829,
      -0.8152389231864404,1.808916611395354,-4.276280903149705,1.327623829722839,
      -4.689458981957318,2.303409358216713,-0.6686820957216264,8.401259643873185,
      -5.327640035897593,-2.380733253738293,4.596086431435174,-4.67688716273589,
      -1.613859854887045,-0.8313011663168295,-0.3532695050699889,0.01692167342423888,
      0.7838078136607817,1.115432125134806,-0.6396258416621832,-0.4063403163660048,
      1.633507742829366,4.420179509817743,-13.10100871536221,20.06968509845678,
      -10.70127780752898,5.23777942309397,2.406990467529688,0.7247013815918657,
      -2.116258714513936,3.071067981878374,-2.093285349716039,-4.18119984946874,
      2.303409358216714,-4.689458981957322,-0.6686820957216258,-4.676887162735897,
      4.596086431435186,-2.380733253738301,-5.327640035897589,8.401259643873187,
      4.420179509817745,1.633507742829365,-0.4063403163660045,-0.6396258416621829,
      1.115432125134805,0.7838078136607807,0.01692167342423913,-0.353269505069989,
      -0.8313011663168312,-1.613859854887043,5.23777942309397,-10.70127780752899,
      20.06968509845678,-13.10100871536222,-2.093285349716035,3.071067981878372,
      -2.116258714513934,0.7247013815918641,2.406990467529694,-4.181199849468742,
      -0.4798837984147298,4.990795106388394,-0.4798837984147283,1.048637951819924,
      -1.237281132017881,1.223365632584441,-0.3639631629214217,-6.458245423578346,
      -6.458245423578344,-0.3639631629214246,1.223365632584442,-1.237281132017879,
      1.048637951819922,0.2399023612066666,0.1476405704559507,0.09371675625975429,
      0.1476405704559513,0.2399023612066678,-0.8152389231864428,1.808916611395358,
      -4.27628090314971,13.89283735202935,-4.276280903149705,1.808916611395352,
      -0.8152389231864391,-0.4900563063939816,-0.490056306393983,1.32762382972284,
      -0.6686820957216246,-4.689458981957318,2.303409358216711,1.115432125134803,
      -0.6396258416621817,-0.4063403163660055,1.633507742829363,4.420179509817745,
      8.40125964387318,-5.32764003589758,-2.3807332537383,4.596086431435179,
      -4.676887162735889,-1.613859854887039,-0.8313011663168305,-0.3532695050699883,
      0.016921673424238,0.78380781366078,-2.116258714513931,3.071067981878368,
      -2.093285349716028,-13.10100871536221,20.06968509845676,-10.70127780752896,
      5.23777942309396,2.406990467529691,0.724701381591865,-4.181199849468744,
      -0.6686820957216236,2.30340935821671,-4.689458981957321,0.7838078136607789,
      0.01692167342423904,-0.353269505069988,-0.8313011663168306,-1.61385985488704,
      -4.676887162735886,4.596086431435173,-2.380733253738294,-5.327640035897591,
      8.401259643873187,4.420179509817743,1.633507742829364,-0.4063403163660066,
      -0.6396258416621801,1.115432125134801,-2.116258714513929,0.7247013815918617,
      2.406990467529691,5.237779423093964,-10.70127780752897,20.06968509845676,
      -13.10100871536221,-2.093285349716029,3.071067981878367,-4.18119984946874,
      -0.4798837984147294,-0.4798837984147293,4.990795106388396,0.2399023612066681,
      0.1476405704559513,0.09371675625975392,0.147640570455951,0.2399023612066676,
      1.048637951819923,-1.23728113201788,1.223365632584443,-0.3639631629214256,
      -6.458245423578345,-6.458245423578348,-0.3639631629214194,1.223365632584439,
      -1.237281132017879,1.048637951819923,-0.8152389231864435,-0.4900563063939818,
      -0.4900563063939812,-0.8152389231864419,1.808916611395355,-4.276280903149707,
      13.89283735202936,-4.27628090314971,1.808916611395358,1.327623829722839,
      2.303409358216711,-0.6686820957216234,-4.689458981957324,-1.613859854887043,
      -0.8313011663168315,-0.353269505069988,0.01692167342423824,0.7838078136607798,
      1.115432125134801,-0.639625841662179,-0.4063403163660082,1.633507742829368,
      4.420179509817745,8.401259643873196,-5.327640035897597,-2.380733253738291,
      4.596086431435175,-4.676887162735889,5.23777942309397,2.406990467529692,
      0.7247013815918628,-2.11625871451393,3.071067981878368,-2.093285349716032,
      -13.10100871536222,20.06968509845678,-10.70127780752898,-4.181199849468741,
      -4.689458981957318,-0.6686820957216246,2.303409358216713,4.420179509817742,
      1.633507742829366,-0.4063403163660054,-0.6396258416621825,1.115432125134803,
      0.7838078136607793,0.01692167342423931,-0.3532695050699873,-0.831301166316834,
      -1.61385985488704,-4.676887162735893,4.596086431435182,-2.380733253738301,
      -5.327640035897584,8.401259643873182,-13.10100871536221,-2.093285349716037,
      3.071067981878372,-2.11625871451393,0.7247013815918595,2.406990467529697,
      5.237779423093967,-10.70127780752898,20.06968509845677,-4.181199849468742,
      2.740410900541109,2.740410900541107,2.740410900541116,-3.896719679725134,
      0.8525687056196629,2.629196981978395,0.8525687056196668,-3.896719679725134,
      -3.89671967972513,0.8525687056196584,2.6291969819784,0.8525687056196664,
      -3.896719679725144,-3.896719679725143,0.8525687056196644,2.629196981978399,
      0.8525687056196618,-3.896719679725134,9.218530840490745,-7.999447648758338,
      -7.999447648758347,9.218530840490741,-7.99944764875833,-7.999447648758361,
      9.218530840490766,-7.999447648758353,-7.99944764875834,23.49717556815212
      };

  double* table[10] = {
      e2,e3,e4,e5,e6,NULL,f3,f4,f5,f6};
  int nb = (type == apf::Mesh::TRIANGLE) ?
      curved_face_internal[CURVED_BEZIER][P-1] : P-1;
  int ni = (type == apf::Mesh::TRIANGLE) ?
      curved_face_total[CURVED_BEZIER][P-1] : P+1;
  c.allocate(ni*nb);
  for( int i = 0; i < nb; ++i)
    for( int j = 0; j < ni; ++j)
      c[i*ni+j] = table[5*(type-1)+(P-2)][i*ni+j];
}

/* end of Bezier, begining of Gregory Surface */

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
      values.allocate(curved_face_total[CURVED_GREGORY][3]);
      double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

      apf::ModelEntity* g = m->toModel(e);
      if (m->getModelType(g) != m->getDimension()){
        for(int i = 0; i < 5; ++i)
          for(int j = 0; j < 5-i; ++j)
            values[map[j*5+i-j*(j-1)/2]] =
                binomial(4,i)*binomial(4-i,j)
                *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],4-i-j);

        for(int i = 0; i < 3; ++i){
          double x = xii[index[i][0]] + xii[index[i][1]];
          double bernstein = values[12+pairs[i][0]];
          values[12+pairs[i][1]] = 0.;

          if(x < blendingTol) continue;
          values[12+pairs[i][0]] = bernstein*xii[index[i][0]]/x;
          values[12+pairs[i][1]] = bernstein*xii[index[i][1]]/x;
        }
      } else
        BlendedTriangleGetValues(CURVED_GREGORY,m,e,xi,values);

    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(curved_face_total[CURVED_GREGORY][3]);
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

        double x;
        apf::Vector3 xv, gx;
        apf::NewArray<double> v;

        getValues(m,e,xi,v);

        for(int i = 0; i < 3; ++i){
          x  = xii[index[i][0]] + xii[index[i][1]];
          gx = gxii[index[i][0]] + gxii[index[i][1]];

          apf::Vector3 bernstein = grads[12+pairs[i][0]];
          grads[12+pairs[i][1]].zero();

          if(x < blendingTol) continue;

          grads[12+pairs[i][0]] = bernstein*xii[index[i][0]]/x;
          grads[12+pairs[i][1]] = bernstein*xii[index[i][1]]/x;

        }
      } else
        BlendedTriangleGetLocalGradients(CURVED_GREGORY,m,e,xi,grads);

    }
    int countNodes() const
    {
      return curved_face_total[CURVED_GREGORY][3];
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
      values.allocate(blended_tet_total[CURVED_GREGORY][3]);
      BlendedTetrahedronGetValues(CURVED_GREGORY,m,e,xi,values);
    }
    void getLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
        apf::Vector3 const& xi,
        apf::NewArray<apf::Vector3>& grads) const
    {
      grads.allocate(blended_tet_total[CURVED_GREGORY][3]);
      BlendedTetrahedronGetLocalGradients(CURVED_GREGORY,m,e,xi,grads);
    }
    int countNodes() const {return blended_tet_total[CURVED_GREGORY][3];}
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
  int getOrder() {return std::max(4,(int)B);}
protected:
  std::string name;
};

static void setOrder(const int order, const int blendOrder)
{
  P = order;
  B = blendOrder;
}

apf::FieldShape* getBezier(int dimension, int order, int blendOrder)
{
  crv::setOrder(order,blendOrder);

  static crv::BezierCurve bezierCurve;
  static crv::BezierSurface bezierSurface;

  //  bezier::setOrder(order,blendOrder);

  if(dimension == 2)
    return &bezierCurve;
  else
    return &bezierSurface;
}

// Third order is a modified bezier.
// Technically a Bezier but points set using
// maCurveMesh G1 methodology
// Exists for debugging purposes

apf::FieldShape* getGregory(int order, int blendOrder)
{
  static GregorySurface4 gregorySurface;

  if(order == 4){
    setOrder(order, blendOrder);
    return &gregorySurface;
  } if(order == 3)
    return getBezier(3,order,blendOrder);

  return NULL;
}

void getTransformationCoefficients(int dim, int type,
    apf::NewArray<double>& c){
  if(dim == 2)
    getBezierCurveInterPtsToCtrlPts(c);
  else
    getBezierShapeInterPtsToCtrlPts(type,c);
}

} // namespace crv
