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

static double const blendingTol = 1.e-12;

/* default is no blending */
static int B = 0;

void setBlendingOrder(const int b)
{
  B = b;
}

int getBlendingOrder()
{
  return B;
}

void BlendedTriangleGetValues(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values)
{
  // Triangular Blending
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  for(int i = 0; i < 3; ++i)
    values[i] = -pow(xii[i],B);
  // zero the rest, the face node weight is always zero
  int n = m->getShape()->getEntityShape(apf::Mesh::TRIANGLE)->countNodes();
  for(int i = 3; i < n; ++i)
    values[i] = 0.0;

  double x, xiix;
  apf::Vector3 xv;
  apf::NewArray<double> v;

  int const (*tev)[2] = apf::tri_edge_verts;
  int nE = m->getShape()->countNodesOn(apf::Mesh::EDGE);

  apf::MeshEntity* edges[3];
  m->getDownward(e,1,edges);
  for(int i = 0; i < 3; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];

    if(x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix-1.0;
    m->getShape()->getEntityShape(apf::Mesh::EDGE)
            ->getValues(m,edges[i],xv,v);

    for(int j = 0; j < 2; ++j)
      values[tev[i][j]]   += v[j]*pow(x,B);
    for(int j = 0; j < nE; ++j)
      values[3+i*nE+j] = v[2+j]*pow(x,B);
  }
}

void BlendedTriangleGetLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads)
{

  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};
  apf::Vector3 gxii[3] = {apf::Vector3(-1,-1,0),
      apf::Vector3(1,0,0),apf::Vector3(0,1,0)};

  for(int i = 0; i < 3; ++i)
    grads[i] = gxii[i]*-pow(xii[i],B-1.)*B;

  int n = m->getShape()->getEntityShape(apf::Mesh::TRIANGLE)->countNodes();
  for(int i = 3; i < n; ++i)
    grads[i].zero();

  double x, xiix;
  apf::Vector3 xv, gx;

  apf::NewArray<double> v;
  apf::NewArray<apf::Vector3> gv;
  apf::MeshEntity* edges[3];

  int const (*tev)[2] = apf::tri_edge_verts;
  int nE = m->getShape()->countNodesOn(apf::Mesh::EDGE);

  m->getDownward(e,1,edges);
  for(int i = 0; i < 3; ++i){
    x = xii[tev[i][0]]+xii[tev[i][1]];
    if(x < blendingTol)
      xiix = 0.5;
    else
      xiix = xii[tev[i][1]]/x;

    xv[0] = 2.0*xiix-1.0;

    gx = gxii[tev[i][0]]+gxii[tev[i][1]];

    m->getShape()->getEntityShape(apf::Mesh::EDGE)
              ->getValues(m,edges[i],xv,v);

    m->getShape()->getEntityShape(apf::Mesh::EDGE)
              ->getLocalGradients(m,edges[i],xv,gv);

    for(int j = 0; j < 2; ++j)
      grads[tev[i][j]]   += gx*B*pow(x,B-1.)*v[j]
        + (gxii[tev[i][1]]-gx*xiix)*gv[j][0]*2.*pow(x,B-1.);

    for(int j = 0; j < nE; ++j)
      grads[3+i*nE+j] = gx*B*pow(x,B-1.)*v[j+2]
        + (gxii[tev[i][1]]-gx*xiix)*gv[j+2][0]*2.*pow(x,B-1.);
  }
}

void BlendedTetGetValues(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<double>& values)
{
  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};

  for(int i = 0; i < 4; ++i)
    values[i] = pow(xii[i],B);

  int n = m->getShape()->getEntityShape(apf::Mesh::TET)->countNodes();
  for(int i = 4; i < n; ++i)
    values[i] = 0.0;

  double x;
  apf::Vector3 xv, xiix;
  apf::NewArray<double> v;

  apf::MeshEntity* faces[4];

  int nE = m->getShape()->countNodesOn(apf::Mesh::EDGE);
  int nF = m->getShape()->countNodesOn(apf::Mesh::TRIANGLE);

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

    m->getShape()->getEntityShape(apf::Mesh::EDGE)
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

    m->getShape()->getEntityShape(apf::Mesh::TRIANGLE)
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

void BlendedTetGetLocalGradients(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3 const& xi, apf::NewArray<apf::Vector3>& grads)
{
  double xii[4] = {1.-xi[0]-xi[1]-xi[2],xi[0],xi[1],xi[2]};
  apf::Vector3 gxii[4] = {apf::Vector3(-1,-1,-1),apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),apf::Vector3(0,0,1)};

  for(int i = 0; i < 4; ++i)
    grads[i] = gxii[i]*pow(xii[i],B-1.)*B;

  int n = m->getShape()->getEntityShape(apf::Mesh::TET)->countNodes();
  for(int i = 4; i < n; ++i)
    grads[i].zero();

  double x;
  apf::Vector3 xv, xiix, gx;
  apf::Vector3 gxv[2];
  apf::NewArray<double> v;
  apf::NewArray<apf::Vector3> gv;

  apf::MeshEntity* faces[4];

  int nE = m->getShape()->countNodesOn(apf::Mesh::EDGE);
  int nF = m->getShape()->countNodesOn(apf::Mesh::TRIANGLE);

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

    m->getShape()->getEntityShape(apf::Mesh::EDGE)
        ->getValues(m,edges[i],xv,v);
    m->getShape()->getEntityShape(apf::Mesh::EDGE)
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

    m->getShape()->getEntityShape(apf::Mesh::TRIANGLE)
          ->getValues(m,faces[i],xv,v);
    m->getShape()->getEntityShape(apf::Mesh::TRIANGLE)
          ->getLocalGradients(m,faces[i],xv,gv);


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

}
