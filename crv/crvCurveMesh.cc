/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */
#include "crv.h"
#include "crvAdapt.h"
#include "crvBezier.h"
#include "crvShape.h"
#include "crvSnap.h"

#include <cassert>

namespace crv {

void convertInterpolationPoints(int n, int ne,
    apf::NewArray<apf::Vector3>& nodes,
    apf::NewArray<double>& c,
    apf::NewArray<apf::Vector3>& newNodes){

  for(int i = 0; i < ne; ++i)
    newNodes[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      newNodes[i] += nodes[j]*c[i*n+j];

}

void convertInterpolationPoints(apf::Mesh2* m, apf::MeshEntity* e,
    int n, int ne, apf::NewArray<double>& c){

  apf::NewArray<apf::Vector3> l, b(ne);
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  crv::convertInterpolationPoints(n,ne,l,c,b);

  for(int i = 0; i < ne; ++i)
    m->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}

void snapToInterpolate(apf::Mesh2* m, apf::MeshEntity* e)
{
  int type = m->getType(e);
  if(type == apf::Mesh::VERTEX){
    apf::Vector3 p, pt(0,0,0);
    apf::ModelEntity* g = m->toModel(e);
    m->getParamOn(g,e,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,0,pt);
    return;
  }
  apf::FieldShape * fs = m->getShape();
  int non = fs->countNodesOn(type);
  apf::Vector3 p, xi, pt(0,0,0);
  for(int i = 0; i < non; ++i){
    apf::ModelEntity* g = m->toModel(e);
    fs->getNodeXi(type,i,xi);
    if(type == apf::Mesh::EDGE)
      transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);
    else
      transferParametricOnTriSplit(m,e,xi,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,i,pt);
  }
}

void MeshCurver::synchronize()
{
  apf::synchronize(m_mesh->getCoordinateField());
}

void MeshCurver::snapToInterpolate(int dim)
{
  apf::MeshEntity* e;
  apf::MeshIterator* it = m_mesh->begin(dim);
  while ((e = m_mesh->iterate(it))) {
    if(isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
      crv::snapToInterpolate(m_mesh,e);
  }
  m_mesh->end(it);
}

bool InterpolatingCurver::run()
{
  // interpolate points in each dimension
  for(int d = 1; d < 2; ++d)
    snapToInterpolate(d);

  synchronize();

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

void BezierCurver::convertInterpolatingToBezier()
{
  apf::FieldShape * fs = m_mesh->getShape();
  int order = fs->getOrder();

  int md = m_mesh->getDimension();
  int blendingOrder = getBlendingOrder(apf::Mesh::simplexTypes[md]);
  // go downward, and convert interpolating to control points
  int startDim = md - (blendingOrder > 0);

  for(int d = startDim; d >= 1; --d){
    if(!fs->hasNodesIn(d)) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getBezierTransformationCoefficients(order,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n,ne,c);
    }
    m_mesh->end(it);
  }
  // if we have a full representation, we need to place internal nodes on
  // triangles and tetrahedra
  for(int d = 2; d <= md; ++d){
    if(!fs->hasNodesIn(d) ||
        getBlendingOrder(apf::Mesh::simplexTypes[d])) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<double> c;
    getInternalBezierTransformationCoefficients(m_mesh,order,1,
        apf::Mesh::simplexTypes[d],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(!isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n-ne,ne,c);
    }
    m_mesh->end(it);
  }

  synchronize();
}

bool BezierCurver::run()
{
  std::string name = m_mesh->getShape()->getName();
  if(m_order < 1 || m_order > 6){
    fail("trying to convert to unimplemented Bezier order\n");
  }
  // if the currentOrder is NOT one, assume the mesh
  // is already interpolated and curved as one would want it,
  // otherwise, change it down to first order and then go back up
  int currentOrder = m_mesh->getShape()->getOrder();

  // if its already bezier, check what needs to be done, if anything
  if(name == std::string("Bezier")){
    changeMeshOrder(m_mesh,m_order);
    return true;
  } else {
    // if the initial mesh is first order, project the points
    // onto the new shape
    if(currentOrder == 1 || m_order == 1 || currentOrder > m_order)
    {
      apf::changeMeshShape(m_mesh, getBezier(m_order),true);
    } else if(currentOrder > 1){
      // assume it is interpolating, don't project
      apf::changeMeshShape(m_mesh,getBezier(m_order),false);
    }
  }

  if (m_mesh->canSnap()){
    for(int d = 1; d <= 2; ++d)
      snapToInterpolate(d);
    synchronize();
  }

  convertInterpolatingToBezier();

  // curving 1D meshes, while rare, is important in testing
  // do not fix shape if this is the case
  // does not work for blended shapes, yet
  // comment out for now

//  if( m_mesh->getDimension() >= 2 && m_order > 1 && blendingOrder == 0){
//    ma::Input* shapeFixer = configureShapeCorrection(m_mesh);
//    crv::adapt(shapeFixer);
//  }
  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

void GregoryCurver::setCubicEdgePointsUsingNormals()
{
  apf::MeshEntity* e;
  apf::Vector3 p, xi, pt;
  apf::Vector3 points[4];

  apf::MeshIterator* it = m_mesh->begin(1);

  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);

    if(m_mesh->getModelType(g) == 3) continue;
    // set edges using normals
    if(m_mesh->getModelType(g) == 1) {

      apf::Vector3 t[2];
      apf::MeshEntity* v[2];
      m_mesh->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m_mesh->getPoint(v[i],0,points[i*3]);
        m_mesh->getParamOn(g,v[i],p);
        m_mesh->getFirstDerivative(g,p,t[i],t[i]);
        t[i] = t[i].normalize();
      }
      double d = (points[3]-points[0]).getLength();
      apf::Vector3 l = (points[3]-points[0])/d;
      points[1] = points[0] + t[0]*(l*t[0])/fabs(l*t[0])*d/3.;
      points[2] = points[3] - t[1]*(l*t[1])/fabs(l*t[1])*d/3.;
      for(int i = 0; i < 2; ++i)
        m_mesh->setPoint(e,i,points[i+1]);

    } else {
      // set edges using tangents
      apf::Vector3 n[2];
      apf::MeshEntity* v[2];
      m_mesh->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m_mesh->getPoint(v[i],0,points[i*3]);
        m_mesh->getParamOn(g,v[i],p);
        m_mesh->getNormal(g,p,n[i]);
      }
      double d = (points[3]-points[0]).getLength();
      apf::Vector3 l = (points[3]-points[0])/d;
      double a[3] = {n[0]*l,n[1]*l,n[0]*n[1]};

      double r = 6.*(2.*a[0]+a[2]*a[1])/(4.-a[2]*a[2]);
      double s = 6.*(2.*a[1]+a[2]*a[0])/(4.-a[2]*a[2]);
      points[1] = points[0] + (l*6. - n[0]*2.*r+ n[1]*s)   *d/18.;
      points[2] = points[3] - (l*6. + n[0]*r   - n[1]*2.*s)*d/18.;
      for(int i = 0; i < 2; ++i)
        m_mesh->setPoint(e,i,points[i+1]);

    }
  }
  m_mesh->end(it);
}

static void elevateBezierCurves(apf::Mesh2* m)
{

  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    if(isBoundaryEntity(m,e))
      elevateBezierCurve(m,e,3,1);
  }
  m->end(it);
}

void GregoryCurver::setInternalPointsLocally()
{
  apf::Vector3 D[3][4];
  apf::Vector3 W[3][3];
  apf::Vector3 A[3][3];

  double lam[3][2];
  double mu[3][2];
  apf::Vector3 G[6];

  apf::MeshEntity* e;
  apf::MeshIterator* it = m_mesh->begin(2);
  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);
    if(!m_mesh->isOwned(e) || m_mesh->getModelType(g) != 2) continue;

    apf::Vector3 n[3];
    apf::MeshEntity* verts[3];
    apf::MeshEntity* edges[3];
    m_mesh->getDownward(e,0,verts);
    m_mesh->getDownward(e,1,edges);

    // elevated edges
    apf::NewArray<apf::Vector3> q(12);

    for(int i = 0; i < 3; ++i){
      apf::Vector3 param;
      m_mesh->getPoint(verts[i],0,q[i]);
      m_mesh->getParamOn(g,verts[i],param);
      m_mesh->getNormal(g,param,n[i]);
    }

    // elevate the edge points without formally setting them to q
    // compute tangent vectors, W
    for(int i = 0; i < 3; ++i){
      apf::Element* edge =
          apf::createElement(m_mesh->getCoordinateField(),edges[i]);
      apf::NewArray<apf::Vector3> ep;
      apf::getVectorNodes(edge,ep);

      bool flip;
      int which, rotate;
      apf::getAlignment(m_mesh,e,edges[i],which,flip,rotate);

      if(flip){
        W[i][0] = ep[3]-ep[1];
        W[i][1] = ep[2]-ep[3];
        W[i][2] = ep[0]-ep[2];
        q[i*3+3] = ep[1]*0.25+ep[3]*0.75;
        q[i*3+4] = ep[3]*0.5+ep[2]*0.5;
        q[i*3+5] = ep[2]*0.75+ep[0]*0.25;
      } else {
        W[i][0] = ep[2]-ep[0];
        W[i][1] = ep[3]-ep[2];
        W[i][2] = ep[1]-ep[3];
        q[i*3+3] = ep[0]*0.25+ep[2]*0.75;
        q[i*3+4] = ep[2]*0.5+ep[3]*0.5;
        q[i*3+5] = ep[3]*0.75+ep[1]*0.25;
      }

      apf::destroyElement(edge);
    }
    int const (*tev)[2] = apf::tri_edge_verts;

    for(int i = 0; i < 3; ++i){
      A[i][0] = apf::cross(n[tev[i][0]],W[i][0].normalize());
      A[i][2] = apf::cross(n[tev[i][1]],W[i][2].normalize());
      A[i][1] = (A[i][0]+A[i][2]).normalize();
    }

    D[0][0] = q[11] - (q[0]+q[3] )*0.5;
    D[0][3] = q[6]  - (q[1]+q[5] )*0.5;

    D[1][0] = q[5]  - (q[1]+q[6] )*0.5;
    D[1][3] = q[9]  - (q[2]+q[8] )*0.5;

    D[2][0] = q[8]  - (q[2]+q[9] )*0.5;
    D[2][3] = q[3]  - (q[0]+q[11])*0.5;

    for(int i = 0; i < 3; ++i){
      lam[i][0] = D[i][0]*W[i][0]/(W[i][0]*W[i][0]);
      lam[i][1] = D[i][3]*W[i][2]/(W[i][2]*W[i][2]);

      mu[i][0]  = D[i][0]*A[i][0];
      mu[i][1]  = D[i][3]*A[i][2];
    }

    for(int i = 0; i < 3; ++i){
      G[i] = (q[i*3+3]+q[i*3+4])*0.5
          + W[i][1]*2./3.*lam[i][0] + W[i][0]*1./3.*lam[i][1]
          + A[i][1]*2./3.*mu[i][0] + A[i][0]*1./3.*mu[i][1];
      G[3+i] = (q[i*3+4]+q[i*3+5])*0.5
          + W[i][2]*1./3.*lam[i][0] + W[i][1]*2./3.*lam[i][1]
          + A[i][2]*1./3.*mu[i][0] + A[i][1]*2./3.*mu[i][1];
    }
    for(int i = 0; i < 6; ++i)
      m_mesh->setPoint(e,i,G[i]);

  }
  m_mesh->end(it);
}

bool GregoryCurver::run()
{
  if(m_order != 4){
    fail("cannot convert to G1 of this order\n");
  }
  if(m_mesh->getDimension() != 3)
    fail("can only convert to 3D mesh\n");

  apf::changeMeshShape(m_mesh, getGregory(),true);
  int md = m_mesh->getDimension();
  apf::FieldShape * fs = m_mesh->getShape();

  // interpolate points in each dimension
  for(int d = 1; d < 2; ++d)
    snapToInterpolate(d);

  synchronize();

  // go downward, and convert interpolating to control points
  for(int d = md; d >= 1; --d){
    if(!fs->hasNodesIn(d)) continue;

    int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
    apf::NewArray<apf::Vector3> l, b(ne);

    apf::NewArray<double> c;

    getGregoryTransformationCoefficients(apf::Mesh::simplexTypes[d],c);

    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);

    while ((e = m_mesh->iterate(it))) {
      if(m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n,ne,c);
    }
    m_mesh->end(it);
  }

  setCubicEdgePointsUsingNormals();
  setInternalPointsLocally();

  elevateBezierCurves(m_mesh);

  for(int d = 2; d <= md; ++d){
    if(!fs->hasNodesIn(d) ||
        getBlendingOrder(apf::Mesh::simplexTypes[d])) continue;
    int type = apf::Mesh::simplexTypes[d];
    int n = fs->getEntityShape(type)->countNodes();
    int ne = fs->countNodesOn(type);
    apf::NewArray<double> c;
    getGregoryBlendedTransformationCoefficients(1,type,c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      if(!isBoundaryEntity(m_mesh,e) && m_mesh->isOwned(e))
        convertInterpolationPoints(m_mesh,e,n-ne,ne,c);
    }
    m_mesh->end(it);
  }

  synchronize();

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

} // namespace crv
