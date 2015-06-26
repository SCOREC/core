/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/
#include "crv.h"
#include "crvSnap.h"

#include <maSnap.h>
#include <apfField.h>

namespace crv {

void MeshCurver::snapToInterpolateEdge(apf::MeshEntity* e)
{
  apf::FieldShape * fs = m_mesh->getCoordinateField()->getShape();
  int non = fs->countNodesOn(apf::Mesh::EDGE);
  apf::Vector3 p, xi, pt;
  for(int i = 0; i < non; ++i){
    apf::ModelEntity* g = m_mesh->toModel(e);
    fs->getNodeXi(apf::Mesh::EDGE,i,xi);
    ma::transferParametricOnEdgeSplit(m_mesh,e,0.5*(xi[0]+1.),p);
    m_mesh->snapToModel(g,p,pt);
    m_mesh->setPoint(e,i,pt);
  }
}

void MeshCurver::snapToInterpolateTri(apf::MeshEntity* e)
{
  apf::FieldShape * fs = m_mesh->getCoordinateField()->getShape();
  int non = fs->countNodesOn(apf::Mesh::TRIANGLE);
  apf::Vector3 p, xi, pt;
  for(int i = 0; i < non; ++i){
    apf::ModelEntity* g = m_mesh->toModel(e);
    fs->getNodeXi(apf::Mesh::TRIANGLE,i,xi);
    ma::transferParametricOnTriSplit(m_mesh,e,xi,p);
    m_mesh->snapToModel(g,p,pt);
    m_mesh->setPoint(e,i,pt);
  }
}

void MeshCurver::snapToInterpolate(int dim)
{
  int t = (dim == 1) ? apf::Mesh::EDGE : apf::Mesh::TRIANGLE;
  apf::MeshEntity* e;
  apf::Vector3 p, xi, pt;
  apf::MeshIterator* it = m_mesh->begin(dim);
  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);
    if(m_mesh->getModelType(g) == m_mesh->getDimension()) continue;
    if(t == apf::Mesh::EDGE)
      snapToInterpolateEdge(e);
    else
      snapToInterpolateTri(e);
  }
  m_mesh->end(it);
}

void MeshCurver::convertInterpolationPoints(apf::MeshEntity* e,
    int n, int ne, apf::NewArray<double>& c){

  apf::NewArray<apf::Vector3> l, b(ne);
  apf::Element* elem =
      apf::createElement(m_mesh->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  for(int i = 0; i < ne; ++i)
    b[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      b[i] += l[j]*c[i*n+j];

  for(int i = 0; i < ne; ++i)
    m_mesh->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}

bool InterpolatingCurver::run()
{
  // interpolate points in each dimension
  for(int d = 1; d < m_mesh->getDimension(); ++d)
    snapToInterpolate(d);

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

bool BezierCurver::run()
{
  if(m_order < 1 || m_order > 6){
    fail("trying to convert to unimplemented Bezier order\n");
  }

  int md = m_mesh->getDimension();
  apf::changeMeshShape(m_mesh, getBezier(md,m_order,m_blendOrder),true);
  apf::FieldShape * fs = m_mesh->getCoordinateField()->getShape();

  // interpolate points in each dimension
  for(int d = 1; d < md; ++d)
    snapToInterpolate(d);

  // go downward, and convert interpolating to control points
  for(int d = md-1; d >= 1; --d){
    int n = (d == 2)? (m_order+1)*(m_order+2)/2 : m_order+1;
    int ne = fs->countNodesOn(d);

    apf::NewArray<double> c;
    getTransformationCoefficients(md,d,c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);
    while ((e = m_mesh->iterate(it))){
      convertInterpolationPoints(e,n,ne,c);
    }
    m_mesh->end(it);
  }
  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

/* Elevates a bezier curve from order n to order n+r */
static void elevateBezierCurve(apf::Mesh2* m, apf::MeshEntity* edge, int n, int r)
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
    apf::ModelEntity* g = m->toModel(e);
    if(m->getModelType(g) == 3) continue;
    elevateBezierCurve(m,e,3,1);
  }
  m->end(it);
}

void GregoryCurver::setInternalPointsUsingNeighbors()
{
  apf::MeshEntity* e;
  apf::MeshIterator* it = m_mesh->begin(1);
  while ((e = m_mesh->iterate(it))) {
    apf::ModelEntity* g = m_mesh->toModel(e);
    if(m_mesh->getModelType(g) != 2) continue;
    int tag = m_mesh->getModelTag(g);
    apf::Up up;
    m_mesh->getUp(e,up);
    apf::MeshEntity* faces[2];
    int iF = 0;
    for(int i = 0; i < up.n; ++i){
      if(m_mesh->getModelTag(m_mesh->toModel(up.e[i])) == tag)
        faces[iF++] = up.e[i];

    }
    assert(m_mesh->getModelType(m_mesh->toModel(faces[0])) ==
        m_mesh->getModelType(m_mesh->toModel(faces[1])) );
    // now we have the faces
    int which[2], rotate[2];
    bool flip[2];
    apf::getAlignment(m_mesh,faces[0],e,which[0],flip[0],rotate[0]);
    apf::getAlignment(m_mesh,faces[1],e,which[1],flip[1],rotate[1]);
  }
  m_mesh->end(it);
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
    if(m_mesh->getModelType(g) != 2) continue;

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
  if(m_order < 3 || m_order > 4){
    fail("cannot convert to G1 of this order\n");
  }
  if(m_mesh->getDimension() != 3){
    fail("can only convert 3D Mesh to G1 continuous surface\n");
  }

  apf::changeMeshShape(m_mesh, getGregory(m_order,m_blendOrder),true);

  int md = m_mesh->getDimension();

  apf::FieldShape * fs = m_mesh->getCoordinateField()->getShape();

  // interpolate points in each dimension
  for(int d = 1; d < md; ++d)
    snapToInterpolate(d);

  // go downward, and convert interpolating to control points
  for(int d = md-1; d >= 1; --d){

    int n = (d == 2)? (m_order+1)*(m_order+2)/2 : m_order+1;
    int ne = fs->countNodesOn(d);
    apf::NewArray<apf::Vector3> l, b(ne);

    apf::NewArray<double> c;
    getTransformationCoefficients(md,d,c);

    apf::MeshEntity* e;
    apf::MeshIterator* it = m_mesh->begin(d);

    while ((e = m_mesh->iterate(it))) {
      apf::ModelEntity* g = m_mesh->toModel(e);
      // we'll set the points to be G1 continuous later on
      if (m_mesh->getModelType(g) < md) continue;

      // set extra internal points
      if(d == 2 && m_order == 4){

        apf::Element* elem =
            apf::createElement(m_mesh->getCoordinateField(),e);
        apf::getVectorNodes(elem,l);

        for(int i = 0; i < ne; ++i)
          b[i].zero();

        for( int i = 0; i < 3; ++i)
          for( int j = 0; j < n; ++j)
            b[i] += l[j]*c[i*n+j];

        int map[3] = {1,2,0};
        for( int i = 3; i < 6; ++i)
          for( int j = 0; j < n; ++j)
            b[i] += l[j]*c[map[i-3]*n+j];

        for(int i = 0; i < ne; ++i)
          m_mesh->setPoint(e,i,b[i]);

        apf::destroyElement(elem);
      } else
        convertInterpolationPoints(e,n,ne,c);

    }
    m_mesh->end(it);
  }

  if(m_order == 4){
    setCubicEdgePointsUsingNormals();
    setInternalPointsLocally();
    elevateBezierCurves(m_mesh);
  } else
    setCubicEdgePointsUsingNormals();

  m_mesh->acceptChanges();
  m_mesh->verify();
  return true;
}

} // namespace ma
