/******************************************************************************

  Copyright 2013 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

 *******************************************************************************/
#include "maCurveMesh.h"
#include "maSnap.h"
#include "maAdapt.h"
#include <apfField.h>
#include <apfShape.h>
#include <apfMesh.h>
#include <gmi.h>

#include <fstream>
#include <sstream>

namespace ma {

MeshCurver::MeshCurver(Adapt *a, int o)
{
  adapt = a;
  order = o;
}

void MeshCurver::snapToInterpolateEdge(Entity* e)
{
  Mesh* m = adapt->mesh;
  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  int non = fs->countNodesOn(Mesh::EDGE);
  Vector p, xi, pt;
  for(int i = 0; i < non; ++i){
    Model* g = m->toModel(e);
    fs->getNodeXi(Mesh::EDGE,i,xi);
    transferParametricOnEdgeSplit(m,e,0.5*(xi[0]+1.),p);

    m->snapToModel(g,p,pt);
    m->setPoint(e,i,pt);
  }
}

void MeshCurver::snapToInterpolateTri(Entity* e)
{
  Mesh* m = adapt->mesh;
  apf::FieldShape * fs = m->getCoordinateField()->getShape();
  int non = fs->countNodesOn(Mesh::TRIANGLE);
  Vector p, xi, pt;
  for(int i = 0; i < non; ++i){
    Model* g = m->toModel(e);
    fs->getNodeXi(Mesh::TRIANGLE,i,xi);
    transferParametricOnTriSplit(m,e,xi,p);
    m->snapToModel(g,p,pt);
    m->setPoint(e,i,pt);
  }
}

void MeshCurver::snapToInterpolate(int dim)
{
  Mesh* m = adapt->mesh;
  int t = (dim == 1) ? Mesh::EDGE : Mesh::TRIANGLE;
  Entity* e;
  Vector p, xi, pt;
  Iterator* it = m->begin(dim);
  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) == m->getDimension()) continue;
    if(t == Mesh::EDGE)
      snapToInterpolateEdge(e);
    else
      snapToInterpolateTri(e);
  }
  m->end(it);
}

void MeshCurver::convertInterpolationPoints(Entity* e,
    int n, int ne, apf::NewArray<double>& c){

  Mesh* m = adapt->mesh;
  apf::NewArray<Vector> l, b(ne);
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  apf::getVectorNodes(elem,l);

  for(int i = 0; i < ne; ++i)
    b[i].zero();

  for( int i = 0; i < ne; ++i)
    for( int j = 0; j < n; ++j)
      b[i] += l[j]*c[i*n+j];

  for(int i = 0; i < ne; ++i)
    m->setPoint(e,i,b[i]);

  apf::destroyElement(elem);
}

bool BezierCurver::run()
{
  if(order < 1 || order > 6){
    fprintf(stderr,"Warning: cannot convert to Bezier of order %d\n",order);
    return false;
  }
  Mesh* m = adapt->mesh;
  int md = m->getDimension();
  apf::changeMeshShape(m, apf::getBezier(md,order),true);
  apf::FieldShape * fs = m->getCoordinateField()->getShape();

  // interpolate points in each dimension
  for(int d = 1; d < md; ++d)
    snapToInterpolate(d);

  // go downward, and convert interpolating to control points
  for(int d = md-1; d >= 1; --d){
    int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
    int ne = fs->countNodesOn(d);

    apf::NewArray<double> c;
    apf::getTransformationCoefficients(order,md,d,c);
    Entity* e;
    Iterator* it = m->begin(d);
    while ((e = m->iterate(it))) {
      Model* g = m->toModel(e);
      if(m->getModelType(g) == m->getDimension()) continue;
      convertInterpolationPoints(e,n,ne,c);
    }
    m->end(it);
  }
  m->acceptChanges();
  m->verify();

  return true;
}

/* Elevates a bezier curve from order n to order n+r */
static void elevateBezierCurve(Mesh* m, Entity* edge, int n, int r)
{

  apf::Element* elem =
       apf::createElement(m->getCoordinateField(),edge);

  Vector pt;
  apf::NewArray<Vector> p;
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
      pt += p[map[j]]*apf::binomial(n,j)*
          apf::binomial(r,i-j)/
          apf::binomial(n+r,i);
    m->setPoint(edge,i-1,pt);
  }

  apf::destroyElement(elem);
}

void GregoryCurver::setCubicEdgePointsUsingNormals()
{
  Mesh* m = adapt->mesh;
  Entity* e;
  Vector p, xi, pt;
  Vector points[4];

  Iterator* it = m->begin(1);

  apf::NewArray<double> c;
  apf::getTransformationCoefficients(4,3,1,c);

  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) == m->getDimension()) continue;
    if(m->getModelType(g) == 1) {
      Vector t[2];
      Entity* v[2];
      m->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m->getPoint(v[i],0,points[i*3]);
        m->getParamOn(g,v[i],p);
        m->getFirstDerivative(g,p,t[i],t[i]);
        t[i] = t[i].normalize();
      }
        double d = (points[3]-points[0]).getLength();
        Vector l = (points[3]-points[0])/d;
        points[1] = points[0] + t[0]*(l*t[0])/fabs(l*t[0])*d/3.;
        points[2] = points[3] - t[1]*(l*t[1])/fabs(l*t[1])*d/3.;
        for(int i = 0; i < 2; ++i)
          m->setPoint(e,i,points[i+1]);
    } else {
      Vector n[2];
      Entity* v[2];
      m->getDownward(e,0,v);
      for(int i = 0; i < 2; ++i){
        m->getPoint(v[i],0,points[i*3]);
        m->getParamOn(g,v[i],p);
        m->getNormal(g,p,n[i]);
      }
      double d = (points[3]-points[0]).getLength();
      Vector l = (points[3]-points[0])/d;
      double a[3] = {n[0]*l,n[1]*l,n[0]*n[1]};

      double r = 6.*(2.*a[0]+a[2]*a[1])/(4.-a[2]*a[2]);
      double s = 6.*(2.*a[1]+a[2]*a[0])/(4.-a[2]*a[2]);
      points[1] = points[0] + (l*6. - n[0]*2.*r+ n[1]*s)   *d/18.;
      points[2] = points[3] - (l*6. + n[0]*r   - n[1]*2.*s)*d/18.;
      for(int i = 0; i < 2; ++i)
        m->setPoint(e,i,points[i+1]);
    }
  }
  m->end(it);
}

static void elevateBezierCurves(Mesh* m)
{

  Entity* e;
  Iterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) == 3) continue;
    elevateBezierCurve(m,e,3,1);
  }
  m->end(it);
}

void GregoryCurver::setInternalPointsUsingNeighbors()
{
  Mesh* m = adapt->mesh;

  Entity* e;
  Iterator* it = m->begin(1);
  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) != 2) continue;
    int tag = m->getModelTag(g);
    apf::Up up;
    m->getUp(e,up);
    Entity* faces[2];
    int iF = 0;
    for(int i = 0; i < up.n; ++i){
      if(m->getModelTag(m->toModel(up.e[i])) == tag)
        faces[iF++] = up.e[i];

    }
    assert(m->getModelType(m->toModel(faces[0])) ==
           m->getModelType(m->toModel(faces[1])) );
    // now we have the faces
    int which[2], rotate[2];
    bool flip[2];
    apf::getAlignment(m,faces[0],e,which[0],flip[0],rotate[0]);
    apf::getAlignment(m,faces[1],e,which[1],flip[1],rotate[1]);
  }
  m->end(it);
}

void GregoryCurver::setInternalPointsLocally()
{
  Mesh* m = adapt->mesh;

  Vector D[3][4];
  Vector W[3][3];
  Vector A[3][3];

  double lam[3][2];
  double mu[3][2];
  Vector G[6];

  Entity* e;
  Iterator* it = m->begin(2);
  while ((e = m->iterate(it))) {
    Model* g = m->toModel(e);
    if(m->getModelType(g) != 2) continue;

    Vector n[3];
    Entity* verts[3];
    Entity* edges[3];
    m->getDownward(e,0,verts);
    m->getDownward(e,1,edges);

    // elevated edges
    apf::NewArray<Vector> q(12);

    for(int i = 0; i < 3; ++i){
      Vector param;
      m->getPoint(verts[i],0,q[i]);
      m->getParamOn(g,verts[i],param);
      m->getNormal(g,param,n[i]);
    }

    // elevate the edge points without formally setting them to q
    // compute tangent vectors, W
    for(int i = 0; i < 3; ++i){
      apf::Element* edge =
          apf::createElement(m->getCoordinateField(),edges[i]);
      apf::NewArray<Vector> ep;
      apf::getVectorNodes(edge,ep);

      bool flip;
      int which, rotate;
      apf::getAlignment(m,e,edges[i],which,flip,rotate);

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
      m->setPoint(e,i,G[i]);

  }
  m->end(it);
}

bool GregoryCurver::run()
{
  Mesh* m = adapt->mesh;
  if(order != 4){
    fprintf(stderr,"Warning: cannot convert to Gregory of order %d\n",order);
    return false;
  }
  if(m->getDimension() != 3){
    fprintf(stderr,"Warning: can only convert 3D Mesh to Gregory\n");
    return false;
  }

  apf::changeMeshShape(m, apf::getGregory(order),true);

  setCubicEdgePointsUsingNormals();
  setInternalPointsLocally();
  elevateBezierCurves(m);

  return true;
}

double interpolationError(Mesh* m, Entity* e, int n){
  Model* g = m->toModel(e);
  if (m->getModelType(g) == m->getDimension())
	  return 0.;
  int d = apf::getDimension(m,e);
  int nj = (d == 2) ? n : 1;
  Vector pt,pa(0.,0.,0.),cpt,cpa;
  double max = 0.0;
  apf::Element* elem =
      apf::createElement(m->getCoordinateField(),e);
  for (int j = 0; j < nj; ++j){
    pa[1] = 1.*j/(nj-1.);
    for (int i = 0; i < n-j; ++i){
      if(d == 1)
    	  pa[0] = 2.*i/(n-1)-1.0;
      else
    	  pa[0] = 1.*i/(nj-1.);
      apf::getVector(elem,pa,pt);
      m->getClosestPoint(g,pt,cpt,cpa);
      max = std::max((cpt-pt).getLength(),max);
    }
  }
  apf::destroyElement(elem);
  return max;
}

void writePointSet(Mesh* m, int d, int n, const char* prefix)
{
  int nj = (d > 1) ? n : 1;
  int nk = (d == 3) ? n : 1;

  apf::DynamicArray<apf::Vector3> pts(0);

  Iterator* it = m->begin(d);
  Entity* e;
  Vector pa,pt;

  while ((e = m->iterate(it))) {
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);
    for (int k = 0; k < nk; ++k){
      pa[2] = 1.*k/nk;
      for (int j = 0; j < nj; ++j){
        pa[1] = 1.*j/nj;
        for (int i = 0; i < n+1-j-k; ++i){
          if(d == 1) pa[0] = 2.*i/n-1.;
          else pa[0] = 1.*i/nj;
          apf::getVector(elem,pa,pt);
          pts.append(pt);
        }
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);

  std::stringstream ss;
  ss << prefix << "_" << d << "_"
      << m->getCoordinateField()->getShape()->getName() << ".csv";

  std::string fileName = ss.str();
  std::ofstream file(fileName.c_str());
  assert(file.is_open());
  int size = pts.getSize();
  file << "x,y,z\n";
  for(int p = 0; p < size; ++p)
    file << pts[p][0] << ","<< pts[p][1]<< "," << pts[p][2] << "\n";
}

} // namespace ma

