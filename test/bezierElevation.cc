#include <crv.h>
#include <crvBezier.h>
#include <crvSnap.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <mth.h>
#include <mth_def.h>
#include <cassert>

/*
 * This contains all the tests for bezier elevation
 */

void vertFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}
// simple quartic edge
void edgeFunction(double const p[2], double x[3], void*)
{
  x[0] = p[0]*p[0]*p[0]*p[0];
  x[1] = p[0]*p[0];
  x[2] = p[0];
}
void reparam_zero(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 0;
  to[1] = 0;
}
void reparam_one(double const from[2], double to[2], void*)
{
  (void)from;
  to[0] = 1;
  to[1] = 0;
}

agm_bdry add_bdry(gmi_model* m, gmi_ent* e)
{
  return agm_add_bdry(gmi_analytic_topo(m), agm_from_gmi(e));
}

agm_use add_adj(gmi_model* m, agm_bdry b, int tag)
{
  agm* topo = gmi_analytic_topo(m);
  int dim = agm_dim_from_type(agm_bounds(topo, b).type);
  gmi_ent* de = gmi_find(m, dim - 1, tag);
  return agm_add_use(topo, b, agm_from_gmi(de));
}

gmi_model* makeEdgeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edPer = 0;
  double edRan[2] = {0, 1};
  gmi_add_analytic(model, 0, 0, vertFunction, NULL,NULL,NULL);
  gmi_add_analytic(model, 0, 1, vertFunction, NULL,NULL,NULL);
  gmi_ent* ed = gmi_add_analytic(model, 1, 0, edgeFunction, &edPer, &edRan, 0);
  agm_bdry b = add_bdry(model, ed);
  agm_use u0 = add_adj(model, b, 0);
  gmi_add_analytic_reparam(model, u0, reparam_zero, 0);
  agm_use u1 = add_adj(model, b, 1);
  gmi_add_analytic_reparam(model, u1, reparam_one, 0);
  return model;
}

// edges go counter clockwise

// face areas are 1/2 and 19/30
void vert0(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}
// edges go counter clockwise
void edge0(double const p[2], double x[3], void*)
{
  x[0] = p[0];
  x[1] = 0.;
}
void edge1(double const p[2], double x[3], void*)
{
  x[0] = 1.0-p[0]*(p[0]-1.0)*p[0]*(p[0]-1.0);
  x[1] = p[0];
}
void edge2(double const p[2], double x[3], void*)
{
  double u = 1.-p[0];
  x[0] = u;
  x[1] = u;
}

void face0(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}

void make_edge_topo(gmi_model* m, gmi_ent* e, int v0tag, int v1tag)
{
  agm_bdry b = add_bdry(m, e);
  agm_use u0 = add_adj(m, b, v0tag);
  gmi_add_analytic_reparam(m, u0, reparam_zero, 0);
  agm_use u1 = add_adj(m, b, v1tag);
  gmi_add_analytic_reparam(m, u1, reparam_one, 0);
}

gmi_model* makeFaceModel()
{
  gmi_model* model = gmi_make_analytic();
  int edPer = 0;
  double edRan[2] = {0, 1};
  for(int i = 0; i < 3; ++i)
    gmi_add_analytic(model, 0, i, vertFunction,NULL,NULL,NULL);
  gmi_ent* eds[3];
  eds[0] = gmi_add_analytic(model, 1, 0, edge0, &edPer, &edRan, 0);
  eds[1] = gmi_add_analytic(model, 1, 1, edge1, &edPer, &edRan, 0);
  eds[2] = gmi_add_analytic(model, 1, 2, edge2, &edPer, &edRan, 0);
  for(int i = 0; i < 3; ++i)
    make_edge_topo(model, eds[i], i, (i+1) % 3);
  int faPer[2] = {0, 0};
  double faRan[2][2] = {{0,1},{0,1}};
  gmi_add_analytic(model, 2, 0, face0, faPer, faRan, 0);

  return model;
}

apf::Mesh2* createMesh2D()
{
  gmi_model* model = makeFaceModel();
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, true);
  apf::MeshEntity* v[3], *edges[3];
  apf::Vector3 points2D[3] =
  {apf::Vector3(0,0,0),apf::Vector3(1,0,0),apf::Vector3(1,1,0)};

  for (int i = 0; i < 3; ++i){
    v[i] = m->createVertex(m->findModelEntity(0,i),points2D[i],points2D[i]);
  }
  for (int i = 0; i < 3; ++i){
    apf::ModelEntity* edge = m->findModelEntity(1,i);
    apf::MeshEntity* ved[2] = {v[i],v[(i+1) % 3]};
    edges[i] = m->createEntity(apf::Mesh::EDGE,edge,ved);
  }
  apf::ModelEntity* faceModel = m->findModelEntity(2,0);
  m->createEntity(apf::Mesh::TRIANGLE,faceModel,edges);

  m->acceptChanges();
  m->verify();
  return m;
}

static apf::Vector3 points3D[4] =
{apf::Vector3(0,0,0),
    apf::Vector3(1,0,0),
    apf::Vector3(0,1,0),
    apf::Vector3(0,0,1)};


apf::Mesh2* createMesh3D()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, true);

  apf::buildOneElement(m,0,apf::Mesh::TET,points3D);
  apf::deriveMdsModel(m);

  m->acceptChanges();
  m->verify();
  return m;
}

void testEdgeElevation()
{
  for (int o = 1; o <= 5; ++o){

    gmi_model* model = makeEdgeModel();
    apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 1, true);

    apf::Vector3 points[2] = {apf::Vector3(0,0,0),apf::Vector3(1,1,1)};
    apf::MeshEntity* v[2];
    for (int i = 0; i < 2; ++i)
      v[i] = m->createVertex(m->findModelEntity(0,i),points[i],points[i]);

    apf::MeshEntity* edge = m->createEntity(apf::Mesh::EDGE,
        m->findModelEntity(1,0), v);

    m->acceptChanges();
    // curve the mesh
    crv::BezierCurver bc(m,o,0,3);
    bc.run();

    apf::NewArray<apf::Vector3> pts(101);
    apf::NewArray<apf::Vector3> nodes;

    apf::Element* elem = apf::createElement(m->getCoordinateField(),edge);
    apf::getVectorNodes(elem,nodes);
    // need to precompute values, since the edge is elevated and the old one
    // no longer exists
    apf::Vector3 p, pt;
    for (int i = 0; i <= 100; ++i){
      p[0] = 0.02*i-1.;
      apf::getVector(elem,p,pts[i]);
    }
    apf::destroyElement(elem);

    // elevate everything to 10th order
    apf::NewArray<apf::Vector3> elevatedNodes;
    elevatedNodes.allocate(7);
    crv::elevateBezierEdge(o,6-o,nodes,elevatedNodes);
    apf::changeMeshShape(m,crv::getBezier(6),false);

    m->setPoint(v[0],0,points[0]);
    m->setPoint(v[1],0,points[1]);
    for (int j = 1; j <= 5; ++j)
      m->setPoint(edge,j-1,elevatedNodes[j]);

    elem = apf::createElement(m->getCoordinateField(),edge);
    for (int i = 0; i <= 100; ++i){
      p[0] = 0.02*i-1.;
      apf::getVector(elem,p,pt);
      assert(std::abs((pts[i]-pt).getLength()) < 1e-15);
    }
    apf::destroyElement(elem);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

void testTriElevation()
{

  for (int o = 1; o <= 5; ++o){

    apf::Mesh2* m = createMesh2D();
    crv::BezierCurver bc(m,o,0);
    bc.run();
    apf::MeshIterator* it = m->begin(2);
    apf::MeshEntity* tri = m->iterate(it);
    m->end(it);

    apf::Element* elem = apf::createElement(m->getCoordinateField(),tri);

    apf::MeshEntity* verts[3],* edges[3];
    m->getDownward(tri,0,verts);
    m->getDownward(tri,1,edges);

    apf::NewArray<apf::Vector3> pts(66);
    apf::NewArray<apf::Vector3> nodes;

    apf::getVectorNodes(elem,nodes);

    // again, precompute values, since the edge is elevated and the old one
    // no longer exists
    apf::Vector3 p, pt;
    int n = 0;
    for (int j = 0; j <= 10; ++j){
      p[1] = 1.*j/10;
      for (int i = 0; i <= 10-j; ++i){
        p[0] = 1.*i/10;
        apf::getVector(elem,p,pts[n]);
        n++;
      }
    }
    apf::destroyElement(elem);

    // elevate everything to 6th order
    apf::NewArray<apf::Vector3> elevatedNodes;
    elevatedNodes.allocate(28);
    crv::elevateBezierTriangle(o,6-o,nodes,elevatedNodes);
    apf::changeMeshShape(m,crv::getBezier(6),false);

    m->setPoint(verts[0],0,elevatedNodes[0]);
    m->setPoint(verts[1],0,elevatedNodes[1]);
    m->setPoint(verts[2],0,elevatedNodes[2]);

    for (int i = 0; i < 3; ++i)
      for (int j = 1; j <= 5; ++j)
        m->setPoint(edges[i],j-1,elevatedNodes[2+5*i+j]);
    for (int j = 0; j < 10; ++j)
      m->setPoint(tri,j,elevatedNodes[18+j]);

    elem = apf::createElement(m->getCoordinateField(),tri);
    n = 0;
    for (int j = 0; j <= 10; ++j){
      p[1] = 1.*j/10;
      for (int i = 0; i <= 10-j; ++i){
        p[0] = 1.*i/10;
        apf::getVector(elem,p,pt);
        assert(std::abs((pts[n]-pt).getLength()) < 1e-15);
        n++;
      }
    }
    apf::destroyElement(elem);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

void testTetElevation()
{
  gmi_register_null();

  for (int order = 1; order <= 4; ++order){

    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, crv::getBezier(order),true);
    apf::FieldShape* fs = m->getShape();
    crv::BezierCurver bc(m,order,0);
    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = crv::getNumControlPoints(d,order);
      int ni = fs->countNodesOn(d);
      if(ni <= 0) continue;
      apf::NewArray<double> c;
      crv::getBezierTransformationCoefficients(m,order,d,c);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        if(m->getModelType(m->toModel(e)) == m->getDimension()) continue;
        bc.convertInterpolationPoints(e,n,ni,c);
      }
      m->end(it);
    }
    apf::MeshIterator* it = m->begin(3);
    apf::MeshEntity* tet = m->iterate(it);
    m->end(it);

    apf::Element* elem = apf::createElement(m->getCoordinateField(),tet);

    apf::MeshEntity* verts[4],* edges[6],* faces[4];
    m->getDownward(tet,0,verts);
    m->getDownward(tet,1,edges);
    m->getDownward(tet,2,faces);

    apf::NewArray<apf::Vector3> pts(35);
    apf::NewArray<apf::Vector3> nodes;

    apf::getVectorNodes(elem,nodes);
    apf::Vector3 p, pt;
    int n = 0;
    for (int k = 0; k <= 4; ++k){
      p[2] = 0.25*k;
      for (int j = 0; j <= 4-k; ++j){
        p[1] = 0.25*j;
        for (int i = 0; i <= 4-j-k; ++i){
          p[0] = 0.25*i;
          apf::getVector(elem,p,pts[n]);
          n++;
        }
      }
    }
    apf::destroyElement(elem);

    // elevate everything to new order
    int elevatedOrder = 18;

    apf::NewArray<apf::Vector3> elevatedNodes;
    elevatedNodes.allocate(crv::getNumControlPoints(apf::Mesh::TET,elevatedOrder));
    crv::elevateBezierTet(order,elevatedOrder-order,nodes,elevatedNodes);
    apf::changeMeshShape(m,crv::getBezier(elevatedOrder),false);
    for (int i = 0; i < 4; ++i)
      m->setPoint(verts[i],0,elevatedNodes[i]);

    int nE = crv::getNumInternalControlPoints(apf::Mesh::EDGE,elevatedOrder);
    int nF = crv::getNumInternalControlPoints(apf::Mesh::TRIANGLE,elevatedOrder);
    int nT = crv::getNumInternalControlPoints(apf::Mesh::TET,elevatedOrder);

    int which, rotate;
    bool flip;
    for (int i = 0; i < 6; ++i){
      apf::getAlignment(m,tet,edges[i],which,flip,rotate);
      if(!flip)
        for (int j = 0; j < nE; ++j)
          m->setPoint(edges[i],j,elevatedNodes[4+nE*i+j]);
      else
        for (int j = 0; j < nE; ++j)
          m->setPoint(edges[i],nE-1-j,elevatedNodes[4+nE*i+j]);
    }
    for (int i = 0; i < 4; ++i)
      for (int j = 0; j < nF; ++j)
        m->setPoint(faces[i],j,elevatedNodes[4+6*nE+nF*i+j]);
    for (int i = 0; i < nT; ++i)
      m->setPoint(tet,i,elevatedNodes[4+6*nE+4*nF+i]);

    elem = apf::createElement(m->getCoordinateField(),tet);
    n = 0;
    for (int k = 0; k <= 4; ++k){
      p[2] = 0.25*k;
      for (int j = 0; j <= 4-k; ++j){
        p[1] = 0.25*j;
        for (int i = 0; i <= 4-j-k; ++i){
          p[0] = 0.25*i;
          apf::getVector(elem,p,pt);
          assert(std::abs((pts[n]-pt).getLength()) < 1e-15);
          n++;
        }
      }
    }

    apf::destroyElement(elem);
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  testEdgeElevation();
  testTriElevation();
  testTetElevation();
  PCU_Comm_Free();
  MPI_Finalize();
}
