#include <crv.h>
#include <crvBezier.h>
#include <crvSnap.h>
#include <crvMath.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>
#include <ma.h>
#include <mth.h>
#include <mth_def.h>
#include <cassert>

/* This file contains tests for edge and triangle bezier elevation/subdivision.
 * For subdivision, the approach is to create an edge/triangle, split it
 * and then compare the new subdivisions against the original, which should be
 * exactly identical
 *
 * For elevation, an edge/triangle is elevated to order 6, and again,
 * a check is performed
 *
 * There is also a test for the node numbering functions, comparing the tables
 * against the automatically generated values
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

/*
 * Create a mesh with a single edge,
 * Create two edges as an even subdivision,
 * and keep all three around, using the original one
 * to compare correctness of the split.
 *
 */
void testEdgeSubdivision()
{
  for (int o = 1; o <= 6; ++o){

    gmi_model* model = makeEdgeModel();
    apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 1, true);

    apf::ModelEntity* edgeModel = m->findModelEntity(1,0);

    apf::Vector3 points[2] = {apf::Vector3(0,0,0),apf::Vector3(1,1,1)};
    apf::MeshEntity* v[2];
    for (int i = 0; i < 2; ++i)
      v[i] = m->createVertex(m->findModelEntity(0,i),points[i],points[i]);

    apf::MeshEntity* edge = m->createEntity(apf::Mesh::EDGE,edgeModel, v);

    m->acceptChanges();
    m->verify();

    // curve the mesh
    crv::BezierCurver bc(m,o,0,3);
    bc.run();

    apf::Element* elem = apf::createElement(m->getCoordinateField(),edge);
    apf::NewArray<apf::Vector3> nodes;
    apf::NewArray<apf::Vector3> subNodes[2];
    subNodes[0].allocate(o+1);
    subNodes[1].allocate(o+1);
    apf::getVectorNodes(elem,nodes);
    // subdivide the edge's nodes
    crv::subdivideBezierEdge(o,1./4,nodes,subNodes);

    // create the two new edges
    apf::Vector3 p;
    crv::transferParametricOnEdgeSplit(m,edge,1./4,p);

    apf::MeshEntity* v2 = m->createVertex(edgeModel,subNodes[0][o],p);
    apf::MeshEntity* vE[2][2] = {{v[0],v2},{v2,v[1]}};
    apf::MeshEntity* e[2];

    for (int i = 0; i < 2; ++i){
      e[i] = m->createEntity(apf::Mesh::EDGE,edgeModel,vE[i]);
      for (int j = 1; j < o; ++j)
        m->setPoint(e[i],j-1,subNodes[i][j]);
    }


    // compare the two curves to the original one
    apf::Element* elem0 = apf::createElement(m->getCoordinateField(),e[0]);
    apf::Element* elem1 = apf::createElement(m->getCoordinateField(),e[1]);

    apf::Vector3 pt1, pt2, p1, p2;
    for (int i = 0; i <= 100; ++i){
      p1[0] = 0.02*i-1.;
      p2[0] = 0.005*i-1.;
      apf::getVector(elem0,p1,pt1);
      apf::getVector(elem,p2,pt2);
      assert(std::abs((pt2-pt1).getLength()) < 1e-15);

      p2[0] = 0.015*i-0.5;
      apf::getVector(elem1,p1,pt1);
      apf::getVector(elem,p2,pt2);
      assert(std::abs((pt2-pt1).getLength()) < 1e-15);
    }
    apf::destroyElement(elem);
    apf::destroyElement(elem0);
    apf::destroyElement(elem1);

    m->destroyNative();
    apf::destroyMesh(m);
  }
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

/* Create a single triangle, split it, and make sure each triangle
 * exactly replicates its part of the old one
 *
 */
void testTriSubdivision1()
{

  for(int o = 1; o <= 6; ++o){
    apf::Mesh2* m = createMesh2D();
    crv::BezierCurver bc(m,o,0);
    bc.run();

    apf::MeshIterator* it = m->begin(2);
    apf::MeshEntity* e = m->iterate(it);
    m->end(it);

    apf::Element* elem = apf::createElement(m->getCoordinateField(),e);

    apf::MeshEntity* verts[3],* edges[3];
    m->getDownward(e,0,verts);
    m->getDownward(e,1,edges);

    apf::NewArray<apf::Vector3> nodes;
    apf::NewArray<apf::Vector3> subNodes[3];
    for (int t = 0; t < 3; ++t)
      subNodes[t].allocate((o+1)*(o+2)/2);

    apf::getVectorNodes(elem,nodes);
    apf::Vector3 splitpt(0,0.5,0.5);

    crv::subdivideBezierTriangle(o,splitpt,nodes,subNodes);

    apf::MeshEntity* v3 = m->createVertex(m->findModelEntity(2,0),
        subNodes[0][2],splitpt);
    apf::MeshEntity* newFaces[3],* newEdges[3];
    for (int i = 0; i < 3; ++i){
      apf::MeshEntity* vE[2] = {v3,verts[i]};
      newEdges[i] = m->createEntity(apf::Mesh::EDGE,m->findModelEntity(1,i),
          vE);
      for (int j = 0; j < o-1; ++j){
        m->setPoint(newEdges[i],j,subNodes[i][3+2*(o-1)+j]);
      }
    }
    for (int i = 0; i < 3; ++i){
      apf::MeshEntity* eF[3] = {edges[i],newEdges[(i+1) % 3],
          newEdges[i]};
      newFaces[i] = m->createEntity(apf::Mesh::TRIANGLE,m->findModelEntity(2,0),
          eF);
      for (int j = 0; j < (o-1)*(o-2)/2; ++j){
        m->setPoint(newFaces[i],j,subNodes[i][3+3*(o-1)+j]);
      }
    }

    // compare the three faces to the original one
    apf::Element* elems[3] =
      {apf::createElement(m->getCoordinateField(),newFaces[0]),
       apf::createElement(m->getCoordinateField(),newFaces[1]),
       apf::createElement(m->getCoordinateField(),newFaces[2])};

    apf::Vector3 p,pOld,pt,ptOld;
    for (int j = 0; j <= 10; ++j){
      p[1] = 1.*j/10;
      for (int i = 0; i <= 10-j; ++i){
        p[0] = 1.*i/10;
        p[2] = 1.-p[0]-p[1];
        // p[1] is the new split point, rescale from new to old
        for (int t = 0; t < 3; ++t){
          apf::getVector(elems[t],p,pt);
          for (int pi = 0; pi < 3; ++pi)
            pOld[pi] = splitpt[(pi+2) % 3]*p[1];
          pOld[t] += p[0];
          pOld[(t+2) % 3] += p[2];
          apf::getVector(elem,pOld,ptOld);
          assert(std::abs((ptOld-pt).getLength()) < 1e-15);
        }
      }
    }

    apf::destroyElement(elem);
    m->destroy(e);

    for(int t = 0; t < 3; ++t)
      apf::destroyElement(elems[t]);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

/* Create a single triangle, split it into 4, and try not to crash
 *
 */
void testTriSubdivision4()
{
  for(int o = 2; o <= 2; ++o){
    apf::Mesh2* m = createMesh2D();
    crv::BezierCurver bc(m,o,0);
    bc.run();
    apf::MeshIterator* it = m->begin(2);
    apf::MeshEntity* e = m->iterate(it);
    m->end(it);

    apf::Element* elem = apf::createElement(m->getCoordinateField(),e);

    apf::NewArray<apf::Vector3> nodes;
    apf::NewArray<apf::Vector3> subNodes[4];
    for (int t = 0; t < 4; ++t)
      subNodes[t].allocate((o+1)*(o+2)/2);

    apf::getVectorNodes(elem,nodes);
    apf::Vector3 splitpt(0,0.5,0.5);

    crv::subdivideBezierTriangle(o,nodes,subNodes);

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

apf::Mesh2* createMesh3D()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, true);

  apf::Vector3 points3D[4] =
  {apf::Vector3(0,0,0),
      apf::Vector3(1,0,0),
      apf::Vector3(0,1,0),
      apf::Vector3(0,0,1)};

  apf::buildOneElement(m,0,apf::Mesh::TET,points3D);

  apf::deriveMdsModel(m);

  m->acceptChanges();
  m->verify();
  return m;
}

void testTetElevation()
{
  gmi_register_null();

  for (int order = 1; order <= 4; ++order){

    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, crv::getBezier(order),true);
    crv::setBlendingOrder(0);
    apf::FieldShape* fs = m->getShape();
    crv::BezierCurver bc(m,order,0);
    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = crv::getNumControlPoints(d,order);
      int ne = fs->countNodesOn(d);
      apf::NewArray<double> c;
      crv::getTransformationCoefficients(order,d,c);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        if(m->getModelType(m->toModel(e)) == m->getDimension()) continue;
        bc.convertInterpolationPoints(e,n,ne,c);
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

void testNodeIndexing(){
  for(int P = 1; P <= 10; ++P)
    for(int i = 0; i <= P; ++i)
      for(int j = 0; j <= P-i; ++j)
        assert(crv::computeTriNodeIndex(P,i,j)
          == crv::getTriNodeIndex(P,i,j));

  for(int P = 1; P <= 4; ++P)
    for(int i = 0; i <= P; ++i)
      for(int j = 0; j <= P-i; ++j)
        for(int k = 0; k <= P-i-j; ++k)
          assert(crv::computeTetNodeIndex(P,i,j,k)
            == crv::getTetNodeIndex(P,i,j,k));
}

static double const a_data[35][35] = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.4759624886330922,0.0008234344110729,0,0,0.3882815810042113,0.1187824065033745,0.0161500894482495,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0625,0.0625,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0008234344110728,0.4759624886330922,0,0,0.0161500894482495,0.1187824065033744,0.3882815810042112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0.4759624886330922,0.0008234344110729,0,0,0,-1e-16,0.3882815810042113,0.1187824065033745,0.0161500894482495,0,0,0,0,0,0,0,0,0,0,0,0,0,-1e-16,0,0,0,0,0,0,0,0,0,0,0},
{0,0.0625,0.0625,0,0,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0.0008234344110729,0.4759624886330922,0,0,0,0,0.0161500894482495,0.1187824065033745,0.3882815810042113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0008234344110728,0,0.4759624886330922,0,0,0,0,0,0,0,0.3882815810042112,0.1187824065033744,0.0161500894482495,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0625,0,0.0625,0,0,0,0,0,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.4759624886330922,0,0.0008234344110729,0,0,0,0,0,0,0,0.0161500894482495,0.1187824065033745,0.3882815810042113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.4759624886330922,0,0,0.0008234344110729,0,0,0,0,0,0,0,0,0,0.3882815810042113,0.1187824065033745,0.0161500894482495,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0625,0,0,0.0625,0,0,0,0,0,0,0,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.0008234344110728,0,0,0.4759624886330922,0,0,0,0,0,0,0,0,0,0.0161500894482495,0.1187824065033744,0.3882815810042112,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0.4759624886330922,0,0.0008234344110729,0,0,-1e-16,0,0,0,0,0,0,0,0,0,0.3882815810042113,0.1187824065033745,0.0161500894482495,0,0,0,0,0,0,0,-1e-16,0,0,0,0,0,0,0,0},
{0,0.0625,0,0.0625,0,0,0,0,0,0,0,0,0,0,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0.0008234344110729,0,0.4759624886330922,0,0,0,0,0,0,0,0,0,0,0,0,0.0161500894482495,0.1187824065033745,0.3882815810042113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0.4759624886330922,0.0008234344110729,0,0,0,0,0,0,-1e-16,0,0,0,0,0,0,0,0,0.3882815810042113,0.1187824065033745,0.0161500894482495,0,0,0,0,0,0,0,0,0,0,-1e-16,0,0},
{0,0,0.0625,0.0625,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.25,0.375,0.25,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0.0008234344110729,0.4759624886330922,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0161500894482495,0.1187824065033745,0.3882815810042113,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0.09710323671528701,0.0023806134614917,0.0023806134614917,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0.0095224538459668,0.0142836807689502,0.0095224538459668,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0,0,0,0.1824494429614713,0.0721948696201392,0.0721948696201392,0,0,0,0,0,0,0,0,0,0},
{0.0023806134614917,0.09710323671528701,0.0023806134614917,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0.1536942048287659,0.0912247214807357,0.0240649565400464,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0,0,0,0.0721948696201392,0.1824494429614713,0.0721948696201392,0,0,0,0,0,0,0,0,0,0},
{0.0023806134614917,0.0023806134614917,0.09710323671528701,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0.0240649565400464,0.0912247214807357,0.1536942048287659,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0,0,0,0.0721948696201392,0.0721948696201392,0.1824494429614713,0,0,0,0,0,0,0,0,0,0},
{0.09710323671528701,0.0023806134614917,0,0.0023806134614917,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0.1824494429614713,0.0721948696201392,0.0721948696201392,0,0,0,0,0,0,0},
{0.0023806134614917,0.09710323671528701,0,0.0023806134614917,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0.0721948696201392,0.1824494429614713,0.0721948696201392,0,0,0,0,0,0,0},
{0.0023806134614917,0.0023806134614917,0,0.09710323671528701,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0.0721948696201392,0.0721948696201392,0.1824494429614713,0,0,0,0,0,0,0},
{0,0.09710323671528701,0.0023806134614917,0.0023806134614917,0,0,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0.1824494429614713,0.0721948696201392,0.0721948696201392,0,0,0,0},
{0,0.0023806134614917,0.09710323671528701,0.0023806134614917,0,0,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0.0721948696201392,0.1824494429614713,0.0721948696201392,0,0,0,0},
{0,0.0023806134614917,0.0023806134614917,0.09710323671528701,0,0,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0.0721948696201392,0.0721948696201392,0.1824494429614713,0,0,0,0},
{0.09710323671528701,0,0.0023806134614917,0.0023806134614917,0,0,0,0,0,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0,0,0,0,0,0,0.1824494429614713,0.0721948696201392,0.0721948696201392,0},
{0.0023806134614917,0,0.09710323671528701,0.0023806134614917,0,0,0,0,0,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0.0095224538459668,0.0142836807689502,0.0095224538459668,0,0,0,0.1536942048287659,0.0912247214807357,0.0240649565400464,0,0,0,0,0,0,0,0,0,0.0721948696201392,0.1824494429614713,0.0721948696201392,0},
{0.0023806134614917,0,0.0023806134614917,0.09710323671528701,0,0,0,0,0,0,0.0095224538459668,0.0142836807689502,0.0095224538459668,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0.0240649565400464,0.0912247214807357,0.1536942048287659,0,0,0,0,0,0,0,0,0,0.0721948696201392,0.0721948696201392,0.1824494429614713,0},
{0.00390625,0.00390625,0.00390625,0.00390625,0.015625,0.0234375,0.015625,0.015625,0.0234375,0.015625,0.015625,0.0234375,0.015625,0.015625,0.0234375,0.015625,0.015625,0.0234375,0.015625,0.015625,0.0234375,0.015625,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.046875,0.09375}};

static double const a_invdata[35][35] = {{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1.526804207661546,-0.2500000000000001,0,0,3.375676033412428,-1.28732567375034,0.6884538479994579,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1.017869471774364,1.017869471774365,0,0,-2.709419920941257,4.383100898333787,-2.709419920941258,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-0.2499999999999995,-1.526804207661547,0,0,0.6884538479994572,-1.28732567375034,3.375676033412429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-1.526804207661547,-0.2500000000000002,0,0,0,0,3.375676033412427,-1.287325673750341,0.6884538479994581,0,0,0,0,0,0,0,0,0,0,0,0,0,2.385415368501797e-15,0,0,0,0,0,0,0,0,0,0,0},
    {0,1.017869471774365,1.017869471774364,0,0,0,0,-2.709419920941257,4.383100898333788,-2.709419920941257,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.914606690679597e-15,0,0,0,0,0,0,0,0,0,0,0},
    {0,-0.25,-1.526804207661547,0,0,0,0,0.6884538479994575,-1.287325673750341,3.375676033412428,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-0.2499999999999996,0,-1.526804207661547,0,0,0,0,0,0,0,3.375676033412429,-1.28732567375034,0.6884538479994573,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1.017869471774364,0,1.017869471774365,0,0,0,0,0,0,0,-2.709419920941258,4.383100898333788,-2.709419920941257,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1.526804207661547,0,-0.2500000000000002,0,0,0,0,0,0,0,0.6884538479994581,-1.287325673750341,3.375676033412428,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-1.526804207661546,0,0,-0.2500000000000001,0,0,0,0,0,0,0,0,0,3.375676033412428,-1.28732567375034,0.6884538479994579,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1.017869471774364,0,0,1.017869471774365,0,0,0,0,0,0,0,0,0,-2.709419920941257,4.383100898333787,-2.709419920941258,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {-0.2499999999999995,0,0,-1.526804207661547,0,0,0,0,0,0,0,0,0,0.6884538479994572,-1.28732567375034,3.375676033412429,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,-1.526804207661547,0,-0.2500000000000002,0,0,0,0,0,0,0,0,0,0,0,0,3.375676033412427,-1.287325673750341,0.6884538479994581,0,0,0,0,0,0,0,2.385415368501797e-15,0,0,0,0,0,0,0,0},
    {0,1.017869471774365,0,1.017869471774364,0,0,0,0,0,0,0,0,0,0,0,0,-2.709419920941257,4.383100898333788,-2.709419920941257,0,0,0,0,0,0,0,-1.914606690679597e-15,0,0,0,0,0,0,0,0},
    {0,-0.25,0,-1.526804207661547,0,0,0,0,0,0,0,0,0,0,0,0,0.6884538479994575,-1.287325673750341,3.375676033412428,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {0,0,-1.526804207661547,-0.2500000000000002,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3.375676033412427,-1.287325673750341,0.6884538479994581,0,0,0,0,0,0,0,0,0,0,2.385415368501797e-15,0,0},
    {0,0,1.017869471774365,1.017869471774364,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2.709419920941257,4.383100898333788,-2.709419920941257,0,0,0,0,0,0,0,0,0,0,-1.914606690679597e-15,0,0},
    {0,0,-0.25,-1.526804207661547,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.6884538479994575,-1.287325673750341,3.375676033412428,0,0,0,0,0,0,0,0,0,0,0,0,0},
    {1.473866405971785,-0.4873245458152488,-0.4873245458152483,0,-2.157170326583087,-0.789537182578665,1.002268609355066,0.4569922360188253,0.4160673296503343,0.456992236018825,1.002268609355065,-0.7895371825786652,-2.157170326583087,0,0,0,0,0,0,0,0,0,7.066481928037422,-2.00343662222666,-2.00343662222666,0,0,0,0,0,0,0,0,0,0},
    {-0.4873245458152486,1.473866405971785,-0.4873245458152484,0,1.002268609355066,-0.7895371825786656,-2.157170326583088,-2.157170326583087,-0.7895371825786647,1.002268609355065,0.4569922360188253,0.4160673296503344,0.4569922360188249,0,0,0,0,0,0,0,0,0,-2.00343662222666,7.066481928037422,-2.00343662222666,0,0,0,0,0,0,0,0,0,0},
    {-0.4873245458152489,-0.4873245458152486,1.473866405971784,0,0.4569922360188252,0.4160673296503342,0.4569922360188253,1.002268609355065,-0.7895371825786655,-2.157170326583087,-2.157170326583087,-0.7895371825786657,1.002268609355066,0,0,0,0,0,0,0,0,0,-2.00343662222666,-2.003436622226659,7.066481928037422,0,0,0,0,0,0,0,0,0,0},
    {1.473866405971786,-0.4873245458152485,0,-0.4873245458152488,-2.157170326583088,-0.7895371825786649,1.002268609355066,0,0,0,0,0,0,-2.157170326583087,-0.7895371825786649,1.002268609355066,0.456992236018825,0.4160673296503342,0.4569922360188252,0,0,0,0,0,0,7.066481928037422,-2.00343662222666,-2.00343662222666,0,0,0,0,0,0,0},
    {-0.4873245458152486,1.473866405971785,0,-0.4873245458152487,1.002268609355066,-0.7895371825786656,-2.157170326583088,0,0,0,0,0,0,0.456992236018825,0.4160673296503344,0.4569922360188255,-2.157170326583087,-0.789537182578665,1.002268609355065,0,0,0,0,0,0,-2.00343662222666,7.066481928037422,-2.003436622226661,0,0,0,0,0,0,0},
    {-0.4873245458152488,-0.4873245458152487,0,1.473866405971785,0.4569922360188252,0.4160673296503344,0.4569922360188257,0,0,0,0,0,0,1.002268609355066,-0.7895371825786659,-2.157170326583088,1.002268609355066,-0.7895371825786657,-2.157170326583087,0,0,0,0,0,0,-2.003436622226661,-2.003436622226661,7.066481928037424,0,0,0,0,0,0,0},
    {0,1.473866405971785,-0.4873245458152484,-0.4873245458152483,0,0,0,-2.157170326583087,-0.7895371825786648,1.002268609355065,0,0,0,0,0,0,-2.157170326583087,-0.7895371825786647,1.002268609355065,0.4569922360188252,0.4160673296503343,0.4569922360188252,0,-1.524360512849797e-15,0,0,-1.524360512849797e-15,0,7.066481928037422,-2.00343662222666,-2.00343662222666,0,0,0,0},
    {0,-0.4873245458152486,1.473866405971785,-0.4873245458152485,0,0,0,1.002268609355065,-0.7895371825786658,-2.157170326583087,0,0,0,0,0,0,0.4569922360188253,0.4160673296503344,0.4569922360188254,-2.157170326583087,-0.7895371825786648,1.002268609355065,0,0,0,0,0,0,-2.003436622226661,7.066481928037424,-2.003436622226661,0,-1.524360512849797e-15,0,0},
    {0,-0.4873245458152486,-0.4873245458152487,1.473866405971785,0,0,0,0.4569922360188253,0.4160673296503344,0.4569922360188254,0,0,0,0,0,0,1.002268609355065,-0.7895371825786657,-2.157170326583087,1.002268609355066,-0.7895371825786657,-2.157170326583087,0,0,0,0,0,0,-2.003436622226661,-2.003436622226661,7.066481928037424,0,0,0,0},
    {1.473866405971785,0,-0.4873245458152488,-0.4873245458152486,0,0,0,0,0,0,1.002268609355066,-0.7895371825786652,-2.157170326583087,-2.157170326583087,-0.7895371825786647,1.002268609355066,0,0,0,0.4569922360188252,0.4160673296503343,0.4569922360188252,0,0,0,0,0,0,0,0,0,7.066481928037422,-2.00343662222666,-2.00343662222666,0},
    {-0.4873245458152485,0,1.473866405971785,-0.4873245458152487,0,0,0,0,0,0,-2.157170326583087,-0.7895371825786656,1.002268609355065,0.4569922360188249,0.4160673296503343,0.4569922360188255,0,0,0,-2.157170326583087,-0.7895371825786648,1.002268609355065,0,0,0,0,0,0,0,0,0,-2.003436622226659,7.066481928037422,-2.003436622226661,0},
    {-0.487324545815249,0,-0.4873245458152486,1.473866405971785,0,0,0,0,0,0,0.4569922360188253,0.4160673296503343,0.4569922360188254,1.002268609355066,-0.7895371825786661,-2.157170326583088,0,0,0,1.002268609355066,-0.7895371825786659,-2.157170326583087,0,0,0,0,0,0,0,0,0,-2.003436622226661,-2.00343662222666,7.066481928037424,0},
    {-0.6654926381785982,-0.6654926381785983,-0.6654926381785977,-0.6654926381785985,0.6979094812091967,0.4963403688403295,0.6979094812091966,0.6979094812091965,0.4963403688403295,0.697909481209196,0.6979094812091964,0.4963403688403295,0.6979094812091965,0.6979094812091965,0.4963403688403297,0.6979094812091967,0.6979094812091963,0.4963403688403293,0.6979094812091965,0.697909481209196,0.4963403688403293,0.6979094812091966,-1.529804341792051,-1.529804341792051,-1.529804341792051,-1.529804341792051,-1.52980434179205,-1.529804341792051,-1.529804341792051,-1.529804341792051,-1.529804341792051,-1.529804341792051,-1.52980434179205,-1.529804341792051,10.66666666666667}};

void testMatrixInverse(){
  mth::Matrix<double> A(35,35);
  for (int i = 0; i < 35; ++i)
    for (int j = 0; j < 35; ++j)
      A(i,j) = a_data[i][j];
  mth::Matrix<double> Ai(35,35);
  crv::invertMatrix(35,A,Ai);
  mth::Matrix<double> eye(35,35);
  mth::multiply(A,Ai,eye);
  for (int i = 0; i < 35; ++i)
    for (int j = 0; j < 35; ++j){
      assert(fabs(Ai(i,j) - a_invdata[i][j]) < 1e-14);
      assert(fabs(eye(i,j) - (i == j)) < 1e-14);
    }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  testEdgeSubdivision();
  testEdgeElevation();
  testTriSubdivision1();
  testTriSubdivision4();
  testTriElevation();
  testTetElevation();
  testNodeIndexing();
  testMatrixInverse();
  PCU_Comm_Free();
  MPI_Finalize();
}

