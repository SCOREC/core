#include <crv.h>
#include <crvBezier.h>
#include <crvSnap.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>
#include <ma.h>
#include <cassert>

/*
 * This analytic function is an edge in 3D
 */
void vertFunction(double const p[2], double x[3], void*)
{
  (void)p;
  (void)x;
}
// edges go counter clockwise
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

    apf::Vector3 p, pt;
    for (int i = 0; i <= 100; ++i){
      p[0] = 0.02*i-1.;
      apf::getVector(elem,p,pts[i]);
    }
    apf::destroyElement(elem);

    // elevate everything to 6th order
    apf::NewArray<apf::Vector3> elevatedNodes;
    elevatedNodes.allocate(7);
    crv::elevateBezierEdge(o,6-o,nodes,elevatedNodes);
    apf::changeMeshShape(m,crv::getBezier(3,6),false);

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
  x[0] = p[1]+1.-p[0]-p[1];
  x[1] = 1.-p[0]-p[1];
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
    gmi_add_analytic(model, 0, i, vertFunction, NULL,NULL,NULL);
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
void testTriSubdivision()
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
    m->destroy(e);

    apf::destroyElement(elem);
    for(int t = 0; t < 3; ++t)
      apf::destroyElement(elems[t]);

    m->destroyNative();
    apf::destroyMesh(m);
  }

  for(int o = 1; o <= 6; ++o){
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
    apf::MeshEntity* fourSplit[4];

    for (int t = 0; t < 4; ++t){
      apf::Vector3 points[3];
      for (int i = 0; i < 3; ++i)
        points[i] = subNodes[t][i];
      fourSplit[t] = apf::buildOneElement(m,m->findModelEntity(2,0),apf::Mesh::TRIANGLE,points);
      apf::MeshEntity* edges[3];
      m->getDownward(fourSplit[t],1,edges);
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < o-1; ++j){
          m->setPoint(edges[i],j,subNodes[t][3+i*(o-1)+j]);
        }
      for (int j = 0; j < (o-1)*(o-2)/2; ++j){
        m->setPoint(fourSplit[t],j,subNodes[t][3+3*(o-1)+j]);
      }
    }
    m->destroy(e);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

void testTriElevation()
{
  apf::Vector3 points2D[3] =
  {apf::Vector3(0,0,0),apf::Vector3(1,0,0),apf::Vector3(1,1,0)};

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
    apf::changeMeshShape(m,crv::getBezier(3,6),false);

    m->setPoint(verts[0],0,points2D[0]);
    m->setPoint(verts[1],0,points2D[1]);
    m->setPoint(verts[2],0,points2D[2]);

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

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  testEdgeSubdivision();
  testEdgeElevation();
  testTriSubdivision();
  testTriElevation();
  PCU_Comm_Free();
  MPI_Finalize();
}

