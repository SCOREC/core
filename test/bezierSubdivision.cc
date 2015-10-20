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
 * This contains all the tests for bezier subdivision
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

/* Create a single tet, split it, and make sure each tet
 * exactly replicates its part of the old one
 *
 */
void testTetSubdivision1()
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

    apf::MeshEntity* verts[4],* edges[6],* faces[4];
    m->getDownward(tet,0,verts);
    m->getDownward(tet,1,edges);
    m->getDownward(tet,2,faces);
    apf::Element* elem = apf::createElement(m->getCoordinateField(),tet);

    apf::NewArray<apf::Vector3> nodes;
    apf::NewArray<apf::Vector3> subNodes[4];
    for (int t = 0; t < 4; ++t)
      subNodes[t].allocate(crv::getNumControlPoints(apf::Mesh::TET,order));

    apf::getVectorNodes(elem,nodes);
    apf::Vector3 splitpt(0.25,0.35,0.25);

    crv::subdivideBezierTet(order,splitpt,nodes,subNodes);

    apf::MeshEntity* v4 = m->createVertex(m->findModelEntity(2,0),
        subNodes[0][3],splitpt);

    apf::MeshEntity* newFaces[6],* newEdges[4],* newTets[4];
    for (int i = 0; i < 4; ++i){
      apf::MeshEntity* vE[2] = {verts[i],v4};
      newEdges[i] = m->createEntity(apf::Mesh::EDGE,m->findModelEntity(1,i),
          vE);
      for (int j = 0; j < order-1; ++j){
        m->setPoint(newEdges[i],j,subNodes[i][4+3*(order-1)+j]);
      }
    }
    int const tet_tri[6][2] =
    {{0,1},{0,2},{2,3},{3,1},{1,3},{2,1}};
    int nE = (order-1);
    int nF = (order-1)*(order-2)/2;
    // this compensates for alignment issues
    for (int f = 0; f < 6; ++f){
      apf::MeshEntity* eF[3] = {edges[f],newEdges[apf::tet_edge_verts[f][1]],
          newEdges[apf::tet_edge_verts[f][0]]};
      newFaces[f] = m->createEntity(apf::Mesh::TRIANGLE,m->findModelEntity(2,0),
          eF);
      for (int j = 0; j < nF; ++j){
        m->setPoint(newFaces[f],j,subNodes[tet_tri[f][0]][4+6*nE+tet_tri[f][1]*nF+j]);
        if(f == 3 && order == 4){
          int o[3] = {1,0,2};
          m->setPoint(newFaces[f],j,subNodes[tet_tri[f][0]][4+6*nE+tet_tri[f][1]*nF+o[j]]);
        }
      }
    }

    apf::MeshEntity* fT0[4] = {faces[0],newFaces[0],newFaces[1],newFaces[2]};
    newTets[0] = m->createEntity(apf::Mesh::TET,m->findModelEntity(3,0),fT0);
    apf::MeshEntity* fT1[4] = {faces[1],newFaces[0],newFaces[3],newFaces[4]};
    newTets[1] = m->createEntity(apf::Mesh::TET,m->findModelEntity(3,0),fT1);
    apf::MeshEntity* fT2[4] = {faces[2],newFaces[1],newFaces[4],newFaces[5]};
    newTets[2] = m->createEntity(apf::Mesh::TET,m->findModelEntity(3,0),fT2);
    apf::MeshEntity* fT3[4] = {faces[3],newFaces[2],newFaces[5],newFaces[3]};
    newTets[3] = m->createEntity(apf::Mesh::TET,m->findModelEntity(3,0),fT3);
    if(order == 4){
      for (int t = 0; t < 4; ++t){
        apf::Vector3 pt = (points3D[apf::tet_tri_verts[t][0]]
                    + points3D[apf::tet_tri_verts[t][1]]
                    + points3D[apf::tet_tri_verts[t][2]] + splitpt)*0.25;
        m->setPoint(newTets[t],0,pt);
      }
    }
    apf::Element* elems[4] =
      {apf::createElement(m->getCoordinateField(),newTets[0]),
       apf::createElement(m->getCoordinateField(),newTets[1]),
       apf::createElement(m->getCoordinateField(),newTets[2]),
       apf::createElement(m->getCoordinateField(),newTets[3])};
    double totalVolume = apf::measure((apf::MeshElement*)elem);
    double volumeSum = apf::measure((apf::MeshElement*)elems[0]) +
        apf::measure((apf::MeshElement*)elems[1]) +
        apf::measure((apf::MeshElement*)elems[2]) +
        apf::measure((apf::MeshElement*)elems[3]);

    assert(fabs(totalVolume - volumeSum) < 1e-15);
    apf::destroyElement(elem);
    m->destroy(tet);
    for(int t = 0; t < 4; ++t)
      apf::destroyElement(elems[t]);

    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  testEdgeSubdivision();
  testTriSubdivision1();
  testTriSubdivision4();
  testTetSubdivision1();
  PCU_Comm_Free();
  MPI_Finalize();
}
