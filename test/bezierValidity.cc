#include <crv.h>
#include <crvBezier.h>
#include <crvQuality.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>

#include <math.h>
#include <cassert>

/* This test file uses an alternative and more traditional method to
 * compute Jacobian differences, using the property of Bezier's that
 * derivatives of Bezier's are themselves Bezier's multiplied by
 * control point differences. As implementing weights*control point differences
 * is not possible in our current framework, this method is not implemented in
 * the main code, but serves its use for code validation.
 *
 * This test file also contains validity checks.
 * In 2D, orders 3-6 provide invalid meshes
 */

static void testJacobian(apf::Mesh2* m)
{
  int n = 10;

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 xi;
  apf::Matrix3x3 Jac;

  while ((e = m->iterate(it))) {
    apf::MeshElement* me =
        apf::createMeshElement(m,e);
    for (int j = 0; j <= n; ++j){
      xi[1] = 1.*j/n;
      for (int i = 0; i <= n-j; ++i){
        xi[0] = 1.*i/n;
        apf::getJacobian(me,xi,Jac);
        double detJ = (Jac[0][0]*Jac[1][1])-(Jac[1][0]*Jac[0][1]);
        double J = crv::computeTriJacobianDetFromBezierFormulation(m,e,xi);
        assert(fabs(detJ-J) < 1e-14);
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);
}

static void testEdgeGradients(apf::Mesh2* m)
{
  int n = 5;
  int d = m->getShape()->getOrder();
  apf::NewArray<int> map(d+1);
  for(int p = 1; p < d; ++p){
    map[p] = p+1;
  }
  map[0] = 0;
  map[d] = 1;

  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  apf::Vector3 xi;
  apf::Matrix3x3 Jac;
  apf::NewArray<apf::Vector3> nodes;

  while ((e = m->iterate(it))) {
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);
    apf::getVectorNodes(elem,nodes);
    apf::MeshElement* me =
        apf::createMeshElement(m,e);

    for (int i = 0; i <= n; ++i){
      xi[0] = 2.*i/n-1.;
      apf::getJacobian(me,xi,Jac);
      apf::Vector3 J(0,0,0);
      xi[0] = (double)i/n;
      for(int p = 0; p <= d-1; ++p){
        J += (nodes[map[p+1]]-nodes[map[p]])*d
            *crv::binomial(d-1,p)*crv::Bij(p,d-1-p,xi[0],1.-xi[0]);
      }
      assert(fabs(J[0]-Jac[0][0]*2.0) < 1e-14
          && fabs(J[1]-Jac[0][1]*2.0) < 1e-14);
    }
    apf::destroyMeshElement(me);
    apf::destroyElement(elem);
  }
  m->end(it);
}

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
  x[1] = p[0]*(p[0]-1.0);
}
void edge1(double const p[2], double x[3], void*)
{
  x[0] = 1.0-5.0*p[0]*(p[0]-1.0)*p[0]*(p[0]-1.0);
  x[1] = p[0];
}
void edge2(double const p[2], double x[3], void*)
{
  double u = 1.-p[0];
  x[0] = u;
  x[1] = 1.;
}
void edge3(double const p[2], double x[3], void*)
{
  double v = 1.-p[0];
  x[0] = 0;
  x[1] = v;
}

void face0(double const p[2], double x[3], void*)
{
  x[0] = p[0];
  x[1] = p[1];
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

void make_edge_topo(gmi_model* m, gmi_ent* e, int v0tag, int v1tag)
{
  agm_bdry b = add_bdry(m, e);
  agm_use u0 = add_adj(m, b, v0tag);
  gmi_add_analytic_reparam(m, u0, reparam_zero, 0);
  agm_use u1 = add_adj(m, b, v1tag);
  gmi_add_analytic_reparam(m, u1, reparam_one, 0);
}

gmi_model* makeModel()
{
  gmi_model* model = gmi_make_analytic();
  int edPer = 0;
  double edRan[2] = {0, 1};
  for(int i = 0; i < 4; ++i)
    gmi_add_analytic(model, 0, i, vert0, NULL,NULL,NULL);
  gmi_ent* eds[4];
  eds[0] = gmi_add_analytic(model, 1, 0, edge0, &edPer, &edRan, 0);
  eds[1] = gmi_add_analytic(model, 1, 1, edge1, &edPer, &edRan, 0);
  eds[2] = gmi_add_analytic(model, 1, 2, edge2, &edPer, &edRan, 0);
  eds[3] = gmi_add_analytic(model, 1, 3, edge3, &edPer, &edRan, 0);
  for(int i = 0; i < 4; ++i)
    make_edge_topo(model, eds[i], i, (i+1) % 4);
  int faPer[2] = {0, 0};
  double faRan[2][2] = {{0,1},{0,1}};
  gmi_add_analytic(model, 2, 0, face0, faPer, faRan, 0);

  return model;
}

apf::Mesh2* createMesh2D()
{
  gmi_model* model = makeModel();
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 2, true);
  apf::MeshEntity* v[4];
  apf::Vector3 points2D[4] =
  {apf::Vector3(0,0,0),
      apf::Vector3(1,0,0),
      apf::Vector3(1,1,0),
      apf::Vector3(0,1,0)};

  for (int i = 0; i < 4; ++i){
    v[i] = m->createVertex(m->findModelEntity(0,i),points2D[i],points2D[i]);
  }
  for (int i = 0; i < 4; ++i){
    apf::ModelEntity* edge = m->findModelEntity(1,i);
    apf::MeshEntity* ved[2] = {v[i],v[(i+1) % 4]};
    apf::buildElement(m, edge, apf::Mesh::EDGE, ved);
  }

  apf::MeshEntity* ved[2] = {v[0],v[2]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::EDGE, ved);

  apf::MeshEntity* vf0[3] = {v[0],v[1],v[2]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, vf0);
  apf::MeshEntity* vf1[3] = {v[0],v[2],v[3]};
  apf::buildElement(m, m->findModelEntity(2,0), apf::Mesh::TRIANGLE, vf1);
  m->acceptChanges();
  m->verify();
  return m;
}
void checkEntityValidity(int numInvalid, int entity, int order)
{
  if(entity == 1){
    assert(numInvalid == 0);
  } else {
    assert((numInvalid && order != 3) || (numInvalid == 0 && order == 3));
  }
}

void checkValidity(apf::Mesh* m, int order)
{
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  int entityNum = 0;
  while ((e = m->iterate(it))) {
    apf::MeshEntity* entities[3];
    int numInvalid = crv::checkTriValidity(m,e,entities,2);
    checkEntityValidity(numInvalid,entityNum,order);
    numInvalid = crv::checkTriValidity(m,e,entities,3);
    checkEntityValidity(numInvalid,entityNum,order);
    numInvalid = crv::checkTriValidity(m,e,entities,4);
    checkEntityValidity(numInvalid,entityNum,order);

    entityNum++;


  }
  m->end(it);
}

void test2D()
{
  for(int order = 2; order <= 6; ++order){
      apf::Mesh2* m = createMesh2D();
      apf::changeMeshShape(m, crv::getBezier(order),true);
      crv::BezierCurver bc(m,order,0);
      // creates interpolation points based on the edges of the geometry
      bc.snapToInterpolate(1);
      apf::FieldShape* fs = m->getShape();

      // go downward, and convert interpolating to control points
      {
        int n = order+1;
        int ne = fs->countNodesOn(apf::Mesh::EDGE);
        apf::NewArray<double> c;
        crv::getBezierTransformationCoefficients(m,order,apf::Mesh::EDGE,c);
        apf::MeshEntity* e;
        apf::MeshIterator* it = m->begin(1);
        while ((e = m->iterate(it))) {
          bc.convertInterpolationPoints(e,n,ne,c);
        }
        m->end(it);
      }
      if(fs->hasNodesIn(2)) {
        int n = crv::getNumControlPoints(2,order);
        int ne = fs->countNodesOn(apf::Mesh::TRIANGLE);
        apf::NewArray<double> c;
        crv::getInternalBezierTransformationCoefficients(m,order,1,
            apf::Mesh::TRIANGLE,c);
        apf::MeshEntity* e;
        apf::MeshIterator* it = m->begin(2);
        while ((e = m->iterate(it))){
          bc.convertInterpolationPoints(e,n-ne,ne,c);
        }
        m->end(it);
      }
      crv::writeCurvedVtuFiles(m,apf::Mesh::TRIANGLE,50,"curved");
      crv::writeCurvedVtuFiles(m,apf::Mesh::EDGE,10,"curved");

      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(2);
      while ((e = m->iterate(it))){
        crv::getQuality(m,e);
        break;
      }
      m->end(it);

      testJacobian(m);
      testEdgeGradients(m);
      checkValidity(m,order);

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

void test3D()
{
  gmi_register_null();

  for(int order = 2; order <= 4; ++order){
    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, crv::getBezier(order),true);
    apf::FieldShape* fs = m->getShape();
    crv::BezierCurver bc(m,order,0);
    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
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


    // get face 2
    apf::MeshIterator* it = m->begin(3);
    apf::MeshEntity* tet = m->iterate(it);
    m->end(it);

    apf::MeshEntity* faces[4], *edges[3];
    m->getDownward(tet,2,faces);
    apf::MeshEntity* face = faces[2];
    m->getDownward(face,1,edges);
    for (int edge = 0; edge < 3; ++edge){
      int non = m->getShape()->countNodesOn(apf::Mesh::EDGE);
      for (int i = 0; i < non; ++i){
        apf::Vector3 pt;
        m->getPoint(edges[edge],i,pt);
        pt = pt*0.6;
        m->setPoint(edges[edge],i,pt);
      }
    }
    int non = m->getShape()->countNodesOn(apf::Mesh::TRIANGLE);
    for (int i = 0; i < non; ++i){
      apf::Vector3 pt;
      m->getPoint(face,i,pt);
      pt = pt*0.6;
      m->setPoint(face,i,pt);
    }
    for(int d = 2; d <= 3; ++d){
      if(!fs->hasNodesIn(d)) continue;
      int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
      int ne = fs->countNodesOn(apf::Mesh::simplexTypes[d]);
      apf::NewArray<double> c;
      crv::getInternalBezierTransformationCoefficients(m,order,1,
          apf::Mesh::simplexTypes[d],c);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))){
        bc.convertInterpolationPoints(e,n-ne,ne,c);
      }
      m->end(it);
    }
    m->acceptChanges();

    apf::MeshEntity* entities[14];
    int numInvalid = crv::checkTetValidity(m,tet,entities,2);

    if(order == 4){
      assert(numInvalid > 0);
    } else {
      assert(numInvalid == 0);
    }
    numInvalid = crv::checkTetValidity(m,tet,entities,3);
    crv::writeCurvedVtuFiles(m,apf::Mesh::TET,50,"curvedValidity");

    if(order == 4){
      assert(numInvalid > 0);
    } else {
      assert(numInvalid == 0);
    }
    crv::getQuality(m,tet);
    m->destroyNative();
    apf::destroyMesh(m);
  }

}
int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  test2D();
  test3D();
  PCU_Comm_Free();
  MPI_Finalize();
}
