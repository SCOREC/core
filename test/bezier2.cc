#include <crv.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>

#include <math.h>

/* This test file uses an alternative and more traditional method to
 * compute Jacobian differences, using the property of Bezier's that
 * derivatives of Bezier's are themselves Bezier's multiplied by
 * control point differences. As implementing weights*control point differences
 * is not possible in our current framework, this method is not implemented in
 * the main code, but serves its use for code validation.
 */
static int getPointIndex(int P, int I, int J)
{
  static apf::NewArray<int> map[6];
  int m1[] = {2,0,1};
  int m2[] = {2,5,0,4,3,1};
  int m3[] = {2,7,8,0,6,9,3,5,4,1};
  int m4[] = {2,9,10,11,0,8,14,12,3,7,13,4,6,5,1};
  int m5[] = {2,11,12,13,14,0,10,19,20,15,3,9,18,16,4,8,17,5,7,6,1};
  int m6[] = {2,13,14,15,16,17,0,12,24,25,26,18,3,11,23,27,19,4,10,
      22,20,5,9,21,6,8,7,1};

  int* maps[6] = {m1,m2,m3,m4,m5,m6};;
  for(int j = 1; j <= 6; ++j){
    map[j-1].allocate((j+1)*(j+2)/2);
    for(int i = 0; i < (j+1)*(j+2)/2; ++i)
      map[j-1][i] = maps[j-1][i];
  }
  return map[P-1][J*(P+1)+I-J*(J-1)/2];
}

static double getJacobianDeterminant(apf::NewArray<apf::Vector3>& nodes,
    int d, int i1, int j1, int i2, int j2)
{
  int p00 = getPointIndex(d,i1+1,j1);
  int p01 = getPointIndex(d,i1,j1+1);
  int p10 = getPointIndex(d,i2+1,j2);
  int p11 = getPointIndex(d,i2,j2);
  apf::Vector3 c0 = nodes[p01]-nodes[p00];
  apf::Vector3 c1 = nodes[p11]-nodes[p10];
  return c0[0]*c1[1]-c1[0]*c0[1];
}

static int factorial(int i)
{
  static int table[13] = {1,1,2,6,24,120,720,5040,40320,362880,
      3628800,39916800,479001600};
  return table[i];
}

static int CDIJK(int d, int i, int j, int k)
{
  return factorial(d)/factorial(i)/factorial(j)/factorial(k);
}

static double NIJK(apf::NewArray<apf::Vector3>& nodes,
    int d, int I, int J, int K)
{
  int CD = CDIJK(2*(d-1),I,J,K);
  double sum = 0.;
  for(int j1 = 0; j1 <= J; ++j1){
    for(int i1 = 0; i1 <= I; ++i1){
      int k1 = d-1 - i1 - j1;
      if(i1 > d-1 || j1 > d-1 || I-i1 > d-1 || J-j1 > d-1 ||
         i1+j1 > d-1 || I-i1 + J-j1 > d-1)
        continue;
      sum += CDIJK(d-1,i1,j1,k1)
          *CDIJK(d-1,I-i1,J-j1,K-k1)
          *getJacobianDeterminant(nodes,d,i1,j1,I-i1,J-j1)/CD;
    }
  }
  return sum*d*d;
}

static double BIJK(int d, int i, int j, apf::Vector3& xi)
{
  double xii[3] = {1.-xi[0]-xi[1],xi[0],xi[1]};

  return factorial(d)/factorial(i)/factorial(j)/factorial(d-i-j)
      *pow(xii[0],i)*pow(xii[1],j)*pow(xii[2],d-i-j);
}

static double testTriangleJacobian(apf::Mesh* m, apf::MeshEntity* e,
    apf::Vector3& xi)
{
  int d = m->getShape()->getOrder();
  apf::Element* elem = apf::createElement(m->getCoordinateField(),e);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);

  double detJ = 0.;
  for(int I = 0; I <= 2*(d-1); ++I){
    for(int J = 0; J <= 2*(d-1)-I; ++J){
      detJ += BIJK(2*(d-1),I,J,xi)*NIJK(nodes,d,I,J,2*(d-1)-I-J);
    }
  }
  apf::destroyElement(elem);
  return detJ;
}

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
        double J = testTriangleJacobian(m,e,xi);
        assert(fabs(detJ-J) < 1e-14);
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);
}

static double BIJ(int d, int i, apf::Vector3& xi)
{
  return factorial(d)/factorial(i)/factorial(d-i)
      *pow(xi[0],i)*pow(1.-xi[0],d-i);
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
        J += (nodes[map[p+1]]-nodes[map[p]])*BIJ(d-1,p,xi)*d;
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
  x[0] = 1.0-10.0*p[0]*(p[0]-1.0)*p[0]*(p[0]-1.0);
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

static void fixMidPoints(apf::Mesh2* m)
{
  apf::FieldShape * fs = m->getShape();

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 xi;

  while ((e = m->iterate(it))) {
    apf::MeshEntity* v[3];
    m->getDownward(e,0,v);
    for(int i = 0; i < fs->countNodesOn(2); ++i){
      apf::Vector3 pt(0,0,0), vertex;

      fs->getNodeXi(apf::Mesh::TRIANGLE,i,xi);
      xi[2] = 1.-xi[0]-xi[1];
      for(int j = 0; j < 3; ++j){
        m->getPoint(v[j],0,vertex);
        pt += vertex*xi[j];
      }
      m->setPoint(e,i,pt);
    }
  }
  m->end(it);
}

void test2D()
{
  for(int order = 1; order <= 6; ++order){
    for(int blendOrder = 0; blendOrder <= 0; ++blendOrder){
      apf::Mesh2* m = createMesh2D();
      crv::BezierCurver bc(m,order,blendOrder);
      bc.run();
      if(order > 2)
        fixMidPoints(m);
      crv::writeCurvedVtuFiles(m,apf::Mesh::VERTEX,order,"curvedBezier2D");
      crv::writeCurvedVtuFiles(m,apf::Mesh::EDGE,order,"curvedBezier2D");
      crv::writeCurvedVtuFiles(m,apf::Mesh::TRIANGLE,order,"curvedBezier2D");
      testJacobian(m);
      testEdgeGradients(m);
      m->destroyNative();
      apf::destroyMesh(m);
    }
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  test2D();
  PCU_Comm_Free();
  MPI_Finalize();
}
