#include <crv.h>
#include <crvBezier.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apf.h>
#include <PCU.h>
#include <apfDynamicMatrix.h>
#include <cassert>
#include <cstdlib>

/* this test file contains tests for
 * a curved 2D surface mesh
 * and a 3D Planar tetrahedron
 *
 * The basic idea is generate a mesh,
 * curve it, test it, and destroy it,
 * and do this for all orders
 */

/*
 * This analytic surface is a "pringle",
 * defined on [0,1] in R^2 as z = ((2x-1)^2-(2y-1)^2)*exp(xy)
 * matlab >> ezsurf('((2*x-1)^2-(2*y-1)^2)*exp(x*y)',[0,1],[0,1])
 */
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
  x[2] = (2.*p[0]-1.)*(2.*p[0]-1.)-1.;
}
void edge1(double const p[2], double x[3], void*)
{
  x[0] = 1;
  x[1] = p[0];
  x[2] = (1.-(2.*p[0]-1.)*(2.*p[0]-1.))*std::exp(p[0]);
}
void edge2(double const p[2], double x[3], void*)
{
  double u = 1.-p[0];
  x[0] = u;
  x[1] = 1.;
  x[2] = ((2.*u-1.)*(2.*u-1.)-1.)*std::exp(u);
}
void edge3(double const p[2], double x[3], void*)
{
  double v = 1.-p[0];
  x[0] = 0;
  x[1] = v;
  x[2] = 1.-(2.*v-1.)*(2.*v-1.);
}
void face0(double const p[2], double x[3], void*)
{
  x[0] = p[0];
  x[1] = p[1];
  x[2] = ((2.*p[0]-1.)*(2.*p[0]-1.)
      -(2.*p[1]-1.)*(2.*p[1]-1.))
      *std::exp(p[0]*p[1]);
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

/* test whether the bezier edge interpolates the bounding geometric edges
 * at the prescribed nodes. (note that the nodes are control points, not
 * interpolating points)
 */
void testInterpolatedPoints2D(apf::Mesh2* m){
  apf::FieldShape * fs = m->getShape();
  apf::MeshIterator* it = m->begin(1);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* g = m->toModel(e);
    if (m->getModelType(g) == m->getDimension())
      continue;
    apf::Vector3 pt,pa(0.,0.,0.),cpt,cpa;
    apf::Element* elem =
        apf::createElement(m->getCoordinateField(),e);

    for(int i = 0; i < fs->countNodesOn(1); ++i){
      fs->getNodeXi(apf::Mesh::EDGE,i,pa);
      apf::getVector(elem,pa,pt);
      pa[0] = 0.5*(pa[0]+1.);
      m->snapToModel(g,pa,cpt);
      double error = (cpt-pt).getLength();
      if (error > 1.e-7) {
        std::cout << apf::Mesh::typeName[m->getType(e)] <<
            " is not being interpolated correctly by " <<
            fs->getName() << " " << fs->getOrder() <<
            " with error " << error << "\n";
        abort();
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);
}

/* test the sizes of the edges and the face using measure
 * Uses a tolerance to perform a sort of convergence check, that is
 * ensure that the % error in the order,
 * (measured - exact)/exact < exp(-order*C) where C > 0.
 * As order increases, this approaches zero,
 * tightening the restriction on relative error.
 */
void testSize2D(apf::Mesh2* m, int order)
{
  // these are the exact sizes for the outer edges, middle edge, and surface
  double sizes[3] = {3.721402252672,2.323391881216468,2.133984157270};
  for(int d = 1; d <= 2; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::ModelEntity* g = m->toModel(e);
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      // these are checks for the middle edge and first order, which does way
      // worse than second order and higher
      bool correct = true;
      if (m->getModelType(g) == 2 && d == 1){
        correct = (std::fabs(v-1.414213562373) < 1e-10);
      } else if(order == 1 && d == 1){
        correct = (std::fabs(v-1.0) < 1e-10);
      } else if(order == 1 && d == 2){
        correct = (std::fabs(v-0.5) < 1e-10);
      } else if(d == 1){ // the middle edge is a difference size
        correct = ((std::fabs(v-sizes[0])/sizes[0] < std::exp(-0.5*order))
            || (std::fabs(v-sizes[1])/sizes[1] < std::exp(-0.5*order)));
      } else if(d == 2){
        correct = (std::fabs(v-sizes[2])/sizes[2] < std::exp(-0.25*order));
      }
      if(!correct){
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        abort();
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
}

void test2D()
{
  // test all orders for all blending orders
  for(int order = 1; order <= 6; ++order){
    for(int blendOrder = 0; blendOrder <= 3; ++blendOrder){
      apf::Mesh2* m = createMesh2D();
      crv::BezierCurver bc(m,order,blendOrder);
      bc.run();
      testInterpolatedPoints2D(m);
      testSize2D(m,order);
      m->destroyNative();
      apf::destroyMesh(m);
    }
  }
  // test the pure interpolation side of things
  {
    apf::Mesh2* m = createMesh2D();
    apf::changeMeshShape(m,apf::getLagrange(2),true);
    crv::InterpolatingCurver ic(m,2);
    ic.run();
    testInterpolatedPoints2D(m);
    testSize2D(m,2);

    m->destroyNative();
    apf::destroyMesh(m);
  }

}
/* Similar to 2D test, but since tetrahedron is planar,
 * more exact result required
 */
void testSize3D(apf::Mesh2* m)
{
  for(int d = 1; d <= 3; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      // sizes for edges, faces, and volume.
      if((d == 1 && !(std::fabs(v-1.) < 1e-13 || std::fabs(v-1.414213562373) < 1e-11)) ||
         (d == 2 && !(std::fabs(v-0.5) < 1e-13 || std::fabs(v-0.866025403780) < 1e-11)) ||
         (d == 3 && !(std::fabs(v-1./6.) < 1e-13)))
      {
        std::stringstream ss;
        ss << "error: " << apf::Mesh::typeName[m->getType(e)]
           << " size " << v
           << " at " << getLinearCentroid(m, e) << '\n';
        std::string s = ss.str();
        fprintf(stderr, "%s", s.c_str());
        abort();
      }
      apf::destroyMeshElement(me);
    }
    m->end(it);
  }
}
/* This tests the 3D Jacobian of the triangle. This isn't that meaningful,
 * but it ensures the Jacobian calculation reduces to the exact same result
 * as the linear tetrahedron would produce,
 * Exact Jacobian = [[v1_x-v0_x,v1_y-v0_y,v1_z-v0_z]
 *                   [v2_x-v0_x,v2_y-v0_y,v2_z-v0_z]
 *                   [0,0,0]]
 * checked at (n+1)*(n+2)/2 evenly spaced locations
 */

void test3DJacobianTri(apf::Mesh2* m)
{
  int n = 2;

  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e;
  apf::Vector3 xi,pt;
  apf::Matrix3x3 J;

  while ((e = m->iterate(it))) {
    apf::MeshElement* elem =
        apf::createMeshElement(m,e);
    apf::MeshEntity* v[3];
    apf::Vector3 vpt[3];
    m->getDownward(e,0,v);
    for(int i = 0; i < 3; ++i)
      m->getPoint(v[i],0,vpt[i]);

    apf::Matrix3x3 Je(vpt[1][0]-vpt[0][0],vpt[1][1]-vpt[0][1],vpt[1][2]-vpt[0][2],
        vpt[2][0]-vpt[0][0],vpt[2][1]-vpt[0][1],vpt[2][2]-vpt[0][2],0.,0.,0.);

    for (int j = 0; j <= n; ++j){
      xi[1] = 1.*j/n;
      for (int i = 0; i <= n-j; ++i){
        xi[0] = 1.*i/n;
        apf::getJacobian(elem,xi,J);
        for(int a = 0; a < 3; ++a)
          for(int b = 0; b < 3; ++b)
            if(std::fabs(Je[a][b]-J[a][b]) > 1.e-12)
            {
              std::cout << "xi " << xi << std::endl;
              std::cout << J << std::endl;
              std::cout << Je << std::endl;
              std::stringstream ss;
              ss << "2D Jacobian is incorrect for "
                  << m->getShape()->getName() << m->getShape()->getOrder()
                  << "\n";
              std::string s = ss.str();
              fprintf(stderr, "%s", s.c_str());
              abort();
            }
      }
    }

    apf::destroyMeshElement(elem);
  }
  m->end(it);
}
/* tests the full Jacobian, ensuring it reduces to the Jacobian of the linear
 * tetrahedron,
 * Exact Jacobian = [[v1_x-v0_x,v1_y-v0_y,v1_z-v0_z]
 *                   [v2_x-v0_x,v2_y-v0_y,v2_z-v0_z]
 *                   [v3_x-v0_x,v3_y-v0_y,v3_z-v0_z]]
 * checked at (n+1)*(n+2)*(n+3)/6 evenly spaced locations
 */
void test3DJacobian(apf::Mesh2* m)
{
  int n = 2;

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::Vector3 xi,pt;
  apf::Matrix3x3 J;
  while ((e = m->iterate(it))) {
    apf::MeshElement* elem =
        apf::createMeshElement(m,e);
    apf::MeshEntity* v[4];
    apf::Vector3 vpt[4];
    m->getDownward(e,0,v);
    for(int i = 0; i < 4; ++i)
      m->getPoint(v[i],0,vpt[i]);
    apf::Matrix3x3 Je(vpt[1][0]-vpt[0][0],vpt[1][1]-vpt[0][1],vpt[1][2]-vpt[0][2],
        vpt[2][0]-vpt[0][0],vpt[2][1]-vpt[0][1],vpt[2][2]-vpt[0][2],
        vpt[3][0]-vpt[0][0],vpt[3][1]-vpt[0][1],vpt[3][2]-vpt[0][2]);

    for (int k = 0; k <= n; ++k){
      xi[2] = 1.*k/n;
      for (int j = 0; j <= n-k; ++j){
        xi[1] = 1.*j/n;
        for (int i = 0; i <= n-j-k; ++i){
          xi[0] = 1.*i/n;
          apf::getJacobian(elem,xi,J);
          for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
              if(std::fabs(Je[a][b]-J[a][b]) > 1.e-12)
              {
                std::cout << "xi " << xi << std::endl;
                std::cout << J << std::endl;
                std::cout << Je << std::endl;
                std::stringstream ss;
                ss << "3D Jacobian is incorrect for "
                    << m->getShape()->getName() << m->getShape()->getOrder()
                    <<"\n";
                std::string s = ss.str();
                fprintf(stderr, "%s", s.c_str());
                abort();
              }
        }
      }
    }
    apf::destroyMeshElement(elem);
  }
  m->end(it);
}

/* tests an alternate method for calculating the Jacobian,
 * by creating the Bezier representation of the Jacobian, and
 * evaluating it. This is the basis for the validity checks
 */
static void testAlternateTetJacobian(apf::Mesh2* m)
{
  int n = 2;

  apf::MeshIterator* it = m->begin(3);
  apf::MeshEntity* e;
  apf::Vector3 xi;
  apf::Matrix3x3 Jac;

  while ((e = m->iterate(it))) {
    apf::MeshElement* me =
        apf::createMeshElement(m,e);
    for (int k = 0; k <= n; ++k){
      xi[2] = 1.*k/n;
      for (int j = 0; j <= n-k; ++j){
        xi[1] = 1.*j/n;
        for (int i = 0; i <= n-j-k; ++i){
          xi[0] = 1.*i/n;
          apf::getJacobian(me,xi,Jac);
          double detJ = apf::getDeterminant(Jac);
          double J2 = crv::computeTetJacobianDetFromBezierFormulation(m,e,xi);
          assert(std::fabs(detJ-J2) < 1e-13);
        }
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);
}

/* Tests 3D with blending functions. This can go as high as the order of
 * triangles implemented. There are no nodes inside the tetrahedron
 *
 */
void test3DBlended()
{
  gmi_register_null();
  apf::Mesh2* mbase = createMesh3D();

  testSize3D(mbase);
  test3DJacobian(mbase);
  test3DJacobianTri(mbase);

  mbase->destroyNative();
  apf::destroyMesh(mbase);

  for(int order = 1; order <= 6; ++order){
    for(int blendOrder = 1; blendOrder <= 3; ++blendOrder){
      apf::Mesh2* m = createMesh3D();
      apf::changeMeshShape(m, crv::getBezier(order),true);
      crv::setBlendingOrder(apf::Mesh::TYPES,blendOrder);
      apf::FieldShape * fs = m->getShape();
      crv::BezierCurver bc(m,order,blendOrder);
      // go downward, and convert interpolating to control points
      for(int d = 2; d >= 1; --d){
        int n = fs->getEntityShape(apf::Mesh::simplexTypes[d])->countNodes();
        int ni = fs->countNodesOn(d);
        if(ni <= 0) continue;

        apf::NewArray<double> c;
        crv::getBezierTransformationCoefficients(order,d,c);
        apf::MeshEntity* e;
        apf::MeshIterator* it = m->begin(d);
        while ((e = m->iterate(it))) {
          if(m->getModelType(m->toModel(e)) == m->getDimension()) continue;
          crv::convertInterpolationPoints(m,e,n,ni,c);
        }
        m->end(it);
      }
      m->acceptChanges();
      testSize3D(m);
      test3DJacobian(m);
      test3DJacobianTri(m);
      m->destroyNative();
      apf::destroyMesh(m);
    }
  }
}

/* Tests Full Bezier tetrahedra */
void test3DFull()
{
  gmi_register_null();

  for(int order = 1; order <= 6; ++order){
    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, crv::getBezier(order),true);
    apf::FieldShape* fs = m->getShape();
    crv::setBlendingOrder(apf::Mesh::TYPES,0);
    crv::BezierCurver bc(m,order,0);
    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
      int ni = fs->countNodesOn(d);
      if(ni <= 0) continue;
      apf::NewArray<double> c;
      crv::getBezierTransformationCoefficients(order,d,c);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        crv::convertInterpolationPoints(m,e,n,ni,c);
      }
      m->end(it);
    }
    if(!fs->hasNodesIn(3)) continue;
    int n = fs->getEntityShape(apf::Mesh::simplexTypes[3])->countNodes();
    int ne = fs->countNodesOn(apf::Mesh::simplexTypes[3]);
    apf::NewArray<double> c;
    crv::getInternalBezierTransformationCoefficients(m,order,1,
        apf::Mesh::simplexTypes[3],c);
    apf::MeshEntity* e;
    apf::MeshIterator* it = m->begin(3);
    while ((e = m->iterate(it))){
      crv::convertInterpolationPoints(m,e,n-ne,ne,c);
    }
    m->end(it);

    m->acceptChanges();
    if(order <= 6)
      testSize3D(m);
    test3DJacobian(m);
    test3DJacobianTri(m);
    if(order <= 8){
      testAlternateTetJacobian(m);
    }

    if(order == 4){
      // put a field on the mesh to make sure nothing fails

      apf::Field* f1 =
          apf::createField(m,"field1",apf::SCALAR,apf::getLagrange(2));
      apf::Field* f2 =
          apf::createField(m,"field2",apf::VECTOR,apf::getLagrange(2));
      apf::Field* f3 =
          apf::createField(m,"field3",apf::MATRIX,apf::getLagrange(2));

      for(int d = 0; d <= 3; ++d){
        int ne = apf::getLagrange(2)->countNodesOn(apf::Mesh::simplexTypes[d]);
        apf::MeshEntity* e;
        apf::MeshIterator* it = m->begin(d);
        while ((e = m->iterate(it))) {
          for(int i = 0; i < ne; ++i){
            apf::setScalar(f1,e,i,1.2345);
            apf::Vector3 V(1.2345,1.2345,1.2345);
            apf::Matrix3x3 M(1.2345,1.2345,1.2345,
                1.2345,1.2345,1.2345,
                1.2345,1.2345,1.2345);

            apf::setVector(f2,e,i,V);
            apf::setMatrix(f3,e,i,M);
          }
        }
        m->end(it);
      }

      // write the field
//      crv::writeCurvedVtuFiles(m,apf::Mesh::EDGE,2,"curved");
      crv::writeCurvedVtuFiles(m,apf::Mesh::TET,2,"curved");
    }
//    crv::writeControlPointVtuFiles(m,"curved");
    m->destroyNative();
    apf::destroyMesh(m);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  test2D();
  test3DBlended();
  test3DFull();
  PCU_Comm_Free();
  MPI_Finalize();
}
