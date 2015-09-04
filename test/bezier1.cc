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

/*
 * This analytic function is a "pringle",
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
      fs->getNodeXi(1,i,pa);
      apf::getVector(elem,pa,pt);
      pa[0] = 0.5*(pa[0]+1.);
      m->snapToModel(g,pa,cpt);
      double error = (cpt-pt).getLength();
      if (error > 1.e-7) {
        std::cout << apf::Mesh::typeName[m->getType(e)] <<
            " is not being interpolated correctly by " <<
            fs->getName() <<
            " with error " << error << "\n";
        abort();
      }
    }
    apf::destroyElement(elem);
  }
  m->end(it);
}

void testSize2D(apf::Mesh2* m)
{
  double sizes[3] = {3.721402252672,2.323391881216468,2.133984157270};
  int order = m->getShape()->getOrder();
  for(int d = 1; d <= 2; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::ModelEntity* g = m->toModel(e);
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      //			      printf("%s has size %.16f\n",apf::Mesh::typeName[m->getType(e)],v);
      // these are checks for the middle edge and first order
      bool correct = true;
      if (m->getModelType(g) == 2 && d == 1){
        correct = (fabs(v-1.414213562373)< 1e-10);
      } else if(order == 1 && d == 1){
        correct = (fabs(v-1.0 )< 1e-10);
      } else if(order == 1 && d == 2){
        correct = (fabs(v-0.5)< 1e-10);
      } else if(d == 1){        // lazy convergence check
        correct = ((fabs(v-sizes[0])/sizes[0] < (1.-0.5*exp(-(double)order)))
            || (fabs(v-sizes[1])/sizes[1] < (1.-0.5*exp(-(double)order))));
      } else if(d == 2){
        correct = (fabs(v-sizes[2])/sizes[2] < (1.-0.5*exp(-(double)order)));
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
  for(int order = 1; order <= 6; ++order){
    for(int blendOrder = 0; blendOrder <= 3; ++blendOrder){
      apf::Mesh2* m = createMesh2D();
      crv::BezierCurver bc(m,order,blendOrder);
      bc.run();
      testInterpolatedPoints2D(m);
      testSize2D(m);
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
    testSize2D(m);

    m->destroyNative();
    apf::destroyMesh(m);
  }

}

void testSize3D(apf::Mesh2* m)
{
  for(int d = 1; d <= 3; ++d){
    apf::MeshIterator* it = m->begin(d);
    apf::MeshEntity* e;
    while ((e = m->iterate(it))) {
      apf::MeshElement* me = apf::createMeshElement(m,e);
      double v = apf::measure(me);
      //      printf("%s has size %.16f\n",apf::Mesh::typeName[m->getType(e)],v);
      if((d == 1 && !(fabs(v-1.) < 1e-8  || fabs(v-1.414213562373) < 1e-8)) ||
          (d == 2 && !(fabs(v-0.5) < 1e-8 || fabs(v-0.866025403780) < 1e-8)) ||
          (d == 3 && !(fabs(v-1./6) < 1e-8)))
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

void test3DJacobianTri(apf::Mesh2* m)
{
  int n = 3;

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
            if(fabs(Je[a][b]-J[a][b]) > 1.e-13)
            {
              std::cout << "xi " << xi << std::endl;
              std::cout << J << std::endl;
              std::cout << Je << std::endl;
              std::stringstream ss;
              ss << "2D Jacobian is incorrect for "
                  << m->getShape()->getName() << "\n";
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

void test3DJacobian(apf::Mesh2* m)
{
  int n = 11;

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
//          if(crv::getBlendingOrder() == 0){
//          printf("%d %f %f\n",m->getShape()->getOrder(),apf::getJacobianDeterminant(J,3),
//              crv::computeAlternativeTetJacobianDet(m,e,xi));
//          assert(std::fabs(apf::getJacobianDeterminant(J,3)-
//              crv::computeAlternativeTetJacobianDet(m,e,xi)) < 1e-13);
//          }
          for(int a = 0; a < 3; ++a)
            for(int b = 0; b < 3; ++b)
              if(std::fabs(Je[a][b]-J[a][b]) > 1.e-13)
              {
                std::cout << "xi " << xi << std::endl;
                std::cout << J << std::endl;
                std::cout << Je << std::endl;
                std::stringstream ss;
                ss << "3D Jacobian is incorrect for "
                    << m->getShape()->getName() << "\n";
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

static void testAlternateTetJacobian(apf::Mesh2* m)
{
  int n = 5;

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
          double J2 = crv::computeAlternateTetJacobianDet(m,e,xi);
          assert(std::fabs(detJ-J2) < 1e-13);
        }
      }
    }
    apf::destroyMeshElement(me);
  }
  m->end(it);
}

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
      apf::changeMeshShape(m, crv::getBezier(3,order),true);
      crv::setBlendingOrder(blendOrder);
      apf::FieldShape * fs = m->getShape();
      crv::BezierCurver bc(m,order,blendOrder);
      // go downward, and convert interpolating to control points
      for(int d = 2; d >= 1; --d){
        int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
        int ne = fs->countNodesOn(d);
        apf::NewArray<double> c;
        crv::getTransformationCoefficients(3,order,d,c);
        apf::MeshEntity* e;
        apf::MeshIterator* it = m->begin(d);
        while ((e = m->iterate(it))) {
          if(m->getModelType(m->toModel(e)) == m->getDimension()) continue;
          bc.convertInterpolationPoints(e,n,ne,c);
        }
        m->end(it);
      }
      m->acceptChanges();
      testSize3D(m);
      if(blendOrder < 4){
        test3DJacobian(m);
        test3DJacobianTri(m);
      }
      m->destroyNative();
      apf::destroyMesh(m);
    }
  }
}

void test3DFull()
{
  gmi_register_null();

  for(int order = 1; order <= 4; ++order){
    apf::Mesh2* m = createMesh3D();
    apf::changeMeshShape(m, crv::getBezier(3,order),true);
    crv::setBlendingOrder(0);
    apf::FieldShape* fs = m->getShape();
    crv::BezierCurver bc(m,4,0);
    // go downward, and convert interpolating to control points
    for(int d = 2; d >= 1; --d){
      int n = (d == 2)? (order+1)*(order+2)/2 : order+1;
      int ne = fs->countNodesOn(d);
      apf::NewArray<double> c;
      crv::getTransformationCoefficients(3,order,d,c);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        if(m->getModelType(m->toModel(e)) == m->getDimension()) continue;
        bc.convertInterpolationPoints(e,n,ne,c);
      }
      m->end(it);
    }
    m->acceptChanges();
    testSize3D(m);
    test3DJacobian(m);
    test3DJacobianTri(m);
    testAlternateTetJacobian(m);
    crv::writeCurvedVtuFiles(m,apf::Mesh::TET,10,"curvedTest");
    // check values sum to 1.0 in the tet fieldshape
    apf::DynamicMatrix A;
    crv::getTransformationMatrix(m,apf::Mesh::TET,A);
    apf::EntityShape* es = fs->getEntityShape(apf::Mesh::TET);
    int n = es->countNodes();
    for(int i = 0; i < n; ++i){
      double sum = 0.;
      for(int j = 0; j < n; ++j)
        sum += A(i,j);
      assert(std::abs(sum - 1.0) < 1e-15);
    }

    // put a field on the mesh for fun
    int types[4] = {apf::Mesh::VERTEX, apf::Mesh::EDGE,
        apf::Mesh::TRIANGLE,apf::Mesh::TET};

    apf::Field* f1 =
        apf::createField(m,"field1",apf::SCALAR,apf::getLagrange(2));
    apf::Field* f2 =
        apf::createField(m,"field2",apf::VECTOR,apf::getLagrange(2));
    apf::Field* f3 =
        apf::createField(m,"field3",apf::MATRIX,apf::getLagrange(2));

    for(int d = 0; d <= 3; ++d){
      int ne = apf::getLagrange(2)->countNodesOn(types[d]);
      apf::MeshEntity* e;
      apf::MeshIterator* it = m->begin(d);
      while ((e = m->iterate(it))) {
        for(int i = 0; i < ne; ++i){
          apf::setScalar(f1,e,i,0.1*(rand() % 10));
          apf::Vector3 V(0.1*(rand() % 10),0.1*(rand() % 10),0.1*(rand() % 10));
          apf::Matrix3x3 M(0.1*(rand() % 10),0.1*(rand() % 10),
              0.1*(rand() % 10),0.1*(rand() % 10),0.1*(rand() % 10),
              0.1*(rand() % 10),0.1*(rand() % 10),0.1*(rand() % 10),
              0.1*(rand() % 10));

          apf::setVector(f2,e,i,V);
          apf::setMatrix(f3,e,i,M);
        }
      }
      m->end(it);
    }

    // write the field
    crv::writeCurvedVtuFiles(m,apf::Mesh::EDGE,5,"curved");
    crv::writeCurvedVtuFiles(m,apf::Mesh::TRIANGLE,5,"curved");
    crv::writeCurvedVtuFiles(m,apf::Mesh::TET,5,"curved");

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
