#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <PCU.h>

namespace {

apf::Vector3 points[4] = 
{apf::Vector3(0,0,0),
 apf::Vector3(1,0,0),
 apf::Vector3(0,1,0),
 apf::Vector3(0,1,2)};

apf::Mesh2* createMesh()
{
  gmi_model* model = gmi_load(".null");
  return apf::makeEmptyMdsMesh(model,3,false);
}

apf::MeshElement* createTet(apf::Mesh2* m, apf::Vector3 const* points)
{
  apf::MeshEntity* e = apf::buildOneElement(
      m,m->findModelEntity(3,0),apf::Mesh::TET,points);
  return apf::createMeshElement(m,e);
}

void setFieldValues(apf::Field* f, int dim, double v)
{
  apf::Mesh* m = apf::getMesh(f);
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(dim);
  while ((e = m->iterate(it)))
    apf::setScalar(f,e,0,v);
}

apf::Field* createTestField(apf::Mesh2* m)
{
  apf::Field* f = apf::createHierarchicField(m,"test",apf::SCALAR);
  for (int d=0; d < 2; ++d)
    setFieldValues(f,d,1.0);
  return f;
}

void testInterpolation(apf::Element* e)
{
  double v = 0.350480947161671;
  apf::Vector3 p(0.25,0.25,0.25);
  assert( (apf::getScalar(e,p) - v) < 1e-12 );
}

void testShapeGrads(apf::Element* e)
{
  apf::Vector3 p(0.25,0.25,0.25);
  apf::NewArray<apf::Vector3> grads;
  apf::getShapeGrads(e,p,grads);
}

void testHierarchic(apf::Mesh2* m)
{
  apf::MeshElement* me = createTet(m,points);
  apf::Field* f = createTestField(m);
  apf::Element* e = apf::createElement(f,me);
  testInterpolation(e);
  testShapeGrads(e);
  apf::destroyField(f);
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_null();
  apf::Mesh2* m = createMesh();
  testHierarchic(m);
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

