#include <apfShape.h>
#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>

void testPrismNodeValues()
{
  apf::EntityShape* prishp =
    apf::getLagrange(1)->getEntityShape(apf::Mesh::PRISM);
  assert( prishp->countNodes() == 6 );
  int i = 0;
  for (int xi2 = -1; xi2 < 2; xi2 += 2)
  for (int xi1 = 0; xi1 < 2; ++xi1)
  for (int xi0 = 0; xi0 < 2; ++xi0) {
    if (xi0 == 1 && xi1 == 1)
      continue;
    apf::NewArray<double> values;
    apf::Vector3 xi(xi0,xi1,xi2);
    prishp->getValues(xi, values);
    for (int j = 0; j < 6; ++j) {
      if (j == i)
        assert(values[j] == 1);
      else
        assert(values[j] == 0);
    }
    ++i;
  }
}

void testPrismVolume()
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false);
  apf::Vector3 points[6] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,1),
   apf::Vector3(1,0,1),
   apf::Vector3(0,1,1)};
  apf::MeshEntity* e = apf::buildOneElement(
      m, m->findModelEntity(3,0), apf::Mesh::PRISM, points);
  apf::MeshElement* me = apf::createMeshElement(m, e);
  double v = apf::measure(me);
  apf::destroyMeshElement(me);
  assert( fabs(v - 0.5) < 1e-10 );
  m->destroyNative();
  apf::destroyMesh(m);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  testPrismNodeValues();
  testPrismVolume();
  PCU_Comm_Free();
  MPI_Finalize();
}
