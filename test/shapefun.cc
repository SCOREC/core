#include <apfShape.h>
#include <apfMesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>

void testNodeValues(int type, apf::Vector3 const* nodes, int nnodes)
{
  apf::EntityShape* shp =
    apf::getLagrange(1)->getEntityShape(type);
  assert( shp->countNodes() == nnodes );
  for (int i = 0; i < nnodes; ++i) {
    apf::NewArray<double> values;
    shp->getValues(nodes[i], values);
    for (int j = 0; j < nnodes; ++j) {
      if (j == i)
        assert(values[j] == 1);
      else
        assert(values[j] == 0);
    }
  }
}

void testPrismNodeValues()
{
  apf::Vector3 nodes[6];
  int i = 0;
  for (int xi2 = -1; xi2 < 2; xi2 += 2)
  for (int xi1 = 0; xi1 < 2; ++xi1)
  for (int xi0 = 0; xi0 < 2; ++xi0) {
    if (xi0 == 1 && xi1 == 1)
      continue;
    nodes[i] = apf::Vector3(xi0,xi1,xi2);
    ++i;
  }
  testNodeValues(apf::Mesh::PRISM, nodes, 6);
}

void testPyramidNodeValues()
{
  apf::Vector3 nodes[5] =
  {apf::Vector3(-1,-1,-1),
   apf::Vector3( 1,-1,-1),
   apf::Vector3( 1, 1,-1),
   apf::Vector3(-1, 1,-1),
   apf::Vector3( 0, 0, 1)};
  testNodeValues(apf::Mesh::PYRAMID, nodes, 5);
}

void testVolume(int type, apf::Vector3 const* points, double volume)
{
  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(model, 3, false);
  apf::MeshEntity* e = apf::buildOneElement(
      m, m->findModelEntity(3,0), type, points);
  apf::MeshElement* me = apf::createMeshElement(m, e);
  double v = apf::measure(me);
  apf::destroyMeshElement(me);
  assert( fabs(v - volume) < 1e-10 );
  m->destroyNative();
  apf::destroyMesh(m);
}

void testPrismVolume()
{
  apf::Vector3 points[6] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,1),
   apf::Vector3(1,0,1),
   apf::Vector3(0,1,1)};
  testVolume(apf::Mesh::PRISM, points, 0.5);
}

void testPyramidVolume()
{
  apf::Vector3 points[5] =
  {apf::Vector3(0,0,0),
   apf::Vector3(1,0,0),
   apf::Vector3(1,1,0),
   apf::Vector3(0,1,0),
   apf::Vector3(0,0,3)};
  testVolume(apf::Mesh::PYRAMID, points, 1);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  testPrismNodeValues();
  testPyramidNodeValues();
  testPrismVolume();
  testPyramidVolume();
  PCU_Comm_Free();
  MPI_Finalize();
}
