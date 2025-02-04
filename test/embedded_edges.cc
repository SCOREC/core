#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <cassert>
#include <iostream>
// this test checks that the destruct function works
// with meshes that have a lower dimension than the manifold
// which tey reside in. E.g. a truss or beam in 3D space

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==2);
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_mesh();
  gmi_register_null();
  apf::Gid* conn;
  double* coords;
  int nelem;
  int etype;
  int nverts;

  gmi_model* model = gmi_load(".null");
  apf::Mesh2* m = apf::loadMdsMesh(model, argv[1], &PCUObj);
  apf::deriveMdsModel(m);
  int dim = m->getDimension();
  extractCoords(m, coords, nverts);
  destruct(m, conn, nelem, etype, 1);
  assert(etype == apf::Mesh::EDGE);
  m->destroyNative();
  apf::destroyMesh(m);

  std::cout<<m->typeDimension[apf::Mesh::EDGE]<<std::endl;

  gmi_model* model2 = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(model2, dim, false, &PCUObj);
  apf::GlobalToVert outMap;
  apf::construct(m, conn, nelem, etype, outMap);
  delete [] conn;
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, coords, nverts, outMap);
  delete [] coords;
  outMap.clear();

  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

