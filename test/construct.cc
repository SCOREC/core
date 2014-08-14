#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  PCU_Protect();
  gmi_register_mesh();
  gmi_register_null();
  int* conn;
  int nelem;
  int etype;
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  int dim = m->getDimension();
  destruct(m, conn, nelem, etype);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_model* model = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(model, dim, false);
  apf::construct(m, conn, nelem, etype);
  delete [] conn;
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  m->verify();
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

