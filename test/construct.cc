#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apf.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
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

  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&PCUObj);
  int dim = m->getDimension();
  extractCoords(m, coords, nverts);
  destruct(m, conn, nelem, etype);
  m->destroyNative();
  apf::destroyMesh(m);

  gmi_model* model = gmi_load(".null");
  m = apf::makeEmptyMdsMesh(model, dim, false, &PCUObj);
  apf::GlobalToVert outMap;
  apf::construct(m, conn, nelem, etype, outMap);
  delete [] conn;
  apf::alignMdsRemotes(m);
  apf::deriveMdsModel(m);
  apf::setCoords(m, coords, nverts, outMap);
  delete [] coords;
  outMap.clear();
  m->verify();

  //apf::writeVtkFiles("after", m);

  m->destroyNative();
  apf::destroyMesh(m);
  }
  MPI_Finalize();
}

