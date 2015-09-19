#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <PCU.h>
#include <cstdlib>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  apf::Numbering* cdn = apf::createNumbering(m, "class_dim",
      m->getShape(), 1);
  apf::Numbering* cin = apf::createNumbering(m, "class_id",
      m->getShape(), 1);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it))) {
    apf::ModelEntity* g = m->toModel(v);
    apf::number(cdn, v, 0, 0, m->getModelType(g));
    apf::number(cin, v, 0, 0, m->getModelTag(g));
  }
  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


