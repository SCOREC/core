#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <lionPrint.h>
#include <cstdlib>
#include <string.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <in .msh> <out .smb> nMVskip \n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();
  int nMVskip=atoi(argv[4]);
  if(false) {
    apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(argv[1]), argv[2], nMVskip);
    // if input model is null derive a basic model for verify to pass.
    if (std::string(argv[1]).compare(".null") == 0)
      apf::deriveMdsModel(m);
  }
  apf::Mesh2* m = apf::loadMdsDmgFromGmsh(argv[1], argv[2], nMVskip);

  m->verify();
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

