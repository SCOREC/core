#include <apf.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#include <cstdlib>
#include <string.h>

apf::Field* convert_my_tag(apf::Mesh* m, std::string name) {
  apf::MeshTag* t = m->findTag(name.c_str());
  apf::MeshEntity* elm;
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::Field* f = apf::createField(m, name.c_str(), apf::SCALAR, 
                                   apf::getConstant(m->getDimension()));
  int vals[1];
  double vals_d[1];
  while ((elm = m->iterate(it))) {
    m->getIntTag(elm, t, vals);
    vals_d[0] = vals[0];
    apf::setComponents(f, elm, 0, vals_d);
  }
  m->end(it);
  return f;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  if ( argc != 4 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <in .msh> <out .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_null();
  gmi_register_mesh();
  apf::Mesh2* m = apf::loadMdsFromGmsh(gmi_load(argv[1]), argv[2]);
  // if input model is null derive a basic model for verify to pass.
  if (std::string(argv[1]).compare(".null") == 0)
    apf::deriveMdsModel(m);
  m->verify();
  convert_my_tag(m,"physical_type");
  apf::writeVtkFiles("foo.vtu", m);
  m->writeNative(argv[3]);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

