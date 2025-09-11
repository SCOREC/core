#include <apf.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <lionPrint.h>
#ifdef PUMI_HAS_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <cstdlib>

static void number_dim(apf::Mesh* m, apf::FieldShape* shape, int dim, std::string const& prefix) {
  std::string name = prefix + "_class_dim";
  apf::Numbering* cdn = apf::createNumbering(m, name.c_str(),
      shape, 1);
  name = prefix + "_class_id";
  apf::Numbering* cin = apf::createNumbering(m, name.c_str(),
      shape, 1);
  apf::MeshIterator* it = m->begin(dim);
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::ModelEntity* g = m->toModel(e);
    apf::number(cdn, e, 0, 0, m->getModelType(g));
    apf::number(cin, e, 0, 0, m->getModelTag(g));
  }
  m->end(it);
}

int main(int argc, char** argv)
{
  pcu::Init(&argc,&argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  if (!(argc == 4 || argc == 5)) {
    if ( !PCUObj.Self() ) {
      printf("Usage: %s <model> <mesh> <out prefix>\n", argv[0]);
      printf("       %s <model> <mesh> <dim> <out prefix>\n", argv[0]);
    }
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef PUMI_HAS_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  char const* modelpath = argv[1];
  char const* meshpath = argv[2];
  apf::Mesh2* m = apf::loadMdsMesh(modelpath, meshpath, &PCUObj);
  int dim;
  char const* vtkpath;
  if (argc == 5) {
    dim = atoi(argv[3]);
    vtkpath = argv[4];
  } else {
    dim = m->getDimension();
    vtkpath = argv[3];
  }
  number_dim(m, m->getShape(), 0, "vert");
  number_dim(m, apf::getConstant(dim), dim, "cell");
  apf::writeVtkFiles(vtkpath, m, dim);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef PUMI_HAS_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}


