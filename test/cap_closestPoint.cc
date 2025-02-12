#include <apf.h>
#include <apfCAP.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <CapstoneModule.h>
#include <PCU.h>

int main (int argc, char* argv[]) {
#ifndef SCOREC_NO_MPI
  MPI_Init(&argc, &argv);
#else
  (void) argc, (void) argv;
#endif
  pcu::PCU* PCUObj = new pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  gmi_register_cap();

  PCU_ALWAYS_ASSERT(argc == 2);
  std::string creFile(argv[1]);
  // 1. Load model.
  CapstoneModule cs("cap_inClosureOf", "Geometry Database : SMLIB",
    "Mesh Database : Create", "Attribution Database : Create");
  cs.load_files(v_string(1, creFile));
  // 2. CreateMesh.
  apf::Mesh2* m = apf::createMesh(cs.get_mesh(), cs.get_geometry(), PCUObj);

  PCU_ALWAYS_ASSERT(m->canGetClosestPoint());

  apf::ModelEntity* face = m->findModelEntity(2, 1);
  apf::Vector3 from(0.3, 0.4, -0.8), to, p;
  m->getClosestPoint(face, from, to, p);
  PCU_ALWAYS_ASSERT((to - apf::Vector3(0.3, 0.4, -0.5)).getLength() < 0.0001);

  apf::ModelEntity* edge = m->findModelEntity(1, 1);
  from = apf::Vector3(0.6, 0.34, 0.6);
  m->getClosestPoint(edge, from, to, p);
  PCU_ALWAYS_ASSERT((to - apf::Vector3(0.5, 0.34, 0.5)).getLength() < 0.0001);

  apf::destroyMesh(m);
  delete PCUObj;
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}
