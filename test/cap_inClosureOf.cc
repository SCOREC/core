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
  pcu::PCU *PCUObj = new pcu::PCU(MPI_COMM_WORLD);
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
  // 3. Get region 1 ModelEntity*.
  apf::ModelEntity* rgn = m->findModelEntity(3, 1);
  PCU_ALWAYS_ASSERT(rgn);
  // 4. Assert each model entity is in the closure of that region.
  //FIXME: gmi_iter
  for (int d = 0; d < 3; ++d) {
    apf::MeshIterator* it = m->begin(d);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      apf::ModelEntity* ge = m->toModel(e);
      PCU_ALWAYS_ASSERT(m->isInClosureOf(ge, rgn));
    }
  }
  // 5. Test face 1 for edges 3, 5, 8, 12 and verts 3, 4, 7, 8.
  apf::ModelEntity* f1 = m->findModelEntity(2, 1);
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(1, 3), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(1, 5), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(1, 8), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(1, 12), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(0, 3), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(0, 4), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(0, 7), f1));
  PCU_ALWAYS_ASSERT(m->isInClosureOf(m->findModelEntity(0, 8), f1));

  apf::destroyMesh(m);
  delete PCUObj;
#ifndef SCOREC_NO_MPI
  MPI_Finalize();
#endif
}
