#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main (int argc, char* argv[]) {
  pcu::Init(&argc, &argv);
  { // pcu object scope
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  gmi_cap_start();
  gmi_register_cap();

  PCU_ALWAYS_ASSERT(argc == 2);
  const char* creFile(argv[1]);
  // 1. Load model.
  gmi_model* model = gmi_load(creFile);
  // 2. CreateMesh.
  apf::Mesh2* m = apf::createCapMesh(model, &PCUObj);
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
  gmi_cap_stop();
  } // pcu object scope
  pcu::Finalize();
}
