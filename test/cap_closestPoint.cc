#include <apf.h>
#include <apfMesh2.h>
#include <apfCAP.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <pcu_util.h>
#include <PCU.h>

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
  gmi_model* model = gmi_cap_load(creFile);
  // 2. CreateMesh.
  apf::Mesh2* m = apf::createCapMesh(model, &PCUObj);

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
  gmi_cap_stop();
  } // pcu object scope
  pcu::Finalize();
}
