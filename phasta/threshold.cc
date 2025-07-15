#include <lionPrint.h>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <phRestart.h>
#include <gmi_mesh.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif

static double process_element(apf::Vector3 x[4], double sol[4][9])
{
  (void)sol;
  double volume = (apf::cross((x[1]-x[0]),(x[2]-x[0]))*(x[3]-x[0]))/6;
  return volume;
}

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  {
  pcu::PCU pcu_obj;
#ifdef HAVE_SIMMETRIX
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  lion_set_verbosity(1);
  ph::Input in;
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2], &pcu_obj);
  m->verify();
  in.restartFileName = argv[3];
  in.timeStepNumber = 0;
  in.ensa_dof = 9;
  ph::readAndAttachFields(in, m);
  apf::Field* solf = m->findField("solution");
  double total_volume = 0;
  apf::MeshIterator* it = m->begin(m->getDimension());
  apf::MeshEntity* e;
  while ((e = m->iterate(it))) {
    apf::MeshEntity* v[4];
    m->getDownward(e, 0, v);
    apf::Vector3 x[4];
    double sol[4][9];
    for (int i = 0; i < 4; ++i) {
      m->getPoint(v[i], 0, x[i]);
      apf::getComponents(solf, v[i], 0, sol[i]);
    }
    total_volume += process_element(x, sol);
  }
  m->end(it);
  lion_oprint(1,"volume %f\n", total_volume);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
#endif
  }
  pcu::Finalize();
}
