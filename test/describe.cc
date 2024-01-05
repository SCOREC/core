#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <lionPrint.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif
#include <pcu_util.h>

static void print_stats(const char* name, double value)
{
  double min, max, avg;
  min = PCU_Min_Double(value);
  max = PCU_Max_Double(value);
  avg = PCU_Add_Double(value);
  avg /= PCU_Comm_Peers();
  double imb = max / avg;
  if (!PCU_Comm_Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

static void list_tags(apf::Mesh* m)
{
  if (PCU_Comm_Self())
    return;
  apf::DynamicArray<apf::MeshTag*> tags;
  m->getTags(tags);
  for (size_t i = 0; i < tags.getSize(); ++i)
    printf("tag: \"%s\"\n",m->getTagName(tags[i]));
}

int main(int argc, char** argv)
{
  PCU_ALWAYS_ASSERT(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  print_stats("kernel used before", PCU_GetMem());
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  m->verify();
  print_stats("kernel heap", PCU_GetMem());
  Parma_PrintPtnStats(m, "");
  list_tags(m);
  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
