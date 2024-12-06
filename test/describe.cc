#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
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


static void print_stats(const char* name, double value, pcu::PCU *PCUObj)

{
  double min, max, avg;
  min = PCUObj->Min<double>(value);
  max = PCUObj->Max<double>(value);
  avg = PCUObj->Add<double>(value);
  avg /= PCUObj->Peers();
  double imb = max / avg;
  if (!PCUObj->Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

static void list_tags(apf::Mesh* m)
{
  if (m->getPCU()->Self())
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
  {
  pcu::PCU pcu_obj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  print_stats("kernal used before", pcu::GetMem(), &pcu_obj);
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);
  m->verify();
  print_stats("kernel heap", pcu::GetMem(), &pcu_obj);
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
  }
  MPI_Finalize();
}
