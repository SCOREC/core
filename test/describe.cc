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
#include <memory>

#ifdef __bgq__
#include <spi/include/kernel/memory.h>

static double get_peak(pcu::PCU*)
{
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return heap;
}

#elif defined (__linux__)

static double get_peak(pcu::PCU*)
{
#if defined(__GNUG__) && defined(PUMI_HAS_MALLINFO2)
  return mallinfo2().arena;
#elif defined(__GNUG__)
  return mallinfo().arena;
#endif
}

#else

static double get_peak(pcu::PCU *PCUObj)
{
  if(!PCUObj->Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}

#endif

static void print_stats(const char* name, double value, pcu::PCU *PCUObj)
{
  double min, max, avg;
  min = PCUObj->Min(value);
  max = PCUObj->Max(value);
  avg = PCUObj->Add(value);
  avg /= PCUObj->Peers();
  double imb = max / avg;
  if (!PCUObj->Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

#if defined(__linux__)

static double get_chunks(pcu::PCU*)
{
#if defined(__GNUG__) && defined(PUMI_HAS_MALLINFO2)
  struct mallinfo2 m = mallinfo2();
#elif defined(__GNUG__)
  struct mallinfo m = mallinfo();
#endif
  return m.uordblks + m.hblkhd;
}

#else
static double get_chunks(pcu::PCU *PCUObj)
{
  if(!PCUObj->Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}
#endif

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
  auto pcu_obj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  print_stats("malloc used before", get_chunks(pcu_obj.get()), pcu_obj.get());
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],pcu_obj.get());
  m->verify();
  print_stats("kernel heap", get_peak(pcu_obj.get()), pcu_obj.get());
  print_stats("malloc used", get_chunks(pcu_obj.get()), pcu_obj.get());
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
