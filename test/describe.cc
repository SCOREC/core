#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <parma.h>
#include <SimUtil.h>
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif
#include <cassert>

#ifdef __bgq__
#include <spi/include/kernel/memory.h>

static double get_peak()
{
  uint64_t heap;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  return heap;
}

#elif defined (__linux__)

static double get_peak()
{
  return mallinfo().arena;
}

#else

static double get_peak()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}

#endif

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

#if defined(__linux__)

static double get_chunks()
{
  struct mallinfo m = mallinfo();
  return m.uordblks + m.hblkhd;
}

#else
static double get_chunks()
{
  if(!PCU_Comm_Self())
    printf("%s:%d: OS Not supported\n", __FILE__, __LINE__);
  return(-1.0);
}
#endif

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
  assert(argc==3);
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_mesh();
  gmi_register_sim();
  print_stats("malloc used before", get_chunks());
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  print_stats("kernel heap", get_peak());
  print_stats("malloc used", get_chunks());
  Parma_PrintPtnStats(m, "");
  list_tags(m);
  m->destroyNative();
  apf::destroyMesh(m);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
