#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <PCU.h>
#include <parma.h>
#if defined(__linux__)
#include <malloc.h>
#else
#include <cstdlib>
#endif

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
  min = value;
  PCU_Min_Doubles(&min, 1);
  max = value;
  PCU_Max_Doubles(&max, 1);
  avg = value;
  PCU_Add_Doubles(&avg, 1);
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
  gmi_register_mesh();
  print_stats("malloc used before", get_chunks());
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);
  print_stats("kernel heap", get_peak());
  print_stats("malloc used", get_chunks());
  print_stats("elements", m->count(m->getDimension()));
  print_stats("vertices", m->count(0));
  Parma_PrintPtnStats(m, "");
  list_tags(m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}




