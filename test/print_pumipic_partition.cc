#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <lionPrint.h>
#include <parma.h>

#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <apfZoltan.h>
#include <pcu_util.h>
#include <cstdlib>
#include <fstream>

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  if ( argc != 5 && argc != 6) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <number of ranks> <partition file prefix>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  if (PCU_Comm_Peers() > 1) {
    if ( !PCU_Comm_Self() )
      printf("This tool must be run in serial.\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();


  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  int num_ranks = atoi(argv[3]);
  //Partition the mesh (Taken from zsplit.cc)
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, num_ranks);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;

  //Write the partition out
  char filename[256];
  sprintf(filename , "%s_%d.ptn", argv[4], num_ranks);
  printf("Writing partition to %s\n", filename);
  std::ofstream out(filename);
  apf::MeshIterator* mitr = m->begin(m->getDimension());
  while (apf::MeshEntity* ent = m->iterate(mitr)) {
    if (plan->has(ent))
      out << plan->sending(ent) << '\n';
    else
      out << 0 << '\n';
  }
  m->end(mitr);

  delete plan;

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
