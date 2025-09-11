#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <lionPrint.h>
#include <parma.h>

#ifdef PUMI_HAS_SIMMETRIX
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
  pcu::Init(&argc,&argv);
  {  
  pcu::PCU pcu_obj;
  lion_set_verbosity(1);
  if ( argc != 5 && argc != 6) {
    if ( !pcu_obj.Self() )
      printf("Usage: %s <model> <mesh> <number of output parts> <partition file prefix>\n", argv[0]);
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
  if (pcu_obj.Peers() > 1) {
    if ( !pcu_obj.Self() )
      printf("This tool must be run in serial.\n");
    pcu::Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef PUMI_HAS_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();


  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);

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
  snprintf(filename , 256, "%s_%d.ptn", argv[4], num_ranks);
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
#ifdef PUMI_HAS_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  }
  pcu::Finalize();
}
