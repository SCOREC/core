#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <ma.h>

static void fixMatches(apf::Mesh2* m)
{
  if (m->hasMatching()) {
    if (apf::alignMdsMatches(m))
      printf("fixed misaligned matches\n");
    else
      printf("matches were aligned\n");
    assert( ! apf::alignMdsMatches(m));
  }
}

static void fixPyramids(apf::Mesh2* m)
{
  if (m->getDimension() != 3)
    return; /* no pyramids exist in 2D */
  if (apf::countEntitiesOfType(m, apf::Mesh::HEX))
    return; /* meshadapt can't even look at hexes */
  ma::Input* in = ma::configureIdentity(m);
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}

static void postConvert(apf::Mesh2* m)
{
  fixMatches(m);
  fixPyramids(m);
  m->verify();
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0] << " <model file> <simmetrix mesh> <scorec mesh>\n";
    return EXIT_FAILURE;
  }
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);
  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl = gmi_load(argv[1]);
  pGModel simModel = gmi_export_sim(mdl);
  pParMesh sim_mesh = PM_load(argv[2], sthreadNone, simModel, progress);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  postConvert(mesh);
  mesh->writeNative(argv[3]);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
