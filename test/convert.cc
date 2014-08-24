#include <PCU.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <SimParasolidKrnl.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>

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
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);
  SimModel_start();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel simModel;
  simModel = GM_load(argv[1], NULL, progress);

  pParMesh sim_mesh = PM_load(argv[2],sthreadNone,0,progress);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  
  gmi_register_sim();
  gmi_model* mdl = gmi_import_sim(simModel);
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  postConvert(mesh);
  mesh->writeNative(argv[3]);

  mesh->destroyNative();
  apf::destroyMesh(mesh);
  GM_release(simModel);

  Progress_delete(progress);
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}
