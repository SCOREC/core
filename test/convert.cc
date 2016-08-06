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
#include <apfNumbering.h>
#include <ma.h>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iostream>

static void attachOrder(apf::Mesh* m)
{
  apf::numberOverlapDimension(m, "sim_order", m->getDimension());
}

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

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  const char* gmi_path = NULL;
  const char* sms_path = NULL;
  const char* smb_path = NULL;
  bool should_fix_pyramids = true;
  bool should_attach_order = false;
  bool found_bad_arg = false;

  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "--no-pyramid-fix")) {
      should_fix_pyramids = false;
    } else if (!strcmp(argv[i], "--attach-order")) {
      should_attach_order = true;
    } else if (!gmi_path) {
      gmi_path = argv[i];
    } else if (!sms_path) {
      sms_path = argv[i];
    } else if (!smb_path) {
      smb_path = argv[i];
    } else {
      if(!PCU_Comm_Self())
        std::cerr << "bad argument \"" << argv[i] << "\"\n";
      found_bad_arg = true;
    }
  }

  if (!gmi_path || !sms_path || !smb_path || found_bad_arg) {
    if(!PCU_Comm_Self()) {
      std::cout << "usage: " << argv[0] << " [options] <model file> <simmetrix mesh> <scorec mesh>\n";
      std::cout << "options:\n";
      std::cout << "  --no-pyramid-fix           Disable quad-connected pyramid tetrahedronization\n";
      std::cout << "  --attach-order             Attach the Simmetrix element order as a Numbering\n";
    }
    return EXIT_FAILURE;
  }

  if (should_attach_order && should_fix_pyramids) {
    if (!PCU_Comm_Self())
      std::cout << "disabling pyramid fix because --attach-order was given\n";
    should_fix_pyramids = false;
  }

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl = gmi_load(gmi_path);
  pGModel simModel = gmi_export_sim(mdl);
  pParMesh sim_mesh = PM_load(sms_path, sthreadNone, simModel, progress);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  if (should_attach_order) attachOrder(simApfMesh);
  
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  fixMatches(mesh);
  if (should_fix_pyramids) fixPyramids(mesh);
  mesh->verify();
  mesh->writeNative(smb_path);

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
