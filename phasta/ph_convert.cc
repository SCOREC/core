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
#include <ph.h>
#include <phRestart.h>
#include <phInput.h>
#include <apfGeometry.h>
#include <cassert>
#include <cstdlib>
#include <iostream>

namespace {
  static FILE* openFileRead(ph::Input&, const char* path) {
    return fopen(path, "r");
  }

  static FILE* openFileWrite(ph::Output&, const char* path) {
    return fopen(path, "w");
  }
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
  ma::Input* in = ma::configureIdentity(m);
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}

static void fixCoords(apf::Mesh2* m)
{
  PCU_Comm_Begin();
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  apf::Vector3 x;
  apf::Vector3 p;
  while ((e = m->iterate(it)))
    if (m->isShared(e) && m->isOwned(e)) {
      apf::Copies remotes;
      m->getRemotes(e, remotes);
      m->getPoint(e, 0, x);
      m->getParam(e, p);
      APF_ITERATE(apf::Copies, remotes, rit) {
        PCU_COMM_PACK(rit->first, rit->second);
        PCU_COMM_PACK(rit->first, x);
        PCU_COMM_PACK(rit->first, p);
      }
    }
  m->end(it);
  PCU_Comm_Send();
  double max_x_diff = 0;
  double max_p_diff = 0;
  apf::Vector3 max_x_diff_point;
  apf::Vector3 max_p_diff_point;
  int x_diffs = 0;
  int p_diffs = 0;
  while (PCU_Comm_Receive()) {
    apf::Vector3 ox, op;
    PCU_COMM_UNPACK(e);
    PCU_COMM_UNPACK(ox);
    PCU_COMM_UNPACK(op);
    m->getPoint(e, 0, x);
    m->getParam(e, p);
    if (!(apf::areClose(p, op, 0.0))) {
      ++p_diffs;
      double p_diff = (op - p).getLength();
      if (p_diff > max_p_diff) {
        max_p_diff = p_diff;
        max_p_diff_point = x;
      }
      m->setParam(e, op);
    }
    if (!(apf::areClose(x, ox, 0.0))) {
      ++x_diffs;
      double x_diff = (ox - x).getLength();
      if (x_diff > max_x_diff) {
        max_x_diff = x_diff;
        max_x_diff_point = x;
      }
      m->setPoint(e, 0, ox);
    }
  }
  double global_max[2];
  global_max[0] = max_x_diff;
  global_max[1] = max_p_diff;
  PCU_Max_Doubles(global_max, 2);
  long global_diffs[2];
  global_diffs[0] = x_diffs;
  global_diffs[1] = p_diffs;
  PCU_Add_Longs(global_diffs, 2);
  /* admittedly not the best way of checking
     which processor had the max */
  if (global_diffs[0] && (global_max[0] == max_x_diff))
    fprintf(stderr, "%ld spatial mismatches corrected, max distance %e\n",
        global_diffs[0], global_max[0]);
  if (global_diffs[1] && (global_max[1] == max_p_diff))
    fprintf(stderr, "%ld parametric mismatches corrected, max distance %e\n",
        global_diffs[1], global_max[1]);
}

static void postConvert(apf::Mesh2* m)
{
  fixMatches(m);
  fixPyramids(m);
  fixCoords(m);
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
  gmi_sim_start();
  gmi_register_sim();
  PCU_Protect();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl = gmi_load(argv[1]);
  pGModel simModel = gmi_export_sim(mdl);
  pParMesh sim_mesh = PM_load(argv[2], sthreadNone, simModel, progress);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  ph::buildMapping(simApfMesh);
  
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  postConvert(mesh);
  mesh->writeNative(argv[3]);
  std::string restartPath = ph::setupOutputDir();
  ph::Input phIn;
  phIn.openfile_read = openFileRead;
  ph::Output phOut;
  phOut.openfile_write = openFileWrite;
  phIn.ensa_dof = 5;
  phIn.timeStepNumber = 0;
  phIn.displacementMigration = false;
  phIn.dwalMigration = false;
  phIn.buildMapping = true;
  ph::detachAndWriteSolution(phIn, phOut, mesh, restartPath);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}

