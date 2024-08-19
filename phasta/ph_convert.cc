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
#include <ph.h>
#include <phRestart.h>
#include <phInput.h>
#include <apfGeometry.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <getopt.h>


static void attachOrder(apf::Mesh* m)
{
  apf::numberOverlapDimension(m, "sim_order", m->getDimension());
}
namespace {
  static FILE* openFileRead(ph::Input&, const char* path, pcu::PCU*) {
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
      lion_oprint(1,"fixed misaligned matches\n");
    else
      lion_oprint(1,"matches were aligned\n");
    PCU_ALWAYS_ASSERT( ! apf::alignMdsMatches(m));
  }
}

static void fixPyramids(apf::Mesh2* m)
{
  if (m->getDimension() != 3)
    return; /* no pyramids exist in 2D */
  if (apf::countEntitiesOfType(m, apf::Mesh::HEX))
    return; /* meshadapt can't even look at hexes */
  ma::Input* in = ma::makeAdvanced(ma::configureIdentity(m));
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}
const char* gmi_path = NULL;
const char* gmi_native_path = NULL;
const char* sms_path = NULL;
const char* smb_path = NULL;
int should_log = 0;
int should_fix_pyramids = 1;
int should_attach_order = 0;
bool found_bad_arg = false;

void getConfig(int argc, char** argv, pcu::PCU *pcu_obj) {

  opterr = 0;

  static struct option long_opts[] = {
    {"no-pyramid-fix", no_argument, &should_fix_pyramids, 0},
    {"attach-order", no_argument, &should_attach_order, 1},
    {"enable-log", no_argument, &should_log, 2},
    {"native-model", required_argument, 0, 'n'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "[options] <model file> <simmetrix mesh> <scorec mesh>\n"
    "Convert a Simmetrix mesh classified on a GeomSim, and/or a Parasolid/ACIS native model, to a PUMI mesh.\n"
    "During conversion, the ordering of mesh entities in the Simmetrix mesh is attached to the PUMI mesh.\n"
    "This information primarily supports transferring size fields from the PUMI mesh to the Simmetrix mesh to\n"
    "drive Simmetrix mesh adaptation.\n\n"
    "<options> <default value> <brief description>:\n"
    "  --no-pyramid-fix                on   Disable quad-connected pyramid tetrahedronization\n"
    "  --attach-order                  off   Attach the Simmetrix element order as a Numbering\n"
    "  --enable-log                    off  Enable Simmetrix logging\n"
    "  --native-model=/path/to/model        Load the native Parasolid or ACIS model that the GeomSim model uses\n";

  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 0: // pyramid fix flag
      case 1: // attach order flag
      case 2: // enable simmetrix logging
        break;
      case 'n':
        gmi_native_path = optarg;
        break;
      case '?':
        if (!pcu_obj->Self())
          lion_oprint(1,"warning: skipping unrecognized option\n");
        break;
      default:
        if (!pcu_obj->Self())
          lion_oprint(1,"Usage %s %s", argv[0], usage);
        exit(EXIT_FAILURE);
    }
  }

  if(argc-optind != 3) {
    if (!pcu_obj->Self())
      lion_oprint(1,"Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  gmi_path = argv[i++];
  sms_path = argv[i++];
  smb_path = argv[i++];

  if (!pcu_obj->Self()) {
    lion_oprint(1,"fix_pyramids %d attach_order %d enable_log %d\n",
            should_fix_pyramids, should_attach_order, should_log);
    lion_oprint(1,"native-model \'%s\' model \'%s\' simmetrix mesh \'%s\' output mesh \'%s\'\n",
      gmi_native_path, gmi_path, sms_path, smb_path);
  }
}


static void fixCoords(apf::Mesh2* m)
{
  m->getPCU()->Begin();
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
        m->getPCU()->Pack(rit->first, rit->second);
        m->getPCU()->Pack(rit->first, x);
        m->getPCU()->Pack(rit->first, p);
      }
    }
  m->end(it);
  m->getPCU()->Send();
  double max_x_diff = 0;
  double max_p_diff = 0;
  apf::Vector3 max_x_diff_point;
  apf::Vector3 max_p_diff_point;
  int x_diffs = 0;
  int p_diffs = 0;
  while (m->getPCU()->Receive()) {
    apf::Vector3 ox, op;
    m->getPCU()->Unpack(e);
    m->getPCU()->Unpack(ox);
    m->getPCU()->Unpack(op);
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
  m->getPCU()->Max<double>(global_max, 2);
  long global_diffs[2];
  global_diffs[0] = x_diffs;
  global_diffs[1] = p_diffs;
  m->getPCU()->Add<long>(global_diffs, 2);
  /* admittedly not the best way of checking
     which processor had the max */
  if (global_diffs[0] && (global_max[0] == max_x_diff))
    lion_eprint(1, "%ld spatial mismatches corrected, max distance %e\n",
        global_diffs[0], global_max[0]);
  if (global_diffs[1] && (global_max[1] == max_p_diff))
    lion_eprint(1, "%ld parametric mismatches corrected, max distance %e\n",
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
  {
  auto pcu_obj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  getConfig(argc, argv, pcu_obj.get());
  if( should_log )
    Sim_logOn("convert.sim.log");

  if (should_attach_order && should_fix_pyramids) {
    if (!pcu_obj.get()->Self())
      std::cout << "disabling pyramid fix because --attach-order was given\n";
    should_fix_pyramids = false;
  }

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  gmi_model* mdl;
  if( gmi_native_path )
    mdl = gmi_sim_load(gmi_native_path,gmi_path);
  else
    mdl = gmi_load(gmi_path);
  pGModel simModel = gmi_export_sim(mdl);
  double t0 = pcu::Time();
  pParMesh sim_mesh = PM_load(sms_path, simModel, progress);
  double t1 = pcu::Time();
  if(!pcu_obj.get()->Self())
    lion_eprint(1, "read and created the simmetrix mesh in %f seconds\n", t1-t0);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh, pcu_obj.get());
  double t2 = pcu::Time();
  if(!simApfMesh->getPCU()->Self())
    lion_eprint(1, "created the apf_sim mesh in %f seconds\n", t2-t1);
  if (should_attach_order) attachOrder(simApfMesh);
  ph::buildMapping(simApfMesh);

  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  double t3 = pcu::Time();
  if(!mesh->getPCU()->Self())
    lion_eprint(1, "created the apf_mds mesh in %f seconds\n", t3-t2);

  apf::printStats(mesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  postConvert(mesh);
  mesh->writeNative(smb_path);
  std::string restartPath = ph::setupOutputDir(mesh->getPCU());
  ph::Input phIn;
  phIn.openfile_read = openFileRead;
  ph::Output phOut;
  phOut.openfile_write = openFileWrite;
  phIn.ensa_dof = 5;
  phIn.timeStepNumber = 0;
  phIn.displacementMigration = false;
  phIn.dwalMigration = false;
  phIn.buildMapping = true;
  ph::BCs bcs;
  ph::generateOutput(phIn, bcs, mesh, phOut);
  ph::detachAndWriteSolution(phIn, phOut, mesh, restartPath);

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
  if( should_log )
    Sim_logOff();
  }
  MPI_Finalize();
}
