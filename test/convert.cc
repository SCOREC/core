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
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include <getopt.h>


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
    PCU_ALWAYS_ASSERT( ! apf::alignMdsMatches(m));
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

const char* gmi_path = NULL;
const char* gmi_native_path = NULL;
const char* sms_path = NULL;
const char* smb_path = NULL;
int should_fix_pyramids = 1;
int should_attach_order = 0;
bool found_bad_arg = false;

void getConfig(int argc, char** argv) {

  opterr = 0;

  static struct option long_opts[] = {
    {"no-pyramid-fix", no_argument, &should_fix_pyramids, 0},
    {"attach-order", no_argument, &should_attach_order, 1},
    {"native-model", required_argument, 0, 'n'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "[options] <model file> <simmetrix mesh> <scorec mesh>\n"
    "options:\n"
    "  --no-pyramid-fix                Disable quad-connected pyramid tetrahedronization\n"
    "  --attach-order                  Attach the Simmetrix element order as a Numbering\n"
    "  --native-model=/path/to/model   Load the native Parasolid or ACIS model that the GeomSim model uses\n";

  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 0: // pyramid fix flag
      case 1: // attach order flag
        break;
      case 'n':
        gmi_native_path = optarg;
        break;
      case '?':
        if (!PCU_Comm_Self())
          printf ("warning: skipping unrecognized option\n");
        break;
      default:
        if (!PCU_Comm_Self())
          printf("Usage %s %s", argv[0], usage);
        exit(EXIT_FAILURE);
    }
  }

  if(argc-optind != 3) {
    if (!PCU_Comm_Self())
      printf("Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  gmi_path = argv[i++];
  sms_path = argv[i++];
  smb_path = argv[i++];

  if (!PCU_Comm_Self()) {
    printf ("fix_pyramids %d attach_order %d\n", should_fix_pyramids, should_attach_order);
    printf ("native-model \'%s\' model \'%s\' simmetrix mesh \'%s\' output mesh \'%s\'\n",
      gmi_native_path, gmi_path, sms_path, smb_path);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  getConfig(argc, argv);

  if (should_attach_order && should_fix_pyramids) {
    if (!PCU_Comm_Self())
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
  double t0 = PCU_Time();
  pParMesh sim_mesh = PM_load(sms_path, simModel, progress);
  double t1 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "read and created the simmetrix mesh in %f seconds\n", t1-t0);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  double t2 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_sim mesh in %f seconds\n", t2-t1);
  if (should_attach_order) attachOrder(simApfMesh);
  
  apf::Mesh2* mesh = apf::createMdsMesh(mdl, simApfMesh);
  double t3 = PCU_Time();
  if(!PCU_Comm_Self())
    fprintf(stderr, "created the apf_mds mesh in %f seconds\n", t3-t2);

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
  SimModel_stop();
  MS_exit();
  PCU_Comm_Free();
  MPI_Finalize();
}
