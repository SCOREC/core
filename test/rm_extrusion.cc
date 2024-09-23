#include <lionPrint.h>
#include <MeshSim.h>
#include <SimPartitionedMesh.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apf_simConfig.h>
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
#include <SimAdvMeshing.h>
#include <getopt.h>



const char* gmi_path = NULL;
const char* gmi_native_path = NULL;
const char* sms_path = NULL;
const char* smsNew_path = NULL;
int should_log = 0;
bool found_bad_arg = false;

#if SIMMODSUITE_MAJOR_VERSION >= 14
void M_removeSurfaceExtrusionConstraints(pUnstructuredMesh, pPList);
#else
void M_removeSurfaceExtrusionConstraints(pMesh, pPList);
#endif


void getConfig(int argc, char** argv, pcu::PCU *PCUObj) {

  opterr = 0;

  static struct option long_opts[] = {
    {"enable-log", no_argument, &should_log, 2},
    {"native-model", required_argument, 0, 'n'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "[options] <model file> <simmetrix mesh> <modified simmetrix mesh>\n"
    "options:\n"
    "  --enable-log                    Enable Simmetrix logging\n"
    "  --native-model=/path/to/model   Load the native Parasolid or ACIS model that the GeomSim model uses\n";

  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 2: // enable simmetrix logging
        break;
      case 'n':
        gmi_native_path = optarg;
        break;
      case '?':
        if (!PCUObj->Self())
          printf ("warning: skipping unrecognized option\n");
        break;
      default:
        if (!PCUObj->Self())
          printf("Usage %s %s", argv[0], usage);
        exit(EXIT_FAILURE);
    }
  }

  if(argc-optind != 3) {
    if (!PCUObj->Self())
      printf("Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  gmi_path = argv[i++];
  sms_path = argv[i++];
  smsNew_path = argv[i++];

  if (!PCUObj->Self()) {
    printf ("enable_log %d\n", should_log);
    printf ("native-model \'%s\' model \'%s\' simmetrix mesh \'%s\' output mesh \'%s\'\n",
      gmi_native_path, gmi_path, sms_path, smsNew_path);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  MS_init();
  SimAdvMeshing_start();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  getConfig(argc, argv, &PCUObj);
  if( should_log )
    Sim_logOn("rm_extrusion.sim.log");

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
  if(!PCUObj.Self())
    fprintf(stderr, "Read model\n");

  double t0 = pcu::Time();
  pMesh sim_mesh = M_load(sms_path, simModel, progress);
  if(!PCUObj.Self())
    fprintf(stderr, "Read mesh\n");

  M_removeSurfaceExtrusionConstraints(sim_mesh, NULL);
  if(!PCUObj.Self())
    fprintf(stderr, "Removed surface extrusion constraints\n");
  M_write(sim_mesh, smsNew_path, 0, progress);
  double t1 = pcu::Time();

  if(!PCUObj.Self())
    fprintf(stderr, "read the mesh, removed the face extrusion attributes, and wrote the mesh %f seconds\n", t1-t0);

  M_release(sim_mesh);
  Progress_delete(progress);
  gmi_destroy(mdl);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  SimAdvMeshing_stop();
  MS_exit();
  if( should_log )
    Sim_logOff();
  }
  MPI_Finalize();
}
