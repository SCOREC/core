#include <lionPrint.h>
#include <MeshSim.h>
#include <SimModel.h>
#include <SimDiscrete.h>
#include <SimAdvMeshing.h>
#include <SimPartitionedMesh.h>
#include <SimModelerAdvUtil.h>
#include <SimModelerUtil.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <apf_simConfig.h>
#include "gmi_sim_config.h"
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <ma.h>
#include <parma.h>
#include <PCU.h>
#include <pcu_util.h>
#include <cstdlib>
#include <memory>

#include <iostream> //cout
#include <getopt.h> //option parser

#ifdef SIM_PARASOLID
#include "SimParasolidKrnl.h"
#endif

#ifdef SIM_ACIS
#include "SimAcisKrnl.h"
#endif


#ifndef AdvMeshing_EXPORT
#define AdvMeshing_EXPORT
#endif

#include "SimAttribute.h"
#include "ModelTypes.h"


namespace {

int should_log = 0;
int disable_volume = 0;
int disable_surface = 0;
int should_fix_pyramids = 1;
int should_attach_order = 0;
std::string modelFile;
std::string nativeModelFile;
std::string surfaceMeshFile;
std::string caseName;
std::string outMeshFile;

void messageHandler(int type, const char* msg)
{
  switch (type) {
    case Sim_WarningMsg:
      if(!PCU_Comm_Self())
        fprintf(stdout, "Warning SimModeler %s\n", msg);
      break;
    case Sim_ErrorMsg:
      if(!PCU_Comm_Self()) 
        fprintf(stdout, "Error SimModeler %s ... exiting\n", msg);
      MPI_Finalize();
      exit(EXIT_SUCCESS); 
      break;
    default:
      // Ignore sim info and debug messages
      break;
  }
}

pParMesh generate(pGModel mdl, std::string meshCaseName, pcu::PCU *PCUObj) {
  pAManager attmngr = SModel_attManager(mdl);
#if SIMMODSUITE_MAJOR_VERSION <= 2024 && SIMMODSUITE_MINOR_VERSION < 240219
  pAManager attmngr = GM_attManager(mdl);
#else
  pAManager attmngr = GM_attManager(mdl,true);
  if(0==PCUObj->Self())
    fprintf(stdout, "Loading mesh case %s...\n", meshCaseName.c_str());
  pACase mcaseFile = AMAN_findCase(attmngr, meshCaseName.c_str());
  PCU_ALWAYS_ASSERT(mcaseFile);

  AttCase_setModel(mcaseFile, mdl);
  AttCase_associate(mcaseFile, NULL);
  pACase mcase = MS_newMeshCase(mdl);
  MeshingOptions meshingOptions;
  MS_processSimModelerMeshingAtts(mcaseFile, mcase, &meshingOptions);
  MS_processSimModelerAdvMeshingAtts(mcaseFile, mcase);
  AttCase_setModel(mcase, mdl);

  pParMesh pmesh;
  if( ! surfaceMeshFile.empty() &&
      disable_surface && !disable_volume ) {
    //load the surface mesh instead of creating it
    pmesh = PM_load(surfaceMeshFile.c_str(), mdl, NULL);
    PM_setTotalNumParts(pmesh, PMU_size()); //enable parallel volume meshing
  } else {
    //create an empty surface mesh
    pmesh = PM_new(0, mdl, PMU_size());
  }

  if( !disable_surface ) {
    const double stime = MPI_Wtime();
    if(0==PCUObj->Self()) {
      printf("Meshing surface..."); fflush(stdout);
    }
    pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, pmesh);
    SurfaceMesher_execute(surfaceMesher, NULL);
    SurfaceMesher_delete(surfaceMesher);
    if(0==PCUObj->Self())
      printf(" %f seconds\n", MPI_Wtime()-stime);
    if( ! surfaceMeshFile.empty() ) {
      if(0==PCUObj->Self())
        printf(" writing surface mesh %s\n", surfaceMeshFile.c_str());
      PM_write(pmesh, surfaceMeshFile.c_str(), NULL);
    }
  }

  if( !disable_volume ) {
    const double vtime = MPI_Wtime();
    if(0==PCUObj->Self()) {
      printf("Meshing volume..."); fflush(stdout);
    }
    pVolumeMesher volumeMesher = VolumeMesher_new(mcase, pmesh);
    VolumeMesher_execute(volumeMesher, NULL);
    VolumeMesher_delete(volumeMesher);
    if(0==PCUObj->Self())
      printf(" %f seconds\n", MPI_Wtime()-vtime);
  }

  return pmesh;
}

void getConfig(int argc, char** argv, pcu::PCU *PCUObj) {
  opterr = 0;

  static struct option long_opts[] = {
    {"no-pyramid-fix", no_argument, &should_fix_pyramids, 0},
    {"attach-order", no_argument, &should_attach_order, 1},
    {"enable-log", no_argument, &should_log, 1},
    {"disable-volume", no_argument, &disable_volume, 1},
    {"disable-surface", no_argument, &disable_surface, 1},
    {"native-model", required_argument, 0, 'n'},
    {"surface-mesh", required_argument, 0, 'm'},
    {0, 0, 0, 0}  // terminate the option array
  };

  const char* usage=""
    "[options] <GeomSim model (.smd)> <mesh case name>\n"
    "options:\n"
    "  --enable-log                            Enable Simmetrix logging\n"
    "  --disable-volume                        Disable volume mesh generation\n"
    "  --disable-surface                       Disable suface mesh generation\n"
    "  --no-pyramid-fix                        Disable quad-connected pyramid tetrahedronization\n"
    "  --attach-order                          Attach the Simmetrix element order as a Numbering\n"
    "  --native-model=/path/to/model           Load the native Parasolid or ACIS model that the GeomSim model uses\n"
    "  --surface-mesh=/path/to/surfaceMesh  read or write the surface mesh - depends on generation mode\n";

  nativeModelFile = "";
  surfaceMeshFile = "";
  int option_index = 0;
  while(1) {
    int c = getopt_long(argc, argv, "", long_opts, &option_index);
    if (c == -1) break; //end of options
    switch (c) {
      case 0: // enable-log|disable-volume|disable-surf
        if (!PCUObj->Self())
          printf ("read arg %d\n", c);
        break;
      case 'n':
        nativeModelFile = std::string(optarg);
        break;
      case 'm':
        surfaceMeshFile = std::string(optarg);
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


  if(argc-optind != 2) {
    if (!PCUObj->Self())
      printf("Usage %s %s", argv[0], usage);
    exit(EXIT_FAILURE);
  }
  int i=optind;
  modelFile = std::string(argv[i++]);
  outMeshFile = caseName = std::string(argv[i++]);
  outMeshFile.append("/");

  if (!PCUObj->Self()) {
    std::cout << "enable_log " << should_log <<
                 " disable_surface " << disable_surface <<
                 " disable_volume " << disable_volume <<
                 " fix_pyramids " << should_fix_pyramids <<
                 " attach_order " << should_attach_order <<
                 " native-model " << nativeModelFile <<
                 " model " << modelFile <<
                 " surface mesh " << surfaceMeshFile <<
                 " case name " << caseName <<
                 " output mesh" << outMeshFile << "\n";
  }
}

static void attachOrder(apf::Mesh* m)
{
  apf::numberOverlapDimension(m, "sim_order", m->getDimension());
}

void fixMatches(apf::Mesh2* m)
{
  if (m->hasMatching()) {
    if (apf::alignMdsMatches(m))
      printf("fixed misaligned matches\n");
    else
      printf("matches were aligned\n");
    PCU_ALWAYS_ASSERT( ! apf::alignMdsMatches(m));
  }
}

void fixPyramids(apf::Mesh2* m)
{
  ma::Input* in = ma::makeAdvanced(ma::configureIdentity(m));
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}

void postConvert(apf::Mesh2* m)
{
  fixMatches(m);
  if (should_fix_pyramids) fixPyramids(m);
  m->verify();
}

#if defined(SIM_ACIS) || defined(SIM_PARASOLID)
bool hasExtension(std::string s, std::string ext) {
  if(s.substr(s.find_last_of(".") + 1) == ext) {
    return true;
  } else {
    return false;
  }
}
#endif

pNativeModel loadNativeModel(pcu::PCU *PCUObj) {
  enum { TEXT_FORMAT = 0 };
  pNativeModel nm = 0;
  if ( nativeModelFile.empty() ) {
    nm = 0;
#ifdef SIM_ACIS
  } else if (hasExtension(nativeModelFile, "sat")) {
    nm = AcisNM_createFromFile(nativeModelFile.c_str(), TEXT_FORMAT);
    if(!PCUObj->Self())
      printf("loaded acis native model\n");
#endif
#ifdef SIM_PARASOLID
  } else if (hasExtension(nativeModelFile, "xmt_txt") ||
             hasExtension(nativeModelFile, "x_t")        ) {
    nm = ParasolidNM_createFromFile(nativeModelFile.c_str(), TEXT_FORMAT);
    if(!PCUObj->Self())
      printf("loaded parasolid native model\n");
#endif
  } else {
    if(!PCUObj->Self())
      printf("native model file has bad extension");
    exit(EXIT_FAILURE);
  }
  return nm;
}

void simStart() {
  SimPartitionedMesh_start(NULL,NULL);
  if(should_log)
    Sim_logOn("generate_sim.log");
  MS_init();
  SimModel_start();
#ifdef SIM_PARASOLID
  SimParasolid_start(1);
#endif
#ifdef SIM_ACIS
  SimAcis_start(1);
#endif
  SimDiscrete_start(0);
  Sim_readLicenseFile(NULL);
  MS_init();
  SimAdvMeshing_start();
  Sim_setMessageHandler(messageHandler);
}

void simStop() {
  SimAdvMeshing_stop();
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  SimDiscrete_stop(0);
#ifdef SIM_ACIS
  SimAcis_stop(1);
#endif
#ifdef SIM_PARASOLID
  SimParasolid_stop(1);
#endif
  MS_exit();
}

} //end unnamed namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));
  lion_set_verbosity(1);
  pcu::Protect();
  getConfig(argc,argv,PCUObj.get());

  if (should_attach_order && should_fix_pyramids) {
    if (!PCUObj.get()->Self())
      std::cout << "disabling pyramid fix because --attach-order was given\n";
    should_fix_pyramids = 0;
  }

  simStart();
  pNativeModel nm = loadNativeModel(PCUObj.get());
  pGModel simModel = GM_load(modelFile.c_str(), nm, NULL);

  const double t0 = MPI_Wtime();
  pParMesh sim_mesh = generate(simModel, caseName, PCUObj.get());
  const double t1 = MPI_Wtime();
  if(!PCUObj.get()->Self())
    printf("Mesh generated in %f seconds\n", t1-t0);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh, PCUObj.get());
  if (should_attach_order) attachOrder(simApfMesh);

  gmi_register_sim();
  gmi_model* model = gmi_import_sim(simModel);
  apf::Mesh2* mesh = apf::createMdsMesh(model, simApfMesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  postConvert(mesh);
  Parma_PrintPtnStats(mesh, "");
  mesh->writeNative(outMeshFile.c_str());

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  simStop();
  Sim_logOff();
  }
  MPI_Finalize();
}
