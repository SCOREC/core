#include <PCU.h>
#include <MeshSim.h>
#include <SimAdvMeshing.h>
#include <SimPartitionedMesh.h>
#include <SimModelerAdvUtil.h>
#include <SimModelerUtil.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include "gmi_sim_config.h"
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <ma.h>
#include <parma.h>
#include <pcu_util.h>
#include <cstdlib>

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

pAManager SModel_attManager(pModel model);

namespace {

std::string modelFile;
std::string nativeModelFile;
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

pParMesh generate(pGModel mdl, std::string meshCaseName) {
  pAManager attmngr = SModel_attManager(mdl);
  if(0==PCU_Comm_Self())
    fprintf(stdout, "Loading mesh case %s...\n", meshCaseName.c_str());
  pACase mcaseFile = AMAN_findCase(attmngr, meshCaseName.c_str());
  PCU_ALWAYS_ASSERT(mcaseFile);

  AttCase_setModel(mcaseFile, mdl);
  pACase mcase = MS_newMeshCase(mdl);
  MeshingOptions meshingOptions;
  MS_processSimModelerMeshingAtts(mcaseFile, mcase, &meshingOptions);
  MS_processSimModelerAdvMeshingAtts(mcaseFile, mcase);
  AttCase_setModel(mcase, mdl);

  pParMesh pmesh = PM_new(0, mdl, PMU_size());

  const double stime = MPI_Wtime();
  if(0==PCU_Comm_Self()) {
    printf("Meshing surface..."); fflush(stdout);
  }
  pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, pmesh);
  SurfaceMesher_execute(surfaceMesher, NULL);
  SurfaceMesher_delete(surfaceMesher);
  if(0==PCU_Comm_Self())
    printf(" %f seconds\n", MPI_Wtime()-stime);

  const double vtime = MPI_Wtime();
  if(0==PCU_Comm_Self()) {
    printf("Meshing volume..."); fflush(stdout);
  }
  pVolumeMesher volumeMesher = VolumeMesher_new(mcase, pmesh);
  VolumeMesher_execute(volumeMesher, NULL);
  VolumeMesher_delete(volumeMesher);
  if(0==PCU_Comm_Self())
    printf(" %f seconds\n", MPI_Wtime()-vtime);

  return pmesh;
}

void getConfig(int argc, char** argv)
{
  if (argc < 3) {
    if(0==PCU_Comm_Self()) {
      printf("Usage: %s <GeomSim model (.smd)> <mesh case name> ", argv[0]);
      printf("       to generate a mesh on a GeomSim model\n");
      printf("   or: %s <SimModeler model (.smd)> <parasolid or acis native model> <mesh case name>\n", argv[0]);
      printf("       to generate a mesh using the specified case name using the SimModeler"
                     "model which references the native parasolid or acis model\n");
    }
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }
  modelFile = argv[1];
  if (argc == 3) {
    nativeModelFile = "";
    outMeshFile = caseName = argv[2];
  } else if (argc == 4) {
    nativeModelFile = argv[2];
    outMeshFile = caseName = argv[3];
  }
  outMeshFile.append("/");
  if(0==PCU_Comm_Self()) {
    printf("Inputs: model \'%s\' native model \'%s\' case \'%s\'\n",
        modelFile.c_str(), nativeModelFile.c_str(), caseName.c_str());
  }
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
  ma::Input* in = ma::configureIdentity(m);
  in->shouldCleanupLayer = true;
  ma::adapt(in);
}

void postConvert(apf::Mesh2* m)
{
  fixMatches(m);
  fixPyramids(m);
  m->verify();
}

bool hasExtension(std::string s, std::string ext) {
  if(s.substr(s.find_last_of(".") + 1) == ext) {
    return true;
  } else {
    return false;
  }
}

pNativeModel loadNativeModel() {
  enum { TEXT_FORMAT = 0 };
  pNativeModel nm = 0;
  if ( nativeModelFile.empty() ) {
    nm = 0;
#ifdef SIM_ACIS
  } else if (hasExtension(nativeModelFile, "sat")) {
    nm = AcisNM_createFromFile(nativeModelFile.c_str(), TEXT_FORMAT);
    if(!PCU_Comm_Self())
      printf("loaded acis native model\n");
#endif
#ifdef SIM_PARASOLID
  } else if (hasExtension(nativeModelFile, "xmt_txt") ||
             hasExtension(nativeModelFile, "x_t")        ) {
    nm = ParasolidNM_createFromFile(nativeModelFile.c_str(), TEXT_FORMAT);
    if(!PCU_Comm_Self())
      printf("loaded parasolid native model\n");
#endif
  } else {
    if(!PCU_Comm_Self())
      printf("native model file has bad extension");
    exit(EXIT_FAILURE);
  }
  return nm;
}

void simStart() {
  SimModel_start();
  SimPartitionedMesh_start(NULL,NULL);
  MS_init();
  SimModel_start();
#ifdef SIM_PARASOLID
  SimParasolid_start(1);
#endif
#ifdef SIM_ACIS
  SimAcis_start(1);
#endif
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
  PCU_Comm_Init();
  PCU_Protect();
  getConfig(argc,argv);

  simStart();
  pNativeModel nm = loadNativeModel();
  pGModel simModel = GM_load(modelFile.c_str(), nm, NULL);

  const double t0 = MPI_Wtime();
  pParMesh sim_mesh = generate(simModel, caseName);
  const double t1 = MPI_Wtime();
  if(!PCU_Comm_Self())
    printf("Mesh generated in %f seconds\n", t1-t0);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);

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
  PCU_Comm_Free();
  MPI_Finalize();
}
