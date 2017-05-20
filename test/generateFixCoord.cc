#include <PCU.h>
#include <MeshSim.h>
#include <SimAdvMeshing.h>
#include <SimPartitionedMesh.h>
#include <SimModelerAdvUtil.h>
#include <SimModelerUtil.h>
#include <SimUtil.h>
#include <apfSIM.h>
#include <apfMDS.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <ma.h>
#include <parma.h>
#include <pcu_util.h>
#include <cstdlib>

#ifndef AdvMeshing_EXPORT
#define AdvMeshing_EXPORT
#endif

#include "SimAttribute.h"
#include "ModelTypes.h"

pAManager SModel_attManager(pModel model);

namespace {

std::string modelFile;
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


static void FixTensorProductCoordinates(apf::Mesh* mesh)
{
  apf::Mesh* m = mesh;
  int n = m->count(0);
  apf::MeshEntity* v;
  int i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((v = m->iterate(it))) {
    apf::Vector3 p;
    m->getPoint(v, 0, p);
    p[0]=p[0]*1.353701006e-3;
    p[2]=p[2]*5.454804022e-4;
    p[1]=p[1]*8.18e-4;
    m->setPoint(v, 0, p);
//    for (int j = 0; j < 3; ++j)
//      x[j * n + i] = p[j]; /* FORTRAN indexing */
//    ++i;
  }
  m->end(it);
  PCU_ALWAYS_ASSERT(i == n);
 // o.arrays.coordinates = x;
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
  if (argc != 3) {
    if(0==PCU_Comm_Self())
      printf("Usage: %s <model (.smd)> <mesh case name>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_SUCCESS);
  }
  modelFile = argv[1];
  outMeshFile = caseName = argv[2];
  outMeshFile.append("/");
  if(0==PCU_Comm_Self())
    printf("Inputs: model %s case %s\n", modelFile.c_str(), caseName.c_str());
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

} //end unnamed namespace

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  getConfig(argc,argv);

  SimModel_start();
  SimPartitionedMesh_start(NULL,NULL);
  SimUtil_start();
  Sim_readLicenseFile(NULL);
  MS_init();
  SimAdvMeshing_start(); 
  Sim_setMessageHandler(messageHandler);

  pGModel simModel = GM_load(modelFile.c_str(), NULL, NULL);

  const double t0 = MPI_Wtime();
  pParMesh sim_mesh = generate(simModel, caseName);
  const double t1 = MPI_Wtime();
  if(0==PCU_Comm_Self())
    printf("Mesh generated in %f seconds\n", t1-t0);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  
  gmi_register_sim();
  gmi_model* model = gmi_import_sim(simModel);
  apf::Mesh2* mesh = apf::createMdsMesh(model, simApfMesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  mesh->verify();
  //postConvert(mesh);
  FixTensorProductCoordinates(mesh);
  Parma_PrintPtnStats(mesh, "");
  mesh->writeNative(outMeshFile.c_str());

  mesh->destroyNative();
  apf::destroyMesh(mesh);

  SimAdvMeshing_stop(); 
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}
