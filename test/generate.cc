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

#ifndef AdvMeshing_EXPORT
#define AdvMeshing_EXPORT
#endif

#include "SimAttribute.h"
#include "ModelTypes.h"

/* The fields in this sturcture are set mainly based on the "Surface Meshing" and "Volume Meshing"
   attributes in the meshing case. */
struct AdvMeshing_EXPORT MeshingOptions {
  bool surfaceRun;                // whether to run surface meshing
  bool surfaceDoFixIntersections; // whether to run fix self intersections
  bool surfaceDoCurve;            // whether to curve surface mesh

  // These correspond as indicated to the fields in the "Surface Meshing" attribute
  int surfaceSmoothingLevel;      // "Smoothing Level"
  int surfaceSmoothingType;       // "Smoothing Type"
  int surfaceFixIntersections;    // "Fix Intersections"
  int surfaceOptimization;        // "Optimization"
  int surfaceEnforceGradation;    // "Enforce Spatial Gradation"
  double surfaceFaceRotationLimit;// "Discrete Face Rotation Limit"
  int surfaceCurveType;           // ignore

  bool volumeRun;                 // whether to run volume meshing
  bool volumeDoStructured;        // always true
  bool volumeDoUnstructured;      // always true
  bool volumeDoCurve;             // whether to run mesh curving

  // These correspond as indicated to the fields in the "Volume Meshing" attribute
  int volumeEnforceSize;          // "Enforce Size"
  int volumeSmoothingLevel;       // "Smoothing Level"
  int volumeSmoothingType;        // "Smoothing Type"
  int volumeOptimization;         // "Optimization"
  int volumeCurveType;            // ignore
};

AdvMeshing_EXPORT void MS_setupSimModelerMeshCase(pACase sourceCase, pACase meshCase, MeshingOptions *options);
pAManager SModel_attManager(pModel model);

pParMesh generate(pGModel mdl, std::string caseName) {
  pAManager attmngr = SModel_attManager(mdl);
  pACase mcaseFile = AMAN_findCase(attmngr, caseName.c_str());
  assert(mcaseFile);

  AttCase_setModel(mcaseFile, mdl);
  pACase mcase = MS_newMeshCase(mdl);
  MeshingOptions meshingOptions;
  MS_setupSimModelerMeshCase(mcaseFile, mcase, &meshingOptions);
  AttCase_setModel(mcase, mdl);

  pParMesh pmesh = PM_new(0, mdl, PMU_size());

  pProgress prog = Progress_new();
  pSurfaceMesher surfaceMesher = SurfaceMesher_new(mcase, pmesh);
  SurfaceMesher_execute(surfaceMesher, prog);

  pVolumeMesher volumeMesher = VolumeMesher_new(mcase, pmesh);
  VolumeMesher_execute(volumeMesher, prog);

  SurfaceMesher_delete(surfaceMesher);
  VolumeMesher_delete(volumeMesher);
  return pmesh;
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  if (argc != 4) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0] << " <model file> <simmetrix mesh> <pumi mesh>\n";
    return 0;
  }
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);
  SimModel_start();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  pGModel simModel;
  simModel = GM_load(argv[1], NULL, progress);

  pParMesh sim_mesh = generate(simModel, argv[2]);
  apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);
  
  gmi_register_sim();
  gmi_model* pumiMdl = gmi_import_sim(simModel);
  apf::Mesh2* pumiApfMesh = apf::createMdsMesh(pumiMdl, simApfMesh);
  apf::destroyMesh(simApfMesh);
  M_release(sim_mesh);
  if (alignMdsMatches(pumiApfMesh))
    printf("fixed misaligned matches\n");
  else
    printf("matches (if any) are aligned ok\n");
  pumiApfMesh->verify();
  pumiApfMesh->writeNative(argv[3]);

  pumiApfMesh->destroyNative();
  apf::destroyMesh(pumiApfMesh);

  Progress_delete(progress);
  SimModel_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
}
