#include <apf.h>
#include <apfMDS.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <pcu_util.h>
#include <apfDynamicVector.h>
#include <apfDynamicMatrix.h>
#include <crv.h>
#include <cassert>
#include <stdlib.h>
#include <sstream>
#include <fstream>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */


static void safe_mkdir(
    const char* path);
static apf::Mesh2* makeEntMesh(
    apf::Mesh2* m,
    apf::MeshEntity* e);
static void makeCavityMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* cavityMeshLinear,
    apf::Mesh2* cavityMeshCurved);

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (PCU_Comm_Peers() > 1) {
    printf("%s should only be used for serial (single part) meshes!\n", argv[0]);
    printf("use the serialize utility to get a serial mesh, and retry!\n");
  }
  if (argc < 5) {
    if (PCU_Comm_Self() == 0) {
      printf("USAGE: %s <model> <mesh> <enttype> <prefix>\n", argv[0]);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  const char* modelFile = argv[1];
  const char* meshFile = argv[2];
  apf::Mesh::Type enttype = atoi(argv[3]);
  const char* prefix = argv[4];

  // load the mesh and get the order
  apf::Mesh2* m = apf::loadMdsMesh(modelFile,meshFile);
  int order = m->getShape()->getOrder();

  int dim = m->getDimension(); // mesh dimension
  int d = apf::Mesh::typeDimension[enttype]; // the dim we iterate on to create cavities

  apf::MeshEntity* e;
  apf::MeshIterator it = m->begin(d);
  int index = 0;

  // for now cavities are defined as follows
  // all the upward adjacent entities of dimension "dim" that
  // are adjacent to the verts of the "e". E.g., in case of edges,
  // this would give us the bi-directional edge collapse cavity.
  while ( (e = m->iterate(it)) ) {
    int etype = m->getType(e);
    int mtype = m->getModelType(m->toModel(e));
    apf::Mesh2* entMesh = makeEntMesh(m, e);
    apf::Mesh2* cavityMeshCurved = 0;
    apf::Mesh2* cavityMeshLinear = 0;
    makeCavityMeshes(m, e, cavityMeshLinear, cavityMeshCurved);

    // write the curved cavity mesh in native format for future retrieval

    // write the entMesh and cavityMeshLinear and cavityMeshCurved in vtk format

    // destroy the entMesh and cavityMeshes
    destroyMesh(entMesh);
    destroyMesh(cavityMeshLinear);
    destroyMesh(cavityMeshCurved);
    index++;
  }
  m->end(it);


  safe_mkdir(inPrefix);

  // rest of the clean up
  m->destroyNative();
  apf::destroyMesh(m);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif

  PCU_Comm_Free();
  MPI_Finalize();
}

static void safe_mkdir(
    const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("Err: could not create directory \"%s\"\n", path);
  }
}

static apf::Mesh2* makeEntMesh(
    apf::Mesh2* m,
    apf::MeshEntity* e)
{
}

static void makeCavityMeshes(
    apf::Mesh2* m,
    apf::MeshEntity* e,
    apf::Mesh2* cavityMeshLinear,
    apf::Mesh2* cavityMeshCurved)
{
  PCU_ALWAYS_ASSERT(!cavityMeshLinear);
  PCU_ALWAYS_ASSERT(!cavityMeshCurved);

  int dim = m->getDimension();

  cavityMeshLinear = apf::makeEmptyMdsMesh(0, dim, false);
  cavityMeshCurved = apf::makeEmptyMdsMesh(0, dim, false);
}
