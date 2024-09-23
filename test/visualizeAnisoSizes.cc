#include <apf.h>
#include <apfMDS.h>
#include <maStats.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <maDBG.h>
#include <lionPrint.h>

#ifdef HAVE_SIMMETRIX
#include <ph.h>
#include <apfSIM.h>
#include <gmi_sim.h>
#include <phastaChef.h>
#include <SimPartitionedMesh.h>
#endif

#include <stdlib.h>
#include <sstream>
#include <fstream>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
// ===============================

void safe_mkdir(const char* path);
double getLargetsSize(
    apf::Mesh2* m,
    apf::Field* sizes);
void visualizeSizeField(
    const char* modelFile,
    const char* meshFile,
    const char* sizeName,
    const char* frameName,
    int sampleSize[2],
    double userScale,
    const char* outputPrefix,
    pcu::PCU *PCUObj);

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if (argc < 8) {
    if (PCUObj.Self() == 0) {
      printf("USAGE1: %s <mesh.smb> <output_prefix> <scale field name>"
          "<frames field name> <n_u> <n_v> <scale>\n", argv[0]);
      printf("USAGE2: %s <mesh.sms> <output_prefix> <scale field name>"
          "<frames field name> <n_u> <n_v> <scale>\n", argv[0]);
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_SIMMETRIX
  SimModel_start();
  Sim_readLicenseFile(0);
  SimPartitionedMesh_start(0, 0);
  gmi_sim_start();
  gmi_register_sim();
#endif

  gmi_register_mesh();
  gmi_register_null();

  const char* meshFile = argv[1];
  const char* inPrefix = argv[2];
  const char* sizeName = argv[3];
  const char* frameName = argv[4];
  /* int n_u = atoi(argv[5]); */
  /* int n_v = atoi(argv[6]); */

  int sampleSize[2] = {atoi(argv[5]),atoi(argv[6])};
  double scale = atof(argv[7]);
  safe_mkdir(inPrefix);
  visualizeSizeField(".null", meshFile, sizeName, frameName, sampleSize, scale, inPrefix, &PCUObj);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  SimModel_stop();
  Sim_unregisterAllKeys();
#endif

  }
  MPI_Finalize();
}

void safe_mkdir(const char* path)
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

double getLargetsSize(
    apf::Mesh2* m,
    apf::Field* sizes)
{
  double maxSize = 0.0;
  apf::MeshEntity* vert;
  apf::MeshIterator* it = m->begin(0);
  while ( (vert = m->iterate(it)) ) {
    apf::Vector3 scales;
    apf::getVector(sizes, vert, 0, scales);
    if (scales[0] > maxSize)
      maxSize = scales[0];
    if (scales[1] > maxSize)
      maxSize = scales[1];
    if (scales[2] > maxSize)
      maxSize = scales[2];
  }
  m->end(it);
  m->getPCU()->Max<double>(&maxSize, 1);
  return maxSize;
}

void visualizeSizeField(
    const char* modelFile,
    const char* meshFile,
    const char* sizeName,
    const char* frameName,
    int sampleSize[2],
    double userScale,
    const char* outputPrefix,
    pcu::PCU *PCUObj)
{
  apf::Mesh2* m;
#ifdef HAVE_SIMMETRIX
  /* if it is a simmetrix mesh */
  if (ph::mesh_has_ext(meshFile, "sms")) {
    pParMesh sim_mesh = PM_load(meshFile, NULL, NULL);
    m = apf::createMesh(sim_mesh, PCUObj);
  } else
#endif
  {
    // load the mesh change to desired order and write as before vtks
    m = apf::loadMdsMesh(modelFile,meshFile,PCUObj);
  }
  m->verify();

  apf::Field* sizes;
  apf::Field* frames;
#ifdef HAVE_SIMMETRIX
  /* if it is a simmetrix mesh */
  if (ph::mesh_has_ext(meshFile, "sms")) {
    if(m->findField("sizes")) apf::destroyField(m->findField("sizes"));
    if(m->findField("frames")) apf::destroyField(m->findField("frames"));
    sizes  = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
    frames = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
    ph::attachSIMSizeField(m, sizes, frames);
  } else
#endif
  {
    char message[512];
    // first find the sizes field
    sizes  = m->findField(sizeName);
    snprintf(message, 512, "Couldn't find a field with name %s in mesh!", sizeName);
    PCU_ALWAYS_ASSERT_VERBOSE(sizes, message);

    // then find the frames field if they exist
    frames = m->findField(frameName);
    snprintf(message, 512, "Couldn't find a field with name %s in mesh!", frameName);
    PCU_ALWAYS_ASSERT_VERBOSE(frames, message);
  }

  // remove every field except for sizes and frames
  int index = 0;
  while (m->countFields() > 2) {
    apf::Field* f = m->getField(index);
    if (f == sizes || f == frames) {
      index++;
      continue;
    }
    m->removeField(f);
    apf::destroyField(f);
  }

  m->verify();

  std::stringstream ss;
  ss << outputPrefix << "/size_field_vis";
  ma_dbg::visualizeSizeField(
      m, sizes, frames, sampleSize, userScale, ss.str().c_str());
  ss.str("");

  m->destroyNative();
  apf::destroyMesh(m);
}
