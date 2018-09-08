#include <apf.h>
#include <apfMDS.h>
#include <maStats.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <PCU.h>

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

static double PI = 3.14159265359;

void safe_mkdir(const char* path);
double getLargetsSize(
    apf::Mesh2* m,
    apf::Field* sizes);
apf::Vector3 getPointOnEllipsoid(
    apf::Vector3 center,
    apf::Vector3 abc,
    apf::Matrix3x3 rotation,
    double scaleFactor,
    double u,
    double v);
void makeEllipsoid(
    apf::Mesh2* msf,
    apf::Mesh2* m,
    apf::Field* sizes,
    apf::Field* frames,
    apf::MeshEntity* vert,
    double scaleFactor,
    int sampleSize[2]);
void visualizeSizeField(
    const char* modelFile,
    const char* meshFile,
    const char* sizeName,
    const char* frameName,
    int sampleSize[2],
    double userScale,
    const char* outputPrefix);

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (argc < 8) {
    if (PCU_Comm_Self() == 0) {
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
  visualizeSizeField(".null", meshFile, sizeName, frameName, sampleSize, scale, inPrefix);

#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  SimModel_stop();
  Sim_unregisterAllKeys();
#endif

  PCU_Comm_Free();
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
  PCU_Max_Doubles(&maxSize, 1);
  return maxSize;
}

apf::Vector3 getPointOnEllipsoid(
    apf::Vector3 center,
    apf::Vector3 abc,
    apf::Matrix3x3 rotation,
    double scaleFactor,
    double u,
    double v)
{
  apf::Vector3 result;
  result[0] = abc[0] * cos(u) * cos(v);
  result[1] = abc[1] * cos(u) * sin(v);
  result[2] = abc[2] * sin(u);

  result = result * scaleFactor;

  result = rotation * result + center;
  return result;
}


void makeEllipsoid(
    apf::Mesh2* msf,
    apf::Mesh2* mesh,
    apf::Field* sizes,
    apf::Field* frames,
    apf::MeshEntity* vert,
    double scaleFactor,
    int sampleSize[2])
{

  apf::Vector3 center;
  mesh->getPoint(vert, 0, center);

  apf::Vector3 abc;
  apf::getVector(sizes, vert, 0, abc);

  apf::Matrix3x3 rotation;
  apf::getMatrix(frames, vert, 0, rotation);


  double U0 = 0.0;
  double U1 = 2 * PI;
  double V0 = -PI/2.;
  double V1 =  PI/2.;
  int n = sampleSize[0];
  int m = sampleSize[1];
  double dU = (U1 - U0) / (n-1);
  double dV = (V1 - V0) / (m-1);

  // make the array of vertex coordinates in the physical space
  std::vector<ma::Vector> ps;
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      double u = U0 + i * dU;
      double v = V0 + j * dV;
      apf::Vector3 pt = getPointOnEllipsoid(center, abc, rotation, scaleFactor, u, v);
      ps.push_back(pt);
    }
  }
  // make the vertexes and set the coordinates using the array
  std::vector<apf::MeshEntity*> vs;
  for (size_t i = 0; i < ps.size(); i++) {
    apf::MeshEntity* newVert = msf->createVert(0);
    msf->setPoint(newVert, 0, ps[i]);
    vs.push_back(newVert);
  }

  PCU_ALWAYS_ASSERT(vs.size() == ps.size());

  apf::MeshEntity* v[3];
  // make the lower/upper t elems
  for (int i = 0; i < n-1; i++) {
    for (int j = 0; j < m-1; j++) {
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 0)];
      v[1] = vs[(i + 0) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(msf, 0, apf::Mesh::TRIANGLE, v);
      // upper triangle
      v[0] = vs[(i + 0) + n * (j + 1)];
      v[1] = vs[(i + 1) + n * (j + 1)];
      v[2] = vs[(i + 1) + n * (j + 0)];
      apf::buildElement(msf, 0, apf::Mesh::TRIANGLE, v);
    }
  }
}


void visualizeSizeField(
    const char* modelFile,
    const char* meshFile,
    const char* sizeName,
    const char* frameName,
    int sampleSize[2],
    double userScale,
    const char* outputPrefix)
{
  apf::Mesh2* m;
#ifdef HAVE_SIMMETRIX
  /* if it is a simmetrix mesh */
  if (ph::mesh_has_ext(meshFile, "sms")) {
    pParMesh sim_mesh = PM_load(meshFile, NULL, NULL);
    m = apf::createMesh(sim_mesh);
  } else
#endif
  {
    // load the mesh change to desired order and write as before vtks
    m = apf::loadMdsMesh(modelFile,meshFile);
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
    sprintf(message, "Couldn't find a field with name %s in mesh!", sizeName);
    PCU_ALWAYS_ASSERT_VERBOSE(sizes, message);

    // then find the frames field if they exist
    frames = m->findField(frameName);
    sprintf(message, "Couldn't find a field with name %s in mesh!", frameName);
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

  // create the size-field visualization mesh
  apf::Mesh2* msf = apf::makeEmptyMdsMesh(gmi_load(".null"), 2, false);

  apf::MeshEntity* vert;
  apf::MeshIterator* it = m->begin(0);
  while ( (vert = m->iterate(it)) )
    if (m->isOwned(vert))
      makeEllipsoid(msf, m, sizes, frames, vert, userScale , sampleSize);
  m->end(it);

  std::stringstream ss;
  ss << outputPrefix << "/size_field_vis";
  apf::writeVtkFiles(ss.str().c_str(), msf);
  ss.str("");

  msf->destroyNative();
  apf::destroyMesh(msf);
  m->destroyNative();
  apf::destroyMesh(m);
}
