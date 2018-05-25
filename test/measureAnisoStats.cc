#include <apf.h>
#include <apfMDS.h>
#include <apfShape.h>
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


void safe_mkdir(const char* path);
void writeTable(const char* outfile,
    const std::vector<std::vector<double> > & table);
void getStats(
    const char* modelFile, const char* meshFile,
    const char* sizeName, const char* frameName,
    const char* outputPrefix);

int main(int argc, char** argv)
{

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if (argc < 5) {
    if (PCU_Comm_Self() == 0) {
      printf("USAGE1: %s <mesh.smb> <output_prefix> <scale field name>"
          "<frames field name>\n", argv[0]);
      printf("USAGE2: %s <mesh.sms> <output_prefix> <scale field name>"
          "<frames field name>\n", argv[0]);
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

  safe_mkdir(inPrefix);

  getStats(".null", meshFile, sizeName, frameName, inPrefix);

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


void writeTable(const char* outFile,
    const std::vector<std::vector<double> > & table)
{
  std::ofstream out;
  out.open(outFile);
  for (size_t i = 0; i < table.size(); i++) {
    std::vector<double> r = table[i];
    for (size_t j = 0; j < r.size() - 1; j++) {
      out << r[j] << " ";
    }
    out << r[r.size() - 1];
    out << std::endl;
  }
  out.close();
}

void getStats(
    const char* modelFile, const char* meshFile,
    const char* sizeName, const char* frameName,
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


  if (PCU_Comm_Self() == 0)
    printf("\n evaluating the statistics! \n");
  // get the stats
  ma::SizeField* sf = ma::makeSizeField(m, sizes, frames, true);
  std::vector<double> el, lq;
  ma::stats(m, sf, el, lq, true);


  // create field for visualizaition
  apf::Field* f_lq = apf::createField(m, "linear_quality", apf::SCALAR, apf::getConstant(m->getDimension()));

  // attach cell-based mesh quality
  int n;
  apf::MeshEntity* r;
  if (m->getDimension() == 3)
    n = apf::countEntitiesOfType(m, apf::Mesh::TET);
  else
    n = apf::countEntitiesOfType(m, apf::Mesh::TRIANGLE);
  size_t i = 0;
  apf::MeshIterator* rit = m->begin(m->getDimension());
  while ((r = m->iterate(rit))) {
    if (! apf::isSimplex(m->getType(r))) {// ignore non-simplex elements
      apf::setScalar(f_lq, r, 0, 100.0); // set as 100
    }
    else {
      apf::setScalar(f_lq, r, 0, lq[i]);
      ++i;
    }
  }
  m->end(rit);
  PCU_ALWAYS_ASSERT(i == (size_t) n);


  std::vector<std::vector<double> > qtable;
  for (size_t i = 0; i < lq.size(); i++) {
    std::vector<double> r;
    r.push_back(lq[i]);
    qtable.push_back(r);
  }

  std::vector<std::vector<double> > etable;
  for (size_t i = 0; i < el.size(); i++) {
    std::vector<double> r;
    r.push_back(el[i]);
    etable.push_back(r);
  }

  std::stringstream ss, sse, ssq, ssm, ssl;
  ss << outputPrefix << "/" << "linear_tables";
  const char* pathName = ss.str().c_str();
  safe_mkdir(pathName);

  ssq << pathName << "/linearQTable_" << PCU_Comm_Self() << ".dat";
  sse << pathName << "/linearETable_" << PCU_Comm_Self() << ".dat";
  ssm << pathName << "/mesh_quality_vis";
  ssl << pathName << "/mesh_edge_length_vis";

  writeTable(ssq.str().c_str(), qtable);
  writeTable(sse.str().c_str(), etable);
  apf::writeVtkFiles(ssm.str().c_str(), m);

  // create field for visualizaition of edge lengths
  apf::Field* f_el = apf::createField(m, "edge_length", apf::SCALAR, apf::getConstant(1));

  // attach edge-based lengths
  apf::MeshEntity* e;
  apf::MeshIterator* it = m->begin(1);
  i = 0;
  while ((e = m->iterate(it))) {
    apf::setScalar(f_el, e, 0, el[i]);
    ++i;
  }
  m->end(it);

  for (int dim = m->getDimension(); dim > 1  ; dim--) {
    it = m->begin(dim);
    while ( (e = m->iterate(it)) ){
      	m->destroy(e);
    }
    m->end(it);
  }
  apf::writeVtkFiles(ssl.str().c_str(), m);


  m->destroyNative();
  apf::destroyMesh(m);
}
