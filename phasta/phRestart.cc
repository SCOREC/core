#include "phRestart.h"
#include <apf.h>
#include "phIO.h"
#include <PCU.h>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace ph {

void attachField(
    apf::Mesh* m,
    const char* fieldname,
    double* data,
    int size)
{
  apf::Field* f = apf::createPackedField(m, fieldname, size);
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it))) {
    apf::setComponents(f, e, 0, data + i);
    i += size;
  }
  m->end(it);
}

void readAndAttachField(
    apf::Mesh* m,
    const char* filename,
    const char* fieldname,
    int size)
{
  double* data;
  ph_read_field(filename, fieldname, &data);
  attachField(m, fieldname, data, size);
  free(data);
}

/* silliest darn field I ever did see */
static double* buildMappingData(apf::Mesh* m)
{
  int n = m->count(0);
  /* malloc instead of new[] for consistency with ph_read_field */
  double* data = (double*)malloc(sizeof(double) * n);
  int self = PCU_Comm_Self();
  for (int i = 0; i < n; ++i) {
    data[i * 2] = self;
    data[i * 2 + 1] = i;
  }
  return data;
}

static int readStepNum()
{
  std::ifstream tinyFile("numstart.dat");
  int step;
  tinyFile >> step;
  return step;
}

static std::string buildRestartFileName(Input& in)
{
  std::stringstream ss;
  int step = readStepNum();
  int rank = PCU_Comm_Self() + 1;
  ss << in.restartFileName << '.' << step << '.' << rank;
  return ss.str();
}

void readAndAttachSolution(Input& in, apf::Mesh* m)
{
  std::string filename = buildRestartFileName(in);
  readAndAttachField(m, filename.c_str(), "solution", in.ensa_dof);
  if (in.displacementMigration)
    readAndAttachField(m, filename.c_str(), "displacement", 3);
  if (in.dwalMigration)
    readAndAttachField(m, filename.c_str(), "dwal", 1);
  if (in.buildMapping && !in.adaptFlag) {
    double* mapping = buildMappingData(m);
    attachField(m, "mapping", mapping, 2);
    free(mapping);
  }
}

}
