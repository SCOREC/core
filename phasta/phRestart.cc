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
  assert(i == size * m->count(0));
}

void detachField(
    apf::Mesh* m,
    const char* fieldname,
    double*& data,
    int& size)
{
  apf::Field* f = m->findField(fieldname);
  size = apf::countComponents(f);
  data = (double*)malloc(sizeof(double) * size * m->count(0));
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
  size_t i = 0;
  while ((e = m->iterate(it))) {
    apf::getComponents(f, e, 0, data + i);
    i += size;
  }
  m->end(it);
  assert(i == size * m->count(0));
}

void readAndAttachField(
    Input& in,
    apf::Mesh* m,
    const char* filename,
    const char* fieldname)
{
  double* data;
  int nodes, vars, step;
  ph_read_field(filename, fieldname, &data,
      &nodes, &vars, &step);
  assert(nodes == static_cast<int>(m->count(0)));
  assert(step == in.timeStepNumber);
  attachField(m, fieldname, data, vars);
  free(data);
}

void detachAndWriteField(
    Input& in,
    apf::Mesh* m,
    FILE* f,
    const char* fieldname)
{
  double* data;
  int size;
  detachField(m, fieldname, data, size);
  ph_write_field(f, fieldname, data, m->count(0), size, in.timeStepNumber);
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

static void readStepNum(Input& in)
{
  std::ifstream tinyFile("numstart.dat");
  int step;
  tinyFile >> step;
  assert(in.timeStepNumber == step);
}

static std::string buildRestartFileName(Input& in)
{
  std::stringstream ss;
  int rank = PCU_Comm_Self() + 1;
  ss << in.restartFileName << '.' << in.timeStepNumber << '.' << rank;
  return ss.str();
}

void readAndAttachSolution(Input& in, apf::Mesh* m)
{
  readStepNum(in);
  std::string filename = buildRestartFileName(in);
  readAndAttachField(in, m, filename.c_str(), "solution");
  if (in.displacementMigration)
    readAndAttachField(in, m, filename.c_str(), "displacement");
  if (in.dwalMigration)
    readAndAttachField(in, m, filename.c_str(), "dwal");
  if (in.buildMapping && !in.adaptFlag) {
    double* mapping = buildMappingData(m);
    attachField(m, "mapping", mapping, 2);
    free(mapping);
  }
}

void detachAndWriteSolution(Input& in, apf::Mesh* m, std::string path)
{
  path += buildRestartFileName(in);
  FILE* f = fopen(path.c_str(), "w");
  ph_write_preamble(f);
  int nodes = m->count(0);
  ph_write_header(f, "number of modes", 0, 1, &nodes);
  ph_write_header(f, "number of variables", 0, 1, &in.ensa_dof);
  detachAndWriteField(in, m, f, "solution");
  if (in.displacementMigration)
    detachAndWriteField(in, m, f, "displacement");
  if (in.dwalMigration)
    detachAndWriteField(in, m, f, "dwal");
  if (in.buildMapping && !in.adaptFlag)
    detachAndWriteField(in, m, f, "mapping");
  fclose(f);
}

}
