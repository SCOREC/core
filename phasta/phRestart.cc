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
  size_t n = m->count(0);
  apf::NewArray<double> c(size);
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    for (int j = 0; j < size; ++j)
      c[j] = data[j * n + i];
    apf::setComponents(f, e, 0, &c[0]);
    ++i;
  }
  m->end(it);
  assert(i == n);
}

void detachField(
    apf::Field* f,
    double*& data,
    int& size)
{
  apf::Mesh* m = apf::getMesh(f);
  size = apf::countComponents(f);
  size_t n = m->count(0);
  apf::NewArray<double> c(size);
  data = (double*)malloc(sizeof(double) * size * m->count(0));
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    apf::getComponents(f, e, 0, &c[0]);
    for (int j = 0; j < size; ++j)
      data[j * n + i] = c[j];
    ++i;
  }
  m->end(it);
  assert(i == n);
}

void detachField(
    apf::Mesh* m,
    const char* fieldname,
    double*& data,
    int& size)
{
  apf::Field* f = m->findField(fieldname);
  detachField(f, data, size);
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

/* silliest darn fields I ever did see */
static double* buildMappingPartId(apf::Mesh* m)
{
  int n = m->count(0);
  /* malloc instead of new[] for consistency with ph_read_field */
  double* data = (double*)malloc(sizeof(double) * n);
  int self = PCU_Comm_Self();
  for (int i = 0; i < n; ++i)
    data[i] = self;
  return data;
}
static double* buildMappingVtxId(apf::Mesh* m)
{
  int n = m->count(0);
  /* malloc instead of new[] for consistency with ph_read_field */
  double* data = (double*)malloc(sizeof(double) * n);
  for (int i = 0; i < n; ++i)
    data[i] = i;
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
  double t0 = MPI_Wtime();
  readStepNum(in);
  std::string filename = buildRestartFileName(in);
  readAndAttachField(in, m, filename.c_str(), "solution");
  if (in.displacementMigration)
    readAndAttachField(in, m, filename.c_str(), "displacement");
  if (in.dwalMigration)
    readAndAttachField(in, m, filename.c_str(), "dwal");
  if (in.buildMapping) {
    double* mapping = buildMappingPartId(m);
    attachField(m, "mapping_partid", mapping, 1);
    free(mapping);
    mapping = buildMappingVtxId(m);
    attachField(m, "mapping_vtxid", mapping, 1);
    free(mapping);
  }
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    printf("solution read and attached in %f seconds\n", t1 - t0);
}

void detachAndWriteSolution(Input& in, apf::Mesh* m, std::string path)
{
  double t0 = MPI_Wtime();
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
  if (in.buildMapping) {
    detachAndWriteField(in, m, f, "mapping_partid");
    detachAndWriteField(in, m, f, "mapping_vtxid");
  }
  fclose(f);
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    printf("solution written in %f seconds\n", t1 - t0);
}

}
