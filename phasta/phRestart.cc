#include <lionPrint.h>
#include "phRestart.h"
#include <apf.h>
#include <apfField.h>
#include "phIO.h"
#include "phiotimer.h"
#include "apfShape.h"
#include "ph.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <pcu_util.h>
#include <cstring>
#ifdef HAVE_SIMMETRIX
#include <apfSIM.h>
#endif

namespace ph {

/* This function is used by phastaChef driver
   to extract a certain field from solution field.
   can also extract sim field from apf field. */
apf::Field* extractField(apf::Mesh* m,
    const char* packedFieldname,
    const char* requestFieldname,
    int firstComp,
    int valueType,
    bool simField)
{
  apf::Field* f = m->findField(packedFieldname);
  if(!f && m->getPCU()->Self() == 0)
    lion_eprint(1, "No packed field \"%s\"", packedFieldname);
  PCU_ALWAYS_ASSERT(f);
  apf::Field* rf = m->findField(requestFieldname);
  if (rf)
    apf::destroyField(rf);
  int numOfComp = 0;
  if (valueType == apf::SCALAR)
    numOfComp = 1;
  else if (valueType == apf::VECTOR)
    numOfComp= 3;
  else
    PCU_ALWAYS_ASSERT(valueType == apf::SCALAR || valueType == apf::VECTOR);
#ifdef HAVE_SIMMETRIX
  if (simField) {
    rf = apf::createSIMFieldOn(m, requestFieldname, valueType);
  } else
#else
  (void)simField;
#endif
  {
    rf = apf::createFieldOn(m, requestFieldname, valueType);
  }
  apf::NewArray<double> inVal(apf::countComponents(f));
  apf::NewArray<double> outVal(numOfComp);
  int endComp = firstComp + numOfComp - 1;
  PCU_ALWAYS_ASSERT(firstComp >= 1);
  PCU_ALWAYS_ASSERT(endComp <= apf::countComponents(f));
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  while ((vtx = m->iterate(it))) {
    apf::getComponents(f, vtx, 0, &inVal[0]);
    int j = 0;
    for (int i = firstComp-1; i < endComp; i++){
      outVal[j] = inVal[i];
      j++;
    }
    PCU_ALWAYS_ASSERT(j == numOfComp);
    apf::setComponents(rf,vtx, 0, &outVal[0]);
  }
  m->end(it);
  return rf;
}

/* default combine 3 fields into 1 */
apf::Field* combineField(apf::Mesh* m,
    const char* packedFieldname,
    const char* inFieldname1,
    const char* inFieldname2,
    const char* inFieldname3)
{
  apf::Field* f1 = m->findField(inFieldname1);
  apf::Field* f2 = m->findField(inFieldname2);
  apf::Field* f3 = m->findField(inFieldname3);
  PCU_ALWAYS_ASSERT(f1);
  PCU_ALWAYS_ASSERT(f2);
  PCU_ALWAYS_ASSERT(f3);
  int in_size1 = apf::countComponents(f1);
  int in_size2 = apf::countComponents(f2);
  int in_size3 = apf::countComponents(f3);
  int out_size = in_size1 + in_size2 + in_size3;
  apf::Field* rf = m->findField(packedFieldname);
  if (rf)
    apf::destroyField(rf);
  rf = apf::createPackedField(m, packedFieldname, out_size);
  apf::NewArray<double> inVal1(in_size1);
  apf::NewArray<double> inVal2(in_size2);
  apf::NewArray<double> inVal3(in_size3);
  apf::NewArray<double> outVal(out_size);
  /* copy data */
  apf::MeshEntity* vtx;
  apf::MeshIterator* it = m->begin(0);
  while ((vtx = m->iterate(it))) {
    apf::getComponents(f1, vtx, 0, &inVal1[0]);
    apf::getComponents(f2, vtx, 0, &inVal2[0]);
    apf::getComponents(f3, vtx, 0, &inVal3[0]);
/* simpler, test later, need algorithm header */
/*
    copy(inVal1, inVal2 + in_size1, outVal);
    copy(inVal2, inVal2 + in_size2, outVal + in_size1);
    copy(inVal3, inVal2 + in_size3, outVal + in_size2);
*/
    int j = 0;
    int i = 0;
    for (i = 0; i < in_size1; i++){
      outVal[j] = inVal1[i];
      j++;
    }
    for (i = 0; i < in_size2; i++){
      outVal[j] = inVal2[i];
      j++;
    }
    for (i = 0; i < in_size3; i++){
      outVal[j] = inVal3[i];
      j++;
    }
    PCU_ALWAYS_ASSERT(j == out_size);
    apf::setComponents(rf, vtx, 0, &outVal[0]);
  }
  m->end(it);
  /* destroy input fields */
  apf::destroyField(f1);
  apf::destroyField(f2);
  apf::destroyField(f3);
  return rf;
}

/* in_size is the number of dofs for the data array
   and out_size is the number of dofs in the field.
   they are unequal when we read a restart file
   and want to add more dofs for the next restart files */

void attachField(
    apf::Mesh* m,
    const char* fieldname,
    double* data,
    int in_size,
    int out_size)
{
  if (!(in_size <= out_size))
    lion_eprint(1, "field \"%s\" in_size %d out_size %d\n", fieldname, in_size, out_size);
  PCU_ALWAYS_ASSERT(in_size <= out_size);
  apf::Field* f = m->findField(fieldname);
  if( f )
    apf::destroyField(f);
  f = apf::createPackedField(m, fieldname, out_size);
  size_t n = m->count(0);
  apf::NewArray<double> c(out_size);
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(0);
  while ((e = m->iterate(it))) {
    for (int j = 0; j < in_size; ++j)
      c[j] = data[j * n + i];
    apf::setComponents(f, e, 0, &c[0]);
    ++i;
  }
  m->end(it);
  PCU_ALWAYS_ASSERT(i == n);
}

bool attachRandField(
    Input& in,
    const char* fieldname,
    double* data,
    int nnodes,
    int nvars)
{
  if(!strcmp(fieldname, "rbParams")) {
    in.nRigidBody = nnodes;
    in.nRBParam   = nvars;
    in.rbParamData.clear();
    for (int i = 0; i < nnodes; i++) {
      for (int j = 0; j < nvars; j++) {
        in.rbParamData.push_back(data[j*nnodes + i]);
      }
    }
    PCU_ALWAYS_ASSERT((size_t) nnodes * nvars == in.rbParamData.size());
    return true;
  }
  return false;
}

void attachCellField(
    apf::Mesh* m,
    const char* fieldname,
    double* data,
    int in_size,
    int out_size)
{
  if (!(in_size <= out_size))
    lion_eprint(1, "field \"%s\" in_size %d out_size %d\n", fieldname, in_size, out_size);
  PCU_ALWAYS_ASSERT(in_size <= out_size);
  apf::Field* f = m->findField(fieldname);
  if( f )
    apf::destroyField(f);
  f = apf::createPackedField(m, fieldname, out_size, apf::getConstant(m->getDimension()));
  size_t n = m->count(m->getDimension());
  apf::NewArray<double> c(out_size);
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it))) {
    for (int j = 0; j < in_size; ++j)
      c[j] = data[j * n + i];
    apf::setComponents(f, e, 0, &c[0]);
    ++i;
  }
  m->end(it);
  PCU_ALWAYS_ASSERT(i == n);
}

/* convenience wrapper, in most cases in_size=out_size */
void attachField(
    apf::Mesh* m,
    const char* fieldname,
    double* data,
    int size)
{
  attachField(m, fieldname, data, size, size);
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
  PCU_ALWAYS_ASSERT(i == n);
  apf::destroyField(f);
}

void detachField(
    apf::Mesh* m,
    const char* fieldname,
    double*& data,
    int& size)
{
  apf::Field* f = m->findField(fieldname);
  PCU_ALWAYS_ASSERT(f);
  detachField(f, data, size);
}

static bool isNodalField(const char* fieldname, int nnodes, apf::Mesh* m)
{
  static char const* const known_nodal_fields[] = {
    "solution",
    "displacement",
    "dwal",
    "mapping_partid",
    "mapping_vtxid",
    "errors",
    "time derivative of solution",
    "motion_coords",
    "mesh_vel",
    "tb_factor",
    "residual",
    "dc_lag",
    "ybar",
    "wss",
    "wssbar",
    "pressure projection vectors",
    "projection vectors",
    "vorticity",
    "BCs",
    "cdelsq"
  };
  static char const* const known_cell_fields[] = {
    "VOF solution",
    "meshQ",
    "meshCFL",
    "VMS_error",
    "err_tri_f",
    "material_type",
    "lesnut"
  };
  static char const* const known_rand_fields[] = {
    "rbParams"
  };
  int known_nodal_field_count =
    sizeof(known_nodal_fields) / sizeof(known_nodal_fields[0]);
  int known_cell_field_count =
    sizeof(known_cell_fields) / sizeof(known_cell_fields[0]);
  int known_rand_field_count =
    sizeof(known_rand_fields) / sizeof(known_rand_fields[0]);
  for (int i = 0; i < known_nodal_field_count; ++i)
    if (!strcmp(fieldname, known_nodal_fields[i])) {
      PCU_ALWAYS_ASSERT(static_cast<size_t>(nnodes) == m->count(0));
      return true;
    }
  for (int i = 0; i < known_cell_field_count; ++i)
    if (!strcmp(fieldname, known_cell_fields[i])) {
      PCU_ALWAYS_ASSERT(static_cast<size_t>(nnodes) == m->count(m->getDimension()));
      return false;
    }
  for (int i = 0; i < known_rand_field_count; ++i)
    if (!strcmp(fieldname, known_rand_fields[i]))
      return false;
  if( !m->getPCU()->Self() ) {
    lion_eprint(1, "unknown restart field name \"%s\"\n", fieldname);
    lion_eprint(1, "please add \"%s\" to isNodalField above line %d of %s\n",
        fieldname, __LINE__, __FILE__);
  }
  if (static_cast<size_t>(nnodes) == m->count(0)) {
    lion_eprint(1, "assuming \"%s\" is a nodal field,\n"
                    "it is the right size...\n", fieldname);
    return true;
  }
  return false;
}

int readAndAttachField(
    Input& in,
    FILE* f,
    apf::Mesh* m,
    int swap)
{
  double* data;
  int nodes, vars, step;
  char hname[1024];
  const char* anyfield = "";
  int ret = ph_read_field(f, m->getPCU()->GetCHandle(), 
      anyfield, swap, &data, &nodes, &vars, &step, hname);
  /* no field was found or the field has an empty data block */
  if(ret==0 || ret==1)
    return ret;
  if (!isNodalField(hname, nodes, m)) {
    if (attachRandField(in, hname, data, nodes, vars)) {
      free(data);
      return 1;
    }
    attachCellField(m, hname, data, vars, vars);
    free(data);
    return 1;
  }
  PCU_ALWAYS_ASSERT(step == in.timeStepNumber);
  int out_size = vars;
  if ( std::string(hname) == std::string("solution") )
    out_size = in.ensa_dof;
  if (m->findField(hname)) {
    if (!m->getPCU()->Self())
      lion_eprint(1, "field \"%s\" already attached to the mesh, "
                      "ignoring request to re-attach...\n", hname);
  } else {
    attachField(m, hname, data, vars, out_size);
  }
  free(data);
  return 1;
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

void detachAndWriteCellField(
    Input& in,
    apf::Mesh* m,
    FILE* file,
    const char* fieldname)
{
  apf::Field* f = m->findField(fieldname);
  PCU_ALWAYS_ASSERT(f);
  double* data;
  int size = apf::countComponents(f);
  size_t n = m->count(m->getDimension());
  apf::NewArray<double> c(size);
  data = (double*)malloc(sizeof(double) * size * m->count(m->getDimension()));
  apf::MeshEntity* e;
  size_t i = 0;
  apf::MeshIterator* it = m->begin(m->getDimension());
  while ((e = m->iterate(it))) {
    apf::getComponents(f, e, 0, &c[0]);
    for (int j = 0; j < size; ++j)
      data[j * n + i] = c[j];
    ++i;
  }
  m->end(it);
  PCU_ALWAYS_ASSERT(i == n);
  apf::destroyField(f);
  ph_write_field(file, fieldname, data, m->count(m->getDimension()), size, in.timeStepNumber);
  free(data);
}

void detachAndWriteRandField(
    Input& in,
    FILE* f,
    const char* fieldname)
{
  if (!strcmp(fieldname, "rbParams")) {
    int nnodes = in.nRigidBody;
    int nvars  = in.nRBParam;
    double* data = (double*) malloc(sizeof(double) * nnodes * nvars);
    size_t iv = 0;
    for (int i = 0; i< nnodes; i++) {
      for (int j = 0; j < nvars; j++) {
        data[j*nnodes + i] = in.rbParamData[iv];
        iv++;
      }
    }
    ph_write_field(f, fieldname, data, nnodes, nvars, in.timeStepNumber);
    free(data);
    in.rbParamData.clear();
  }
}

/* silliest darn fields I ever did see */
static double* buildMappingPartId(apf::Mesh* m)
{
  int n = m->count(0);
  /* malloc instead of new[] for consistency with ph_read_field */
  double* data = (double*)malloc(sizeof(double) * n);
  int self = m->getPCU()->Self();
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

static std::string buildRestartFileName(std::string prefix, int step, pcu::PCU *PCUObj)
{
  std::stringstream ss;
  int rank = PCUObj->Self() + 1;
  ss << prefix << '.' << step << '.' << rank;
  return ss.str();
}

void readAndAttachFields(Input& in, apf::Mesh* m) {
  PCU_t h;
  h.ptr = static_cast<void*>(m->getPCU());
  phastaio_initStats(h);
  double t0 = pcu::Time();
  setupInputSubdir(in.restartFileName, m->getPCU());
  std::string filename = buildRestartFileName(in.restartFileName, in.timeStepNumber, m->getPCU());
  phastaio_setfile(RESTART_READ);
  FILE* f = in.openfile_read(in, filename.c_str(), m->getPCU());
  if (!f) {
    lion_eprint(1,"failed to open \"%s\"!\n", filename.c_str());
    abort();
  }
  int swap = ph_should_swap(f, m->getPCU()->GetCHandle());
  /* stops when ph_read_field returns 0 */
  while( readAndAttachField(in,f,m,swap) ) {}
  PHASTAIO_CLOSETIME(fclose(f);)
  double t1 = pcu::Time();
  if (!m->getPCU()->Self())
    lion_oprint(1,"fields read and attached in %f seconds\n", t1 - t0);
  if(in.printIOtime) phastaio_printStats(h);
}

static void destroyIfExists(apf::Mesh* m, const char* name)
{
  apf::Field* f = m->findField(name);
  if (f)
    apf::destroyField(f);
}

void buildMapping(apf::Mesh* m)
{
  destroyIfExists(m, "mapping_partid");
  double* mapping = buildMappingPartId(m);
  attachField(m, "mapping_partid", mapping, 1);
  free(mapping);
  destroyIfExists(m, "mapping_vtxid");
  mapping = buildMappingVtxId(m);
  attachField(m, "mapping_vtxid", mapping, 1);
  free(mapping);
}

void attachZeroSolution(Input& in, apf::Mesh* m)
{
  int vars = in.ensa_dof;
  int nodes = m->count(0);
  double* data = new double[nodes * vars]();
  attachField(m, "solution", data, vars);
  delete [] data;
}

void detachAndWriteSolution(Input& in, Output& out, apf::Mesh* m, std::string path)
{
  double t0 = pcu::Time();
  path += buildRestartFileName("restart", in.timeStepNumber, m->getPCU());
  phastaio_setfile(RESTART_WRITE);
  FILE* f = out.openfile_write(out, path.c_str());
  if (!f) {
    lion_eprint(1,"failed to open \"%s\"!\n", path.c_str());
    abort();
  }
  ph_write_preamble(f);
  int nodes = m->count(0);
  ph_write_header(f, "number of modes", 0, 1, &nodes);
  ph_write_header(f, "number of variables", 0, 1, &in.ensa_dof);
  apf::Field* errField = m->findField("errors");
  if (errField)
    apf::destroyField(errField);
  if (m->findField("solution"))
    detachAndWriteField(in, m, f, "solution");
  if (m->findField("time derivative of solution"))
    detachAndWriteField(in, m, f, "time derivative of solution");
  if (m->findField("motion_coords"))
    detachAndWriteField(in, m, f, "motion_coords");
  if (m->findField("mesh_vel"))
    detachAndWriteField(in, m, f, "mesh_vel");
  if (m->findField("dc_lag"))
    detachAndWriteField(in, m, f, "dc_lag");
  if (m->findField("pressure projection vectors"))
    detachAndWriteField(in, m, f, "pressure projection vectors");
  if (m->findField("BCs"))
    detachAndWriteField(in, m, f, "BCs");
  if (m->findField("cdelsq"))
    detachAndWriteField(in, m, f, "cdelsq");
  if (in.displacementMigration)
    detachAndWriteField(in, m, f, "displacement");
  if (in.dwalMigration)
    detachAndWriteField(in, m, f, "dwal");
  if (in.buildMapping) {
    detachAndWriteField(in, m, f, "mapping_partid");
    detachAndWriteField(in, m, f, "mapping_vtxid");
  }
  if (m->findField("err_tri_f"))
    detachAndWriteCellField(in, m, f, "err_tri_f");
  if (in.nRigidBody)
    detachAndWriteRandField(in, f, "rbParams");
  /* destroy any remaining fields */
  while(m->countFields())
    apf::destroyField( m->getField(0) );
  PHASTAIO_CLOSETIME(fclose(f);)
  double t1 = pcu::Time();
  if (!m->getPCU()->Self())
    lion_oprint(1,"solution written in %f seconds\n", t1 - t0);
}

} //end namespace ph
