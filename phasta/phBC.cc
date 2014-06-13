#include "phBC.h"
#include <apf.h>
#include <apfMesh.h>
#include <fstream>
#include <sstream>

namespace ph {

BC::BC()
{
  values = 0;
}

BC::~BC()
{
  delete [] values;
}

bool BC::operator<(const BC& other) const
{
  if (dim != other.dim)
    return dim < other.dim;
  return tag < other.tag;
}

static struct { const char* name; int size; } const knownSizes[4] =
{{"initial velocity", 3}
,{"comp1", 4}
,{"comp3", 4}
,{"traction vector", 3}
};

static int getSize(std::string const& name)
{
  for (int i = 0; i < 4; ++i)
    if (name == knownSizes[i].name)
      return knownSizes[i].size;
  return 1;
}

static void readBC(std::string const& line, BCs& bcs)
{
  std::stringstream ss(line);
  std::string name;
  std::getline(ss, name, ':');
  if (!bcs.fields.count(name)) {
    FieldBCs fbcs;
    fbcs.size = getSize(name);
    bcs.fields[name] = fbcs;
  }
  FieldBCs& fbcs = bcs.fields[name];
  BC bc;
  ss >> bc.tag >> bc.dim;
  bc.values = new double[fbcs.size];
  for (int i = 0; i < fbcs.size; ++i)
    ss >> bc.values[i];
  fbcs.bcs.insert(bc);
  bc.values = 0; //ownership of pointer transferred, prevent delete from here
}

void readBCs(const char* filename, BCs& bcs)
{
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line, '\n')) {
    if (line[0] == '#')
      continue;
    readBC(line, bcs);
  }
}

struct KnownBC
{
  const char* name;
  int offset;
  int bit;
  void (*apply)(double* values, int* bits,
    struct KnownBC const& bc, double* inval);
};

static void applyScalar(double* values, int* bits,
    double value, int offset, int bit)
{
  if (offset != -1)
    values[offset] = value;
  if (bit != -1)
    *bits |= (1<<bit);
}

static void applyScalar(double* outval, int* bits,
    KnownBC const& bc, double* inval)
{
  applyScalar(outval, bits, *inval, bc.offset, bc.bit);
}

static void applyMagnitude(double* values, KnownBC const& bc, double v)
{
  values[bc.offset + 3] = v;
}

static void applyComp1(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  int best = 0;
  for (int i = 1; i < 3; ++i)
    if (fabs(inval[i]) > fabs(inval[best]))
      best = i;
  apf::Vector3 v(inval);
  double mag = v.getLength();
  v = v / mag;
  applyScalar(values, bits, v[best], bc.offset + best, bc.bit + best);
  applyMagnitude(values, bc, mag);
}

static void applyComp3(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  apf::Vector3 v(inval);
  double mag = v.getLength();
  v = v / mag;
  for (int i = 0; i < 3; ++i)
    applyScalar(values, bits, v[i], bc.offset + i, bc.bit + i);
  applyMagnitude(values, bc, mag);
}

static void applyVector(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  for (int i = 0; i < 3; ++i)
    values[bc.offset + i] = inval[i];
  *bits |= (1<<bc.bit);
}

static void applySurfID(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  bits[1] = *inval;
}

static KnownBC const essentialBCs[9] = {
  {"density",          0, 0, applyScalar},
  {"temperature",      1, 1, applyScalar},
  {"pressure",         2, 2, applyScalar},
  {"comp 1",           3, 3, applyComp1},
  {"comp 3",           3, 3, applyComp3},
  {"scalar_1",        12, 6, applyScalar},
  {"scalar_2",        13, 7, applyScalar},
  {"scalar_3",        14, 8, applyScalar},
  {"scalar_4",        15, 9, applyScalar},
};

static KnownBC const naturalBCs[10] = {
  {"mass flux",        0, 0, applyScalar},
  {"natural pressure", 1, 1, applyScalar},
  {"traction vector",  2, 2, applyVector},
  {"heat flux",        5, 3, applyScalar},
  {"turbulence wall", -1, 4, applyScalar},
  {"scalar_1 flux",    6, 5, applyScalar},
  {"scalar_2 flux",    7, 6, applyScalar},
  {"scalar_3 flux",    8, 7, applyScalar},
  {"scalar_4 flux",    9, 8, applyScalar},
  {"surf ID",         -1,-1, applySurfID},
};

double* checkForBC(int dim, int tag, BCs& bcs, KnownBC const& kbc)
{
  std::string name(kbc.name);
  if (!bcs.fields.count(name))
    return 0;
  FieldBCs& fbcs = bcs.fields[name];
  BC key;
  key.tag = tag;
  key.dim = dim;
  FieldBCs::Set::iterator it = fbcs.bcs.find(key);
  if (it == fbcs.bcs.end())
    return 0;
  BC& bc = const_cast<BC&>(*it);
  return bc.values;
}

void applyBCs(apf::Mesh* m, apf::MeshEntity* e,
    BCs& appliedBCs,
    KnownBC const* knownBCs,
    int nKnownBCs,
    double* values, int* bits)
{
  apf::ModelEntity* me = m->toModel(e);
  int md = m->getModelType(me);
  int mt = m->getModelTag(me);
  for (int i = 0; i < nKnownBCs; ++i) {
    double* bcvalues = checkForBC(md, mt, appliedBCs, knownBCs[i]);
    if (!bcvalues)
      continue;
    knownBCs[i].apply(values, bits, knownBCs[i], bcvalues);
  }
}

void applyNaturalBCs(apf::Mesh* m, apf::MeshEntity* f,
    BCs& appliedBCs,
    double* values, int* bits)
{
  applyBCs(m, f, appliedBCs, naturalBCs,
      sizeof(naturalBCs) / sizeof(KnownBC), values, bits);
}

void applyEssentialBCs(apf::Mesh* m, apf::MeshEntity* v,
    BCs& appliedBCs,
    double* values, int* bits)
{
  applyBCs(m, v, appliedBCs, essentialBCs,
      sizeof(essentialBCs) / sizeof(KnownBC), values, bits);
}

}
