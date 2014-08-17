#include "phBC.h"
#include <apf.h>
#include <apfMesh.h>
#include <fstream>
#include <sstream>
#include <PCU.h>
#include <gmi.h>

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
  double t0 = MPI_Wtime();
  std::ifstream file(filename);
  std::string line;
  while (std::getline(file, line, '\n')) {
    if (line[0] == '#')
      continue;
    readBC(line, bcs);
  }
  double t1 = MPI_Wtime();
  if (!PCU_Comm_Self())
    printf("\"%s\" loaded in %f seconds\n", filename, t1 - t0);
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

static void applyVector(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  for (int i = 0; i < 3; ++i)
    values[bc.offset + i] = inval[i];
  if (bc.bit != -1)
    *bits |= (1<<bc.bit);
}

static void applySurfID(double* values, int* bits,
    KnownBC const& bc, double* inval)
{
  bits[1] = *inval;
}

static KnownBC const essentialBCs[7] = {
  {"density",          0, 0, applyScalar},
  {"temperature",      1, 1, applyScalar},
  {"pressure",         2, 2, applyScalar},
/*  these guys are too complex to be dealt with
    here, see phConstraint.cc */
//{"comp1",           3, 3, applyComp1},
//{"comp3",           3, 3, applyComp3},
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

static KnownBC const solutionBCs[7] = {
  {"initial pressure",         0,-1, applyScalar},
  {"initial velocity",         1,-1, applyVector},
  {"initial temperature",      4,-1, applyScalar},
  {"initial scalar_1",         5,-1, applyScalar},
  {"initial scalar_2",         6,-1, applyScalar},
  {"initial scalar_3",         7,-1, applyScalar},
  {"initial scalar_4",         8,-1, applyScalar},
};

bool hasBC(BCs& bcs, std::string const& name)
{
  return bcs.fields.count(name);
}

double* getValuesOn(gmi_model* gm, FieldBCs& bcs, gmi_ent* ge)
{
  BC key;
  key.tag = gmi_tag(gm, ge);
  key.dim = gmi_dim(gm, ge);
  FieldBCs::Set::iterator it = bcs.bcs.find(key);
  if (it == bcs.bcs.end())
    return 0;
  BC& bc = const_cast<BC&>(*it);
  return bc.values;
}

/* starting from the current geometric entity,
   try to find an attribute (kbc) attached to
   a geometric entity by searching all upward
   adjacencies.
ex: for a mesh vertex classified on a model
    vertex, the model vertex is first checked
    for the attribute, then the model edges
    adjacent to the model vertex, then the model
    faces adjacent to those model edges.
   this is done by the following recursive function. */
static double* getFirstApplied(gmi_model* gm, gmi_ent* ge,
    FieldBCs& bcs, KnownBC const& kbc)
{
  double* v = getValuesOn(gm, bcs, ge);
  if (v)
    return v;
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i) {
    v = getFirstApplied(gm, up->e[i], bcs, kbc);
    if (v)
      break;
  }
  gmi_free_set(up);
  return v;
}

static bool applyBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs,
    KnownBC const* knownBCs,
    int nKnownBCs,
    double* values, int* bits)
{
  bool appliedAny = false;
  for (int i = 0; i < nKnownBCs; ++i) {
    std::string s(knownBCs[i].name);
    if ( ! hasBC(appliedBCs, s))
      continue;
    FieldBCs& fbcs = appliedBCs.fields[s];
    double* bcvalues = getFirstApplied(gm, ge, fbcs, knownBCs[i]);
    if (!bcvalues)
      continue;
    knownBCs[i].apply(values, bits, knownBCs[i], bcvalues);
    appliedAny = true;
  }
  return appliedAny;
}

bool applyNaturalBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, double* values, int* bits)
{
  return applyBCs(gm, ge, appliedBCs, naturalBCs,
      sizeof(naturalBCs) / sizeof(KnownBC), values, bits);
}

bool applyEssentialBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, double* values, int* bits)
{
  bool didSimple = applyBCs(gm, ge, appliedBCs, essentialBCs,
      sizeof(essentialBCs) / sizeof(KnownBC), values, bits);
  /* the complexity of velocity constraints is delegated to
     the code in phConstraint.cc */
  bool didVelocity = applyVelocityConstaints(gm, appliedBCs,
      ge, values, bits);
  return didSimple || didVelocity;
}

bool applySolutionBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, double* values)
{
  return applyBCs(gm, ge, appliedBCs, solutionBCs,
      sizeof(solutionBCs) / sizeof(KnownBC), values, 0);
}

}
