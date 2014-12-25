#include <PCU.h>
#include "phBC.h"
#include <apf.h>
#include <apfMesh.h>
#include <fstream>
#include <sstream>
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
  double t0 = PCU_Time();
  std::ifstream file(filename);
  assert(file.is_open()); //check if the spj file could be opened successfully
  std::string line;
  while (std::getline(file, line, '\n')) {
    if (line[0] == '#')
      continue;
    readBC(line, bcs);
  }
  double t1 = PCU_Time();
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

static bool applyBit2(int* bits, int bit)
{
  if (bit == -1)
    return true;
  int mask = (1<<bit);
  if (*bits & mask)
    return false;
  *bits |= mask;
  return true;
}

static void applyBit(double*, int* bits,
    KnownBC const& bc, double*)
{
  applyBit2(bits, bc.bit);
}

static void applyScalar2(double* values, int* bits,
    double value, int offset, int bit)
{
  bool alreadySet = ! applyBit2(bits, bit);
  /* only zero values override existing values */
  if (( ! alreadySet) || ( ! value))
    values[offset] = value;
}

static void applyScalar(double* outval, int* bits,
    KnownBC const& bc, double* inval)
{
  applyScalar2(outval, bits, *inval, bc.offset, bc.bit);
}

static void applyVector2(double* values, int* bits,
    double* inval, int offset, int bit)
{
  bool alreadySet = ! applyBit2(bits, bit);
  if (alreadySet)
    /* nonzero vectors cannot override */
    if (inval[0] || inval[1] || inval[2])
      return;
  for (int i = 0; i < 3; ++i)
    values[offset + i] = inval[i];
}

static void applyVector(double* outval, int* bits,
    KnownBC const& bc, double* inval)
{
  applyVector2(outval, bits, inval, bc.offset, bc.bit);
}

static void applySurfID(double*, int* bits,
    KnownBC const&, double* inval)
{
  bits[1] = *inval;
}

static KnownBC const essentialBCs[7] = {
  {"density",          0, 0, applyScalar},
  {"temperature",      1, 1, applyScalar},
  {"pressure",         2, 2, applyScalar},
/*  these conditions are too complex to be dealt with
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
  {"turbulence wall", -1, 4, applyBit},
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
   try to find attribute (kbc) attached to
   geometric entities by searching all upward
   adjacencies.
   we use depth first search starting from
   the model entity of origin, using upward
   adjacencies as graph edges.
   if an attached attribute is found on an
   entity, the search continues without looking
   at its upward adjacencies */
static bool applyBC(gmi_model* gm, gmi_ent* ge,
    FieldBCs& bcs,
    KnownBC const& kbc,
    double* values, int* bits)
{
  double* v = getValuesOn(gm, bcs, ge);
  if (v) {
    kbc.apply(values, bits, kbc, v);
    return true;
  }
  bool appliedAny = false;
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    if ( applyBC(gm, up->e[i], bcs, kbc, values, bits) )
      appliedAny = true;
  gmi_free_set(up);
  return appliedAny;
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
    if ( applyBC(gm, ge, fbcs, knownBCs[i], values, bits) )
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
