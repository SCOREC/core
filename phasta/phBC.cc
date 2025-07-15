#include <lionPrint.h>
#include "phBC.h"
#ifdef HAVE_SIMMETRIX
#include "phAttrib.h"
#include <gmi_sim.h>
#endif
#include "phAxisymmetry.h"
#include <apf.h>
#include <apfMesh.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <gmi.h>
#include <pcu_util.h>

namespace ph {

bool BC::operator<(const BC& other) const
{
  if (dim != other.dim)
    return dim < other.dim;
  return tag < other.tag;
}

ConstantBC::ConstantBC()
{
  value = 0;
}

ConstantBC::~ConstantBC()
{
  delete [] value;
}

double* ConstantBC::eval(apf::Vector3 const& x)
{
  (void)x;
  return value;
}

FieldBCs::~FieldBCs()
{
  while (!bcs.empty()) {
    Set::iterator it = bcs.begin();
    BC* bc = *it;
    bcs.erase(it);
    delete bc;
  }
}

static struct { const char* name; int size; } const knownSizes[7] =
{{"initial velocity", 3}
,{"comp1", 4}
,{"comp3", 4}
,{"comp1_elas", 4}
,{"comp3_elas", 4}
,{"traction vector", 3}
,{"traction vector melas", 3}
};

static int getSize(std::string const& name)
{
  for (int i = 0; i < 7; ++i)
    if (name == knownSizes[i].name)
      return knownSizes[i].size;
  return 1;
}

ConstantBC* makeConstantBC(BCs& bcs, std::string const& name, int dim, int tag,
    int nvals)
{
  if (!bcs.fields.count(name)) {
    FieldBCs fbcs;
    bcs.fields[name] = fbcs;
  }
  FieldBCs& fbcs = bcs.fields[name];
  ConstantBC* bc = new ConstantBC();
  bc->dim = dim;
  bc->tag = tag;
  bc->value = new double[nvals];
  fbcs.bcs.insert(bc);
  return bc;
}

static void readBC(std::string const& line, BCs& bcs)
{
  std::stringstream ss(line);
  std::string name;
  std::getline(ss, name, ':');
  int dim, tag;
  ss >> tag >> dim;
  int nvals = getSize(name);
  ConstantBC* bc = makeConstantBC(bcs, name, dim, tag, nvals);
  for (int i = 0; i < nvals; ++i)
    ss >> bc->value[i];
}

static void readBCsFromSPJ(const char* filename, BCs& bcs)
{
  std::ifstream file(filename);
  PCU_ALWAYS_ASSERT(file.is_open()); //check if the spj file could be opened successfully
  std::string line;
  while (std::getline(file, line, '\n')) {
    if (line[0] == '#')
      continue;
    readBC(line, bcs);
  }
}

void readBCs(gmi_model* m, const char* attFile, bool axisymmetry, BCs& bcs)
{
  /* now load attributes.
     only two cases:
     either its an SPJ file or they came in with the model */
  if (gmi_has_ext(attFile, "spj"))
    readBCsFromSPJ(attFile, bcs);
#ifdef HAVE_SIMMETRIX
  else
    getSimmetrixAttributes(m, bcs);
#endif
  if (axisymmetry)
    attachAllAngleBCs(m, bcs);
}

void loadModelAndBCs(ph::Input& in, gmi_model*& m, BCs& bcs, pcu::PCU *PCUObj)
{
  double t0 = pcu::Time();
  const char* modelfile = in.modelFileName.c_str();
  const char* attribfile = in.attributeFileName.c_str();
  /* loading the model */
  /* case 1: meshmodel */
  if (gmi_has_ext(modelfile, "dmg"))
    m = gmi_load(modelfile);
#ifdef HAVE_SIMMETRIX
  /* cases 2: Simmetrix model (and possibly attributes) file */
  else if (gmi_has_ext(modelfile, "smd"))
    m = gmi_sim_load(0, modelfile);
  /* cases 3&4: assuming native model file */
  else {
    /* case 3: native model and Simmetrix attributes file */
    if (gmi_has_ext(attribfile, "smd"))
      m = gmi_sim_load(modelfile, attribfile);
    /* case 4: just a native model */
    else
      m = gmi_sim_load(modelfile, 0);
  }
#endif
  /* load attributes */
  readBCs(m, attribfile, in.axisymmetry, bcs);
  double t1 = pcu::Time();
  if (!PCUObj->Self())
    lion_oprint(1,"\"%s\" and \"%s\" loaded in %f seconds\n", modelfile, attribfile, t1 - t0);
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

static KnownBC const naturalBCs[11] = {
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
  {"traction vector melas",  10, 9, applyVector},
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

bool haveBC(BCs& bcs, std::string const& name)
{
  return bcs.fields.count(name);
}

double* getBCValue(gmi_model* gm, FieldBCs& bcs, gmi_ent* ge,
    apf::Vector3 const& x)
{
  ConstantBC keyObj;
  keyObj.tag = gmi_tag(gm, ge);
  keyObj.dim = gmi_dim(gm, ge);
  BC* key = &keyObj; /* holding the PGI compiler's hand */
  FieldBCs::Set::iterator it = bcs.bcs.find(key);
  if (it == bcs.bcs.end())
    return 0;
  BC* bc = *it;
  return bc->eval(x);
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
    apf::Vector3 const& x,
    KnownBC const& kbc,
    double* values, int* bits)
{
  double* v = getBCValue(gm, bcs, ge, x);
  if (v) {
    kbc.apply(values, bits, kbc, v);
    return true;
  }
  bool appliedAny = false;
  gmi_set* up = gmi_adjacent(gm, ge, gmi_dim(gm, ge) + 1);
  for (int i = 0; i < up->n; ++i)
    if ( applyBC(gm, up->e[i], bcs, x, kbc, values, bits) )
      appliedAny = true;
  gmi_free_set(up);
  return appliedAny;
}

static bool applyBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs,
    apf::Vector3 const& x,
    KnownBC const* knownBCs,
    int nKnownBCs,
    double* values, int* bits)
{
  bool appliedAny = false;
  for (int i = 0; i < nKnownBCs; ++i) {
    std::string s(knownBCs[i].name);
    if ( ! haveBC(appliedBCs, s))
      continue;
    FieldBCs& fbcs = appliedBCs.fields[s];
    if ( applyBC(gm, ge, fbcs, x, knownBCs[i], values, bits) )
      appliedAny = true;
  }
  return appliedAny;
}

bool applyNaturalBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values, int* bits)
{
  return applyBCs(gm, ge, appliedBCs, x, naturalBCs,
      sizeof(naturalBCs) / sizeof(KnownBC), values, bits);
}

bool applyEssentialBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values, int* bits)
{
  bool didSimple = applyBCs(gm, ge, appliedBCs, x, essentialBCs,
      sizeof(essentialBCs) / sizeof(KnownBC), values, bits);
  /* the complexity of velocity constraints is delegated to
     the code in phConstraint.cc */
  bool didVelocity = applyVelocityConstaints(gm, appliedBCs,
      ge, x, values, bits);
  bool didElastic = applyElasticConstaints(gm, appliedBCs,
      ge, x, values, bits);
  return didSimple || didVelocity || didElastic;
}

bool applySolutionBCs(gmi_model* gm, gmi_ent* ge,
    BCs& appliedBCs, apf::Vector3 const& x, double* values)
{
  return applyBCs(gm, ge, appliedBCs, x, solutionBCs,
      sizeof(solutionBCs) / sizeof(KnownBC), values, 0);
}

}
