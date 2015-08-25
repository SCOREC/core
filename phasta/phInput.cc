#include "phInput.h"
#include <fstream>
#include <map>
#include "ph.h"
#include <cassert>

namespace ph {

static void setDefaults(Input& in)
{
  in.outMeshFileName = "";
  in.numSplit = 10;
  in.tetrahedronize = 0;
  in.localPtn = true; 
  in.recursivePtn = -1;
  in.recursiveUR = 1;
  in.parmaPtn = 0; // No Parma by default
  in.displacementMigration = 0; // Do not migrate displacement field by default
  in.dwalMigration = 0; // Do not migrate dwal field by default
  in.buildMapping = 0; // Do not build the mapping field by default
  in.elementsPerMigration = 1000*1000; // 100k elms per round
  in.threaded = 1;
  in.initBubbles = 0;
  in.formElementGraph = 0;
  in.restartFileName = "restart";
  in.phastaIO = 1;
  in.snap = 0;
  in.splitAllLayerEdges = 0;
  in.filterMatches = 0;
  in.axisymmetry = 0;
  in.elementImbalance = 1.03;
  in.vertexImbalance = 1.05;
}

Input::Input()
{
  setDefaults(*this);
}

typedef std::map<std::string, std::string*> StringMap;
typedef std::map<std::string, int*> IntMap;
typedef std::map<std::string, double*> DblMap;

static void formMaps(Input& in, StringMap& stringMap, IntMap& intMap, DblMap& dblMap)
{
  intMap["globalP"] = &in.globalP;
  intMap["timeStepNumber"] = &in.timeStepNumber;
  intMap["ensa_dof"] = &in.ensa_dof;
  stringMap["restartFileName"] = &in.restartFileName;
  stringMap["attributeFileName"] = &in.attributeFileName;
  stringMap["meshFileName"] = &in.meshFileName;
  stringMap["outMeshFileName"] = &in.outMeshFileName;
  stringMap["modelFileName"] = &in.modelFileName;
  stringMap["outputFormat"] = &in.outputFormat;
  stringMap["partitionMethod"] = &in.partitionMethod;
  intMap["adaptFlag"] = &in.adaptFlag;
  intMap["rRead"] = &in.rRead;
  intMap["rStart"] = &in.rStart;
  intMap["AdaptStrategy"] = &in.adaptStrategy;
  intMap["RecursiveUR"] = &in.recursiveUR;
  intMap["Periodic"] = &in.periodic;
  intMap["prCD"] = &in.prCD;
  intMap["timing"] = &in.timing;
  intMap["internalBCNodes"] = &in.internalBCNodes;
  intMap["WRITEASC"] = &in.writeDebugFiles;
  intMap["phastaIO"] = &in.phastaIO;
  intMap["splitFactor"] = &in.splitFactor;
  intMap["SolutionMigration"] = &in.solutionMigration;
  intMap["DisplacementMigration"] = &in.displacementMigration;
  intMap["isReorder"] = &in.isReorder;
  intMap["numSplit"] = &in.numSplit;
  intMap["Tetrahedronize"] = &in.tetrahedronize;
  intMap["LocalPtn"] = &in.localPtn;
  intMap["RecursivePtn"] = &in.recursivePtn;
  intMap["ParmaPtn"] = &in.parmaPtn;
  intMap["dwalMigration"] = &in.dwalMigration;
  intMap["buildMapping"] = &in.buildMapping;
  intMap["elementsPerMigration"] = &in.elementsPerMigration;
  intMap["threaded"] = &in.threaded;
  intMap["initBubbles"] = &in.initBubbles;
  intMap["formElementGraph"] = &in.formElementGraph;
  intMap["snap"] = &in.snap;
  intMap["splitAllLayerEdges"] = &in.splitAllLayerEdges;
  intMap["filterMatches"] = &in.filterMatches;
  intMap["axisymmetry"] = &in.axisymmetry;
  dblMap["elementImbalance"] = &in.elementImbalance;
  dblMap["vertexImbalance"] = &in.vertexImbalance;
}

template <class T>
static bool tryReading(std::string const& name,
    std::ifstream& f,
    std::map<std::string, T*>& map)
{
  typename std::map<std::string, T*>::iterator it = map.find(name);
  if (it == map.end())
    return false;
  f >> *(it->second);
  return true;
}

static void readInputFile(
    Input& in,
    const char* filename,
    StringMap& stringMap,
    IntMap& intMap,
    DblMap& dblMap)
{
  std::ifstream f(filename);
  if (!f)
    fail("could not open \"%s\"", filename);
  std::string name;
  while (f >> name) {
    if (name[0] == '#') {
      std::getline(f, name, '\n');
      continue;
    }
    if (tryReading(name, f, stringMap))
      continue;
    if (tryReading(name, f, intMap))
      continue;
    if (tryReading(name, f, dblMap))
      continue;
    /* the WEIRD parameter ! */
    if (name == "RecursivePtnStep") {
      if (in.recursivePtn == -1)
        fail("RecursivePtn needs to be set before RecursivePtnStep");
      in.recursivePtnStep.allocate(in.recursivePtn);
      for (int i = 0; i < in.recursivePtn; ++i)
        f >> in.recursivePtnStep[i];
      continue;
    }
    fail("unknown variable \"%s\" in %s\n", name.c_str(), filename);
  }
}

static void validate(Input& in)
{
  assert(in.parmaPtn == 0 || in.parmaPtn == 1);
  assert(in.elementImbalance > 1.0 && in.elementImbalance <= 2.0);
  assert(in.vertexImbalance > 1.0 && in.vertexImbalance <= 2.0);
  assert( ! (in.buildMapping && in.adaptFlag));
}

void Input::load(const char* filename)
{
  setDefaults(*this);
  StringMap stringMap;
  IntMap intMap;
  DblMap dblMap;
  formMaps(*this, stringMap, intMap, dblMap);
  readInputFile(*this, filename, stringMap, intMap, dblMap);
  validate(*this);
}

/* see the comments in phBC.h and
   phOutput.h for explanations about
   these count*BCs functions */
int countNaturalBCs(Input& in)
{
  return in.ensa_dof + 1;
}

int countEssentialBCs(Input& in)
{
  return in.ensa_dof + 7;
}

int countScalarBCs(Input& in)
{
  return in.ensa_dof - 5;
}

}
