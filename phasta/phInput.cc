#include <PCU.h>
#include "phInput.h"
#include <fstream>
#include <map>
#include <set>
#include "ph.h"
#include <pcu_util.h>

/** \file phInput.cc
    \brief The implementation of Chef's interface for execution control */

namespace ph {

static void setDefaults(Input& in)
{
  in.timeStepNumber = 0;
  in.ensa_dof = 0;
  in.ensa_melas_dof = 0; 
  in.outMeshFileName = "";
  in.adaptFlag = 0;
  in.rRead = 0;
  in.rStart = 0;
  in.preAdaptBalanceMethod = "parma";
  in.midAdaptBalanceMethod = "zoltan";
  in.postAdaptBalanceMethod = "zoltan";
  in.prePhastaBalanceMethod = "parma-gap";
  in.adaptStrategy = -1;
  in.adaptErrorThreshold = 1e-6;  //used by adaptStrategy=2 (runFromErrorThreshold)
  in.adaptErrorFieldName = "errors"; //used by adaptStrategy=2 (runFromErrorThreshold)
  in.adaptErrorFieldIndex = 5; //used by adaptStrategy=2 (runFromErrorThreshold)
  in.periodic = 0;
  in.prCD = -1;
  in.timing = 0;
  in.internalBCNodes = 0;
  in.writeDebugFiles = 0;
  in.splitFactor = 1;
  in.partitionMethod = "rib";
  in.localPtn = 1;
  in.solutionMigration = 1;
  in.useAttachedFields = 0;
  in.isReorder = 0;
  in.openfile_read = 0;
  in.tetrahedronize = 0;
  in.recursiveUR = 1;
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
  in.transferParametric = 0;
  in.splitAllLayerEdges = 0;
  in.filterMatches = 0;
  in.axisymmetry = 0;
  in.parmaLoops = 3; //a magical value
  in.parmaVerbosity = 1; //fairly quiet
  in.writeGeomBCFiles = 1;
  in.ramdisk = 0;
  in.meshqCrtn = 0.027; 
  in.elementImbalance = 1.03;
  in.vertexImbalance = 1.05;
  in.rs = 0;
  in.formEdges = 0;
  in.simmetrixMesh = 0;
  in.maxAdaptIterations = 3;
  in.adaptShrinkLimit = 10000;
  in.validQuality = 1.0e-10;
  in.printIOtime = 0;
  in.mesh2geom = 0;
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
  intMap["timeStepNumber"] = &in.timeStepNumber;
  intMap["ensa_dof"] = &in.ensa_dof;
  intMap["ensa_melas_dof"] = &in.ensa_melas_dof;
  stringMap["restartFileName"] = &in.restartFileName;
  stringMap["attributeFileName"] = &in.attributeFileName;
  stringMap["meshFileName"] = &in.meshFileName;
  stringMap["outMeshFileName"] = &in.outMeshFileName;
  stringMap["modelFileName"] = &in.modelFileName;
  stringMap["outputFormat"] = &in.outputFormat;
  stringMap["partitionMethod"] = &in.partitionMethod;
  stringMap["preAdaptBalanceMethod"] = &in.preAdaptBalanceMethod;
  stringMap["midAdaptBalanceMethod"] = &in.midAdaptBalanceMethod;
  stringMap["postAdaptBalanceMethod"] = &in.postAdaptBalanceMethod;
  stringMap["prePhastaBalanceMethod"] = &in.prePhastaBalanceMethod;
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
  intMap["UseAttachedFields"] = &in.useAttachedFields;
  intMap["DisplacementMigration"] = &in.displacementMigration;
  intMap["isReorder"] = &in.isReorder;
  intMap["Tetrahedronize"] = &in.tetrahedronize;
  intMap["LocalPtn"] = &in.localPtn;
  intMap["dwalMigration"] = &in.dwalMigration;
  intMap["buildMapping"] = &in.buildMapping;
  intMap["elementsPerMigration"] = &in.elementsPerMigration;
  intMap["threaded"] = &in.threaded;
  intMap["initBubbles"] = &in.initBubbles;
  intMap["formElementGraph"] = &in.formElementGraph;
  intMap["snap"] = &in.snap;
  intMap["transferParametric"] = &in.transferParametric;
  intMap["splitAllLayerEdges"] = &in.splitAllLayerEdges;
  intMap["filterMatches"] = &in.filterMatches;
  intMap["axisymmetry"] = &in.axisymmetry;
  intMap["parmaLoops"] = &in.parmaLoops;
  intMap["parmaVerbosity"] = &in.parmaVerbosity;
  intMap["writeGeomBCFiles"] = &in.writeGeomBCFiles;
  intMap["ramdisk"] = &in.ramdisk;
  dblMap["validQuality"] = &in.validQuality;
  dblMap["meshqCrtn"] = &in.meshqCrtn;
  dblMap["elementImbalance"] = &in.elementImbalance;
  dblMap["vertexImbalance"] = &in.vertexImbalance;
  dblMap["adaptShrinkLimit"] = &in.adaptShrinkLimit;
  intMap["formEdges"] = &in.formEdges;
  intMap["simmetrixMesh"] = &in.simmetrixMesh;
  intMap["maxAdaptIterations"] = &in.maxAdaptIterations;
  intMap["printIOtime"] = &in.printIOtime;
  intMap["mesh2geom"] = &in.mesh2geom;
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

typedef std::set<std::string> stringset;

static void makeDeprecated(stringset& old)
{
  old.insert("globalP");
  old.insert("numSplit");
  old.insert("ParmaPtn");
  old.insert("RecursivePtn");
  old.insert("RecursivePtnStep");
  old.insert("writePhastaFiles");
}

static bool deprecated(stringset& old, std::string const& name)
{
  if( old.count(name) ) {
    if( !PCU_Comm_Self() )
      fprintf(stderr, "WARNING deprecated input \"%s\" ... "
          "carefully check stderr and stdout for unexpected behavior\n",
          name.c_str());
    return true;
  } else {
    return false;
  }
}

static void readInputFile(
    const char* filename,
    StringMap& stringMap,
    IntMap& intMap,
    DblMap& dblMap)
{
  stringset old;
  makeDeprecated(old);
  std::ifstream f(filename);
  if (!f)
    fail("could not open \"%s\"", filename);
  std::string name;
  while (f >> name) {
    if (name[0] == '#' || deprecated(old,name)) {
      std::getline(f, name, '\n');
      continue;
    }
    if (tryReading(name, f, stringMap))
      continue;
    if (tryReading(name, f, intMap))
      continue;
    if (tryReading(name, f, dblMap))
      continue;
    fail("unknown variable \"%s\" in %s\n", name.c_str(), filename);
  }
}

static void validate(Input& in)
{
  PCU_ALWAYS_ASSERT(in.elementImbalance > 1.0 && in.elementImbalance <= 2.0);
  PCU_ALWAYS_ASSERT(in.vertexImbalance > 1.0 && in.vertexImbalance <= 2.0);
  PCU_ALWAYS_ASSERT( ! (in.buildMapping && in.adaptFlag));
}

void Input::load(const char* filename)
{
  setDefaults(*this);
  StringMap stringMap;
  IntMap intMap;
  DblMap dblMap;
  formMaps(*this, stringMap, intMap, dblMap);
  readInputFile(filename, stringMap, intMap, dblMap);
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
  if(!in.ensa_melas_dof)
    return in.ensa_dof + 7;
  else
    return 3+2+4+7+8; // (assuming 4 scalars to be ON) and 8 is for ec11 ec12 ec13 em1 ec21 ec22 ec23 em2
}

int countScalarBCs(Input& in)
{
  return in.ensa_dof - 5;
}

}
