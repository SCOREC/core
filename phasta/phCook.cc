#include <ph.h>
#include <chef.h>
#include <phstream.h>
#ifdef HAVE_SIMMETRIX
#include <phAttrib.h>
#include <apfSIM.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#endif
#include <phInput.h>
#include <phBC.h>
#include <phRestart.h>
#include <phAdapt.h>
#include <phOutput.h>
#include <phPartition.h>
#include <phFilterMatching.h>
#include "phInterfaceCutter.h"
#include "phiotimer.h" //for phastaio_initStats and phastaio_printStats
#include <parma.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <pcu_io.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <iostream>

static void print_stats(const char* name, double value, pcu::PCU *pcu_obj)
{
  double min, max, avg;
  min = pcu_obj->Min<double>(value);
  max = pcu_obj->Max<double>(value);
  avg = pcu_obj->Add<double>(value);
  avg /= pcu_obj->Peers();
  double imb = max / avg;
  if (!pcu_obj->Self())
    printf("%s: min %f max %f avg %f imb %f\n", name, min, max, avg, imb);
}

#define SIZET(a) static_cast<size_t>(a)

namespace {

pcu::PCU* createGroupComm(int splitFactor, pcu::PCU *PCUObj)
{
  int self = PCUObj->Self();
  int groupRank = self / splitFactor;
  int group = self % splitFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(PCUObj->GetMPIComm(), group, groupRank, &groupComm);
  return new pcu::PCU(groupComm);
}

void loadCommon(ph::Input& in, ph::BCs& bcs, gmi_model*& g, pcu::PCU *PCUObj)
{
  ph::loadModelAndBCs(in, g, bcs, PCUObj);
}

static apf::Mesh2* loadMesh(gmi_model*& g, ph::Input& in, pcu::PCU *PCUObj) {
  const char* meshfile = in.meshFileName.c_str();
  if (ph::mesh_has_ext(meshfile, "sms")) {
    if (in.simmetrixMesh == 0) {
      if (PCUObj->Self()==0)
        lion_eprint(1, "oops, turn on flag: simmetrixMesh\n");
      in.simmetrixMesh = 1;
      in.filterMatches = 0; //not support
    }
  }
  return ph::loadMesh(g, meshfile, PCUObj);
}

void originalMain(apf::Mesh2*& m, ph::Input& in,
    gmi_model* g, apf::Migration*& plan, pcu::PCU *pcu_obj)
{
  if(!m)
    m = loadMesh(g, in, pcu_obj);
  else
    apf::printStats(m);
// Need to set a flag to enable avoiding this when short on time  m->verify();
  if (in.useAttachedFields) 
     lion_eprint(1,"because useAttachedFields set restart not read\n");
  else if (in.solutionMigration && !in.useAttachedFields)
    ph::readAndAttachFields(in, m);
  else
    ph::attachZeroSolution(in, m);
  if (in.buildMapping)
    ph::buildMapping(m);
  apf::setMigrationLimit(SIZET(in.elementsPerMigration), pcu_obj);
  if (in.adaptFlag)
    ph::adapt(in, m);
  if (in.tetrahedronize)
    ph::tetrahedronize(in, m);
  if (in.simmetrixMesh == 0)
    plan = ph::split(in, m);
}

}//end namespace

namespace chef {
  static FILE* openfile_read(ph::Input&, const char* path, pcu::PCU *PCUObj) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = pcu_group_open(PCUObj->GetCHandle(), path, false);)
    return f;
  }

  static FILE* openfile_write(ph::Output& out, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = pcu_group_open(out.mesh->getPCU()->GetCHandle(), path, true);)
    return f;
  }

  static FILE* openstream_write(ph::Output& out, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = openGRStreamWrite(out.grs, path, out.mesh->getPCU());)
    return f;
  }

  static FILE* openstream_read(ph::Input& in, const char* path, pcu::PCU *PCUObj) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      PHASTAIO_OPENTIME(f = openRStreamRead(in.rs, PCUObj);)
    else {
      lion_eprint(1,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }
}

namespace ph {
  void checkBalance(apf::Mesh2* m, ph::Input& in) {
    /* check if balancing was requested */
      Parma_PrintPtnStats(m, "postSplit", false);
      if (in.prePhastaBalanceMethod != "none" && m->getPCU()->Peers() > 1)
        ph::balance(in,m);
  }

  void checkReorder(apf::Mesh2* m, ph::Input& in, int numMasters) {
    /* check if the mesh changed at all */
    if ( (m->getPCU()->Peers()!=numMasters) ||
        in.splitFactor > 1 ||
        in.adaptFlag ||
        in.prePhastaBalanceMethod != "none" ||
        in.tetrahedronize ||
        in.isReorder )
    {
      apf::MeshTag* order = NULL;

      print_stats("malloc used before Bfs", pcu::GetMem(), m->getPCU());

      if (in.isReorder && m->getPCU()->Peers() > 1)
        order = Parma_BfsReorder(m);

      print_stats("malloc used before reorder", pcu::GetMem(), m->getPCU());

      apf::reorderMdsMesh(m,order);

      print_stats("malloc used after reorder", pcu::GetMem(), m->getPCU());

    }
  }

  void balanceAndReorder(apf::Mesh2* m, ph::Input& in, int numMasters) {
    ph::checkBalance(m,in);
    ph::checkReorder(m,in,numMasters);
  }

  void preprocess(apf::Mesh2* m, Input& in, Output& out, BCs& bcs) {
    PCU_t h;
    h.ptr = static_cast<void*>(m->getPCU());
    phastaio_initStats(h);
    if(m->getPCU()->Peers() > 1)
      ph::migrateInterfaceItr(m, bcs);
    if (in.simmetrixMesh == 0)
      ph::checkReorder(m,in,m->getPCU()->Peers());
    if (in.adaptFlag)
      ph::goToStepDir(in.timeStepNumber, m->getPCU(), in.ramdisk);
    std::string path = ph::setupOutputDir(m->getPCU(), in.ramdisk);
    std::string subDirPath = path;
    ph::setupOutputSubdir(subDirPath, m->getPCU(), in.ramdisk);
    ph::enterFilteredMatching(m, in, bcs);
    ph::generateOutput(in, bcs, m, out);
    ph::exitFilteredMatching(m);
    // a path is not needed for inmem
    if ( in.writeRestartFiles ) {
      if(!m->getPCU()->Self()) lion_oprint(1,"write file-based restart file\n");
      // store the value of the function pointer
      FILE* (*fn)(Output& out, const char* path) = out.openfile_write;
      // set function pointer for file writing
      out.openfile_write = chef::openfile_write;
      ph::detachAndWriteSolution(in,out,m,subDirPath); //write restart
      // reset the function pointer to the original value
      out.openfile_write = fn;
    }
    else {
      ph::detachAndWriteSolution(in,out,m,subDirPath); //write restart
    }
    if ( ! in.outMeshFileName.empty() )
      m->writeNative(in.outMeshFileName.c_str());
    if ( in.writeGeomBCFiles ) {
      if(!m->getPCU()->Self()) lion_oprint(1,"write additional geomBC file for visualization\n");
      // store the value of the function pointer
      FILE* (*fn)(Output& out, const char* path) = out.openfile_write;
      // set function pointer for file writing
      out.openfile_write = chef::openfile_write;
      ph::writeGeomBC(out, path, in.timeStepNumber); //write geombc for viz only
      // reset the function pointer to the original value
      out.openfile_write = fn;
    }
    ph::writeGeomBC(out, subDirPath); //write geombc
    if(!m->getPCU()->Self())
      ph::writeAuxiliaryFiles(path, in.timeStepNumber, m->getPCU());
    m->verify();
#ifdef HAVE_SIMMETRIX
    gmi_model* g = m->getModel();
    ph::clearAttAssociation(g,in);
#endif
    if (in.adaptFlag)
      ph::goToParentDir();
    if(in.printIOtime) phastaio_printStats(h);
  }
  void preprocess(apf::Mesh2* m, Input& in, Output& out) {
    gmi_model* g = m->getModel();
    PCU_ALWAYS_ASSERT(g);
    BCs bcs;
    lion_eprint(1, "reading %s\n", in.attributeFileName.c_str());
    ph::readBCs(g, in.attributeFileName.c_str(), in.axisymmetry, bcs);
    if (!in.solutionMigration)
      ph::attachZeroSolution(in, m);
    if (in.buildMapping)
      ph::buildMapping(m);
    preprocess(m,in,out,bcs);
  }
}

namespace {
  struct GroupCode : public Parma_GroupCode {
    ph::Input* input;
    ph::BCs* boundary;
    apf::Mesh2* mesh;
    void run(int) {
      ph::Output groupOut;
      //streaming not supported from group code!
      groupOut.openfile_write = chef::openfile_write;
      ph::checkBalance(mesh,*input);
      ph::preprocess(mesh,*input,groupOut,*boundary);
    }
  };  
}

namespace chef {
  void bake(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& in, ph::Output& out, pcu::PCU *expandedPCUObj) {
    int shrinkFactor=0;
    if(in.splitFactor < 0) {
       shrinkFactor=-1*in.splitFactor; 
       in.splitFactor=1; // this is used in to set readers so if shrinking need to read all
    }
    PCU_ALWAYS_ASSERT(expandedPCUObj->Peers() % in.splitFactor == 0);
    apf::Migration* plan = 0;
    ph::BCs bcs;
    loadCommon(in, bcs, g, expandedPCUObj);
    const int worldRank = expandedPCUObj->Self();
    pcu::PCU *groupPCUObj = createGroupComm(in.splitFactor, expandedPCUObj);
    if ((worldRank % in.splitFactor) == 0)
      originalMain(m, in, g, plan, groupPCUObj);
    if (in.simmetrixMesh == 0)
      m = repeatMdsMesh(m, g, plan, in.splitFactor, expandedPCUObj);
    if (in.simmetrixMesh == 0 && shrinkFactor > 1){
      GroupCode code;
      apf::Unmodulo outMap(expandedPCUObj->Self(), expandedPCUObj->Peers());
      code.mesh=m;
      code.input=&in;
      code.boundary=&bcs;
      Parma_ShrinkPartition(code.mesh, shrinkFactor, code);
    } else {
      ph::checkBalance(m,in);
      ph::preprocess(m,in,out,bcs);
    }
  }
  void cook(gmi_model*& g, apf::Mesh2*& m, pcu::PCU *pcu_obj) {
    ph::Input in;
    in.load("adapt.inp", pcu_obj);
    in.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,in,out,pcu_obj);
    if ((in.writeVTK) == 1)  apf::writeVtkFiles("rendered",m);

  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, pcu::PCU *pcu_obj) {
    ctrl.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,ctrl,out,pcu_obj);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, GRStream* grs, pcu::PCU *pcu_obj) {
    ctrl.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openstream_write;
    out.grs = grs;
    bake(g,m,ctrl,out,pcu_obj);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* rs, pcu::PCU *pcu_obj) {
    ctrl.openfile_read = openstream_read;
    ctrl.rs = rs;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,ctrl,out,pcu_obj);
    return;
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* rs, GRStream* grs, pcu::PCU *pcu_obj) {
    ctrl.openfile_read = openstream_read;
    ctrl.rs = rs;
    ph::Output out;
    out.openfile_write = openstream_write;
    out.grs = grs;
    bake(g,m,ctrl,out,pcu_obj);
    return;
  }

  apf::Field* extractField(apf::Mesh* m,
    const char* packedFieldname,
    const char* requestFieldname,
    int firstComp,
    int numOfComp,
    bool simField) {
    return ph::extractField(m,packedFieldname,
             requestFieldname,firstComp,numOfComp,simField);
  }

  apf::Field* combineField(apf::Mesh* m,
    const char* packedFieldname,
    const char* inFieldname1,
    const char* inFieldname2,
    const char* inFieldname3) {
    return ph::combineField(m,packedFieldname,
             inFieldname1,inFieldname2,inFieldname3);
  }

  void readAndAttachFields(ph::Input& ctrl, apf::Mesh2*& m) {
    ph::readAndAttachFields(ctrl, m);
  }

  void balanceAndReorder(ph::Input& ctrl, apf::Mesh2* m) {
    ph::balanceAndReorder(m,ctrl,m->getPCU()->Peers());
  }

  void balance(ph::Input& ctrl, apf::Mesh2* m) {
    ph::checkBalance(m,ctrl);
  }

  void preprocess(apf::Mesh2*& m, ph::Input& in) {
    ph::Output out;
    out.openfile_write = openfile_write;
    ph::preprocess(m,in,out);
  }

  void preprocess(apf::Mesh2*& m, ph::Input& in, GRStream* grs) {
    ph::Output out;
    out.openfile_write = openstream_write;
    out.grs = grs;
    ph::preprocess(m,in,out);
  }
}