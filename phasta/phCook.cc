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
#include <PCU.h>
#include <pcu_io.h>
#include <pcu_util.h>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <iostream>

#define SIZET(a) static_cast<size_t>(a)

namespace {

void switchToMasters(int splitFactor)
{
  int self = PCU_Comm_Self();
  int groupRank = self / splitFactor;
  int group = self % splitFactor;
  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

void loadCommon(ph::Input& in, ph::BCs& bcs, gmi_model*& g)
{
  ph::loadModelAndBCs(in, g, bcs);
}

static apf::Mesh2* loadMesh(gmi_model*& g, ph::Input& in) {
  const char* meshfile = in.meshFileName.c_str();
  if (ph::mesh_has_ext(meshfile, "sms")) {
    if (in.simmetrixMesh == 0) {
      if (PCU_Comm_Self()==0)
        fprintf(stderr, "oops, turn on flag: simmetrixMesh\n");
      in.simmetrixMesh = 1;
      in.filterMatches = 0; //not support
    }
  }
  return ph::loadMesh(g, meshfile);
}

void originalMain(apf::Mesh2*& m, ph::Input& in,
    gmi_model* g, apf::Migration*& plan)
{
  if(!m)
    m = loadMesh(g, in);
  else
    apf::printStats(m);
  m->verify();
  if (in.solutionMigration && !in.useAttachedFields)
    ph::readAndAttachFields(in, m);
  else
    ph::attachZeroSolution(in, m);
  if (in.buildMapping)
    ph::buildMapping(m);
  apf::setMigrationLimit(SIZET(in.elementsPerMigration));
  if (in.adaptFlag)
    ph::adapt(in, m);
  if (in.tetrahedronize)
    ph::tetrahedronize(in, m);
  if (in.simmetrixMesh == 0)
    plan = ph::split(in, m);
}

}//end namespace

namespace chef {
  static FILE* openfile_read(ph::Input&, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = pcu_group_open(path, false);)
    return f;
  }

  static FILE* openfile_write(ph::Output&, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = pcu_group_open(path, true);)
    return f;
  }

  static FILE* openstream_write(ph::Output& out, const char* path) {
    FILE* f = NULL;
    PHASTAIO_OPENTIME(f = openGRStreamWrite(out.grs, path);)
    return f;
  }

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      PHASTAIO_OPENTIME(f = openRStreamRead(in.rs);)
    else {
      fprintf(stderr,
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
      if (in.prePhastaBalanceMethod != "none" && PCU_Comm_Peers() > 1)
        ph::balance(in,m);
  }

  void checkReorder(apf::Mesh2* m, ph::Input& in, int numMasters) {
    /* check if the mesh changed at all */
    if ( (PCU_Comm_Peers()!=numMasters) ||
        in.adaptFlag ||
        in.prePhastaBalanceMethod != "none" ||
        in.tetrahedronize ||
        in.isReorder )
    {
      apf::MeshTag* order = NULL;
      if (in.isReorder && PCU_Comm_Peers() > 1)
        order = Parma_BfsReorder(m);
      apf::reorderMdsMesh(m,order);
    }
  }

  void balanceAndReorder(apf::Mesh2* m, ph::Input& in, int numMasters) {
    ph::checkBalance(m,in);
    ph::checkReorder(m,in,numMasters);
  }

  void preprocess(apf::Mesh2* m, Input& in, Output& out, BCs& bcs) {
    phastaio_initStats();
    if(PCU_Comm_Peers() > 1)
      ph::migrateInterfaceItr(m, bcs);
    if (in.simmetrixMesh == 0)
      ph::checkReorder(m,in,PCU_Comm_Peers());
    if (in.adaptFlag)
      ph::goToStepDir(in.timeStepNumber,in.ramdisk);
    std::string path = ph::setupOutputDir(in.ramdisk);
    std::string subDirPath = path;
    ph::setupOutputSubdir(subDirPath,in.ramdisk);
    ph::enterFilteredMatching(m, in, bcs);
    ph::generateOutput(in, bcs, m, out);
    ph::exitFilteredMatching(m);
    // a path is not needed for inmem
    ph::detachAndWriteSolution(in,out,m,subDirPath); //write restart
    if ( ! in.outMeshFileName.empty() )
      m->writeNative(in.outMeshFileName.c_str());
    if ( in.writeGeomBCFiles ) {
      if(!PCU_Comm_Self()) printf("write additional geomBC file for visualization\n");
      // store the value of the function pointer
      FILE* (*fn)(Output& out, const char* path) = out.openfile_write;
      // set function pointer for file writing
      out.openfile_write = chef::openfile_write;
      ph::writeGeomBC(out, path, in.timeStepNumber); //write geombc for viz only
      // reset the function pointer to the original value
      out.openfile_write = fn;
    }
    ph::writeGeomBC(out, subDirPath); //write geombc
    if(!PCU_Comm_Self())
      ph::writeAuxiliaryFiles(path, in.timeStepNumber);
    m->verify();
#ifdef HAVE_SIMMETRIX
    gmi_model* g = m->getModel();
    ph::clearAttAssociation(g,in);
#endif
    if (in.adaptFlag)
      ph::goToParentDir();
    if(in.printIOtime) phastaio_printStats();
  }
  void preprocess(apf::Mesh2* m, Input& in, Output& out) {
    gmi_model* g = m->getModel();
    PCU_ALWAYS_ASSERT(g);
    BCs bcs;
    ph::readBCs(g, in.attributeFileName.c_str(), in.axisymmetry, bcs);
    if (!in.solutionMigration)
      ph::attachZeroSolution(in, m);
    if (in.buildMapping)
      ph::buildMapping(m);
    preprocess(m,in,out,bcs);
  }
}

namespace chef {
  void bake(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& in, ph::Output& out) {
    PCU_ALWAYS_ASSERT(PCU_Comm_Peers() % in.splitFactor == 0);
    apf::Migration* plan = 0;
    ph::BCs bcs;
    loadCommon(in, bcs, g);
    const int worldRank = PCU_Comm_Self();
    switchToMasters(in.splitFactor);
    if ((worldRank % in.splitFactor) == 0)
      originalMain(m, in, g, plan);
    switchToAll();
    if (in.simmetrixMesh == 0)
      m = repeatMdsMesh(m, g, plan, in.splitFactor);
    ph::checkBalance(m,in);
    ph::preprocess(m,in,out,bcs);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m) {
    ph::Input in;
    in.load("adapt.inp");
    in.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,in,out);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl) {
    ctrl.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,ctrl,out);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, GRStream* grs) {
    ctrl.openfile_read = openfile_read;
    ph::Output out;
    out.openfile_write = openstream_write;
    out.grs = grs;
    bake(g,m,ctrl,out);
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* rs) {
    ctrl.openfile_read = openstream_read;
    ctrl.rs = rs;
    ph::Output out;
    out.openfile_write = openfile_write;
    bake(g,m,ctrl,out);
    return;
  }
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* rs, GRStream* grs) {
    ctrl.openfile_read = openstream_read;
    ctrl.rs = rs;
    ph::Output out;
    out.openfile_write = openstream_write;
    out.grs = grs;
    bake(g,m,ctrl,out);
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
    ph::balanceAndReorder(m,ctrl,PCU_Comm_Peers());
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

