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
#include <parma.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <PCU.h>
#include <pcu_io.h>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>

#define SIZET(a) static_cast<size_t>(a)

namespace {

#ifdef HAVE_SIMMETRIX
static bool mesh_has_ext(const char* filename, const char* ext)
{
  const char* c = strrchr(filename, '.');
  if (c) {
    ++c; /* exclude the dot itself */
    return !strcmp(c, ext);
  } else {
    return false;
  }
}
#endif

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
  apf::Mesh2* mesh;
  const char* meshfile = in.meshFileName.c_str();
#ifdef HAVE_SIMMETRIX
  /* if it is a simmetrix mesh */
  if (mesh_has_ext(meshfile, "sms")) {
    if (in.simmetrixMesh == 0) {
      if (PCU_Comm_Self()==0)
        fprintf(stderr, "oops, turn on flag: simmetrixMesh\n");
      in.simmetrixMesh = 1;
      in.filterMatches = 0; //not support
    }
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    pGModel simModel = gmi_export_sim(g);
    pParMesh sim_mesh = PM_load(meshfile, sthreadNone, simModel, progress);
    mesh = apf::createMesh(sim_mesh);

    Progress_delete(progress);
  } else
#endif
  /* if it is a SCOREC mesh */
  {
    mesh = apf::loadMdsMesh(g, meshfile);
  }
  return mesh;
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
    return pcu_group_open(path, false);
  }

  static FILE* openfile_write(ph::Output&, const char* path) {
    return pcu_group_open(path, true);
  }

  static FILE* openstream_write(ph::Output& out, const char* path) {
    return openGRStreamWrite(out.grs, path);
  }

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
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
    if(PCU_Comm_Peers() > 1)
      ph::migrateInterfaceItr(m, bcs);
    if (in.simmetrixMesh == 0)
      ph::checkReorder(m,in,PCU_Comm_Peers());
    if (in.adaptFlag)
      ph::goToStepDir(in.timeStepNumber,in.ramdisk);
    std::string path = ph::setupOutputDir(in.ramdisk);
    ph::setupOutputSubdir(path,in.ramdisk);
    ph::enterFilteredMatching(m, in, bcs);
    ph::generateOutput(in, bcs, m, out);
    ph::exitFilteredMatching(m);
    if ( ! in.outMeshFileName.empty() )
      m->writeNative(in.outMeshFileName.c_str());
    // a path is not needed for inmem
    ph::detachAndWriteSolution(in,out,m,path); //write restart
    if ( in.adaptFlag && in.writeGeomBCFiles ) {
      // store the value of the function pointer
      FILE* (*fn)(Output& out, const char* path) = out.openfile_write;
      // set function pointer for file writing
      out.openfile_write = chef::openfile_write;
      ph::writeGeomBC(out, path, in.timeStepNumber); //write geombc for viz only
      // reset the function pointer to the original value
      out.openfile_write = fn;
      in.writeGeomBCFiles = 0;
    }
    ph::writeGeomBC(out, path); //write geombc
    ph::writeAuxiliaryFiles(path, in.timeStepNumber);
    m->verify();
#ifdef HAVE_SIMMETRIX
    gmi_model* g = m->getModel();
    ph::clearAttAssociation(g,in);
#endif
    if (in.adaptFlag)
      ph::goToParentDir();
  }
  void preprocess(apf::Mesh2* m, Input& in, Output& out) {
    gmi_model* g = m->getModel();
    assert(g);
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
    assert(PCU_Comm_Peers() % in.splitFactor == 0);
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

