#include <ph.h>
#include <chef.h>
#include <phstream.h>
#ifdef HAVE_SIMMETRIX
#include <phAttrib.h>
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
#include <apfSIM.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_sim.h>
#include <SimPartitionedMesh.h>
#include <PCU.h>
#include <pcu_io.h>
#include <string>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>

#define SIZET(a) static_cast<size_t>(a)

namespace {
static int mesh_has_ext(const char* filename, const char* ext)
{
  const char* c = strrchr(filename, '.');
  if (!c) {
    if (PCU_Comm_Self()==0)
      fprintf(stderr, "mesh file name with no extension");
    assert(c);
  }
  ++c; /* exclude the dot itself */
  if (!strcmp(c, ext))
    return 1;
  return 0;
}


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

static apf::Mesh2* loadMesh(gmi_model*& g, const char* meshfile) {
  apf::Mesh2* mesh;
  /* if it is a simmetrix mesh */
  if (mesh_has_ext(meshfile, "sms")) {
    pProgress progress = Progress_new();
    Progress_setDefaultCallback(progress);

    pGModel simModel = gmi_export_sim(g);
    pParMesh sim_mesh = PM_load(meshfile, sthreadNone, simModel, progress);
    apf::Mesh* simApfMesh = apf::createMesh(sim_mesh);

    mesh = apf::createMdsMesh(g, simApfMesh);

    apf::destroyMesh(simApfMesh);
    M_release(sim_mesh);
    Progress_delete(progress);
  }
  /* if it is a SCOREC mesh */
  else {
    mesh = apf::loadMdsMesh(g, meshfile);
  }
  return mesh;
}

void originalMain(apf::Mesh2*& m, ph::Input& in,
    gmi_model* g, apf::Migration*& plan)
{
  if(!m)
    m = loadMesh(g, in.meshFileName.c_str());
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
    ph::migrateInterfaceItr(m, bcs);
    if(in.timeStepNumber > 0)
      ph::checkReorder(m,in,PCU_Comm_Peers());
    if (in.adaptFlag)
      ph::goToStepDir(in.timeStepNumber);
    std::string path = ph::setupOutputDir();
    ph::setupOutputSubdir(path);
    ph::enterFilteredMatching(m, in, bcs);
    ph::generateOutput(in, bcs, m, out);
    ph::exitFilteredMatching(m);
    if ( ! in.outMeshFileName.empty() )
      m->writeNative(in.outMeshFileName.c_str());
    // a path is not needed for inmem
    ph::detachAndWriteSolution(in,out,m,path); //write restart
    if (in.adaptFlag && (in.timeStepNumber % in.writeVizFiles == 0) ) {
      // store the value of the function pointer
      FILE* (*fn)(Output& out, const char* path) = out.openfile_write;
      // set function pointer for file writing
      out.openfile_write = chef::openfile_write;
      ph::writeGeomBC(out, path, in.timeStepNumber); //write geombc for viz only
      // reset the function pointer to the original value
      out.openfile_write = fn;
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
    const int numMasters = PCU_Comm_Peers();
    if ((worldRank % in.splitFactor) == 0)
      originalMain(m, in, g, plan);
    switchToAll();
    m = repeatMdsMesh(m, g, plan, in.splitFactor);
    ph::balanceAndReorder(m,in,numMasters);
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
    int numOfComp) {
    return ph::extractField(m,packedFieldname,
             requestFieldname,firstComp,numOfComp);
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

