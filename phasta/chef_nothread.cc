#include <ph.h>
#include <phInput.h>
#include <phBC.h>
#include <phRestart.h>
#include <phAdapt.h>
#include <phOutput.h>
#include <phPartition.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <PCU.h>

#define SIZET(a) static_cast<size_t>(a)

namespace {

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

void afterSplit(apf::Mesh2* m, ph::Input& in, ph::BCs& bcs,
    int numMasters)
{
  std::string path = ph::setupOutputDir();
  ph::setupOutputSubdir(path);
  /* check if the mesh changed at all */
  if ((PCU_Comm_Peers()!=numMasters) ||
      in.adaptFlag ||
      in.tetrahedronize) {
    if (in.parmaPtn && PCU_Comm_Peers() > 1)
      ph::balance(m);
    apf::reorderMdsMesh(m);
  }
  ph::Output o;
  ph::generateOutput(in, bcs, m, o);
  ph::detachAndWriteSolution(in, m, path);
  ph::writeGeomBC(o, path);
  ph::writeAuxiliaryFiles(path, in.timeStepNumber);
  if ( ! in.outMeshFileName.empty() )
    m->writeNative(in.outMeshFileName.c_str());
  m->verify();
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

void loadCommon(ph::Input& in, ph::BCs& bcs,
    gmi_model*& g)
{
  in.load("adapt.inp");
  ph::readBCs(in.attributeFileName.c_str(), bcs);
  g = gmi_load(in.modelFileName.c_str());
}

void originalMain(apf::Mesh2*& m, ph::Input& in,
    gmi_model* g, apf::Migration*& plan)
{
  m = apf::loadMdsMesh(g, in.meshFileName.c_str());
  m->verify();
  if (in.solutionMigration)
    ph::readAndAttachSolution(in, m);
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

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  gmi_model* g = 0;
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  ph::Input in;
  ph::BCs bcs;
  loadCommon(in, bcs, g);
  const int worldRank = PCU_Comm_Self();
  switchToMasters(in.splitFactor);
  const int numMasters = PCU_Comm_Peers();
  if ((worldRank % in.splitFactor) == 0)
    originalMain(m, in, g, plan);
  switchToAll();
  if (in.adaptFlag)
    ph::goToStepDir(in.timeStepNumber);
  m = repeatMdsMesh(m, g, plan, in.splitFactor);
  afterSplit(m, in, bcs, numMasters);
  freeMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

