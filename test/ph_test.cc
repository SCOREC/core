#include <ph.h>
#include <phInput.h>
#include <phBC.h>
#include <phRestart.h>
#include <phAdapt.h>
#include <phOutput.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <gmi_mesh.h>
#include <PCU.h>

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  gmi_register_mesh();
  ph::Input in("adapt.inp");
  apf::Mesh2* m = apf::loadMdsMesh(
      in.modelFileName.c_str(), in.meshFileName.c_str());
  ph::BCs bcs;
  ph::readBCs(in.attributeFileName.c_str(), bcs);
  if (in.solutionMigration)
    ph::readAndAttachSolution(in, m);
  if (in.adaptFlag) {
    ph::adapt(in, m);
    ph::goToStepDir(in.timeStepNumber);
  }
  if (in.tetrahedronize)
    ph::tetrahedronize(in, m);
  std::string path = ph::setupOutputDir();
  ph::setupOutputSubdir(path);
  if (in.phastaIO) {
    if (in.solutionMigration)
      ph::detachAndWriteSolution(in, m, path);
    ph::Output o;
    ph::generateOutput(in, bcs, m, o);
    ph::writeGeomBC(o, path);
  }
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
