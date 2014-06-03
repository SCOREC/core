#include <phInput.h>
#include <phBC.h>
#include <phRestart.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>

int main(int argc, char** argv)
{
  ph::Input in("adapt.inp");
  apf::Mesh2* m = apf::loadMdsMesh(
      in.modelFileName.c_str(), in.meshFileName.c_str());
  ph::BCs bcs;
  ph::readBCs(in.attributeFileName.c_str(), bcs);
  if (in.solutionMigration)
    ph::readAndAttachSolution(in, m);
  m->destroyNative();
  apf::destroyMesh(m);
}
