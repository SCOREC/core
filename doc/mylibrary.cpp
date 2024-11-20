#include <gmi_mesh.h>
#include <apfBox.h>
#include <apfMesh2.h>
#include <apf.h>
void makeMesh(pcu::PCU *PCUObj) {
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(1,1,1,1,1,1,0,PCUObj);
  apf::destroyMesh(m);
}
