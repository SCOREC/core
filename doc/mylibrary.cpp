#include <gmi_mesh.h>
#include <apfBox.h>
#include <apfMesh2.h>
#include <apf.h>
void makeMesh() {
  gmi_register_mesh();
  apf::Mesh2* m = apf::makeMdsBox(1,1,1,1,1,1,0);
  apf::destroyMesh(m);
}
