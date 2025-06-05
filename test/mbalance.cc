#include <iostream>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <gmi_mesh.h>
#include <apfMETIS.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  pcu::Init(&argc, &argv);
  try {
  pcu::PCU pcu_obj;
  if (argc != 4) {
    if (pcu_obj.Self() == 0) {
      std::cerr << "USAGE: <model.dmg> <mesh.smd> <out.smd>" << std::endl;
    }
    throw 1;
  }
  lion_set_verbosity(1);
  gmi_register_mesh();
  //load model and mesh
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2],&pcu_obj);
  apf::Balancer* balancer = apf::makeMETISbalancer(m);
  int imbalance = 1.1;
  balancer->balance(nullptr, imbalance);
  delete balancer;
  m->writeNative(argv[3]);
  // destroy mds
  m->destroyNative();
  apf::destroyMesh(m);
  } catch (...) {
    pcu::Finalize();
    return 1;
  }
  pcu::Finalize();
  return 0;
}
