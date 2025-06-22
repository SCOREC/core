#include <iostream>
#include <exception>
#include <memory>
#include <string>
#include <apf.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <gmi_mesh.h>
#include <apfMETIS.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char* argv[]) {
  pcu::Init(&argc, &argv);
  try {
    pcu::PCU PCU;
    if (argc != 5) {
      if (PCU.Self() == 0)
        std::cerr << "USAGE: <model.dmg> <mesh.smb> <inParts> <out.smd>"
          << std::endl;
      throw std::runtime_error("invalid arguments");
    }
    lion_set_verbosity(1);
    gmi_register_mesh();
    // load model and mesh
    int inParts = std::stoi(argv[3]);
    int group = PCU.Self() / inParts;
    auto loadPCU = PCU.Split(group, 0);
    gmi_model* model = gmi_load(argv[1]);
    apf::Mesh2* m = nullptr;
    if (group == 0) {
      m = apf::loadMdsMesh(model, argv[2], loadPCU.get());
      m->switchPCU(&PCU);
    }
    m = apf::expandMdsMesh(m, model, inParts, &PCU);
    try {
      std::unique_ptr<apf::Balancer> balancer(apf::makeMETISbalancer(m));
      double imbalance = 1.1;
      balancer->balance(nullptr, imbalance);
    } catch (...) {
      std::throw_with_nested(std::runtime_error("balancing failed"));
    }
    m->writeNative(argv[4]);
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
