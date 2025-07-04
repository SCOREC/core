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

namespace {
void print_exception(const std::exception& e, int level = 0);
}

int main(int argc, char* argv[]) {
  int retval = 0;
  pcu::Init(&argc, &argv);
  {
  pcu::PCU PCU;
  try {
    if (argc != 5) {
      if (PCU.Self() == 0)
        std::cerr << "USAGE: <model.dmg> <mesh.smb> <inParts> <outMesh.smb>\n"
          "\nwhere inParts < PCU.Peers()"
          << std::endl;
      throw std::runtime_error("invalid arguments");
    }
    lion_set_verbosity(1);
    gmi_register_mesh();
    // load model and mesh
    int inParts = std::stoi(argv[3]);
    if (inParts >= PCU.Peers()) {
      throw std::runtime_error("inParts >= PCU.Peers()");
    }
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
    m->verify();
    m->writeNative(argv[4]);
    // destroy mds
    m->destroyNative();
    apf::destroyMesh(m);
  } catch (const std::exception& e) {
    if (PCU.Self() == 0) {
      std::cerr << "ERROR: ";
      print_exception(e);
    }
    retval = 1;
  } catch (...) {
    if (PCU.Self() == 0)
      std::cerr << "Unknown exception occurred." << std::endl;
    retval = 1;
  }
  } // PCU object scope
  pcu::Finalize();
  return retval;
}

namespace {

void print_exception(const std::exception& e, int level) {
  std::cerr << std::string(level * 2, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception& nestedE) {
    print_exception(nestedE, level + 1);
  } catch (...) {}
}

} // namespace
