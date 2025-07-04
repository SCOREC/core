#include <iostream>
#include <string>
#include <vector>

#include <apf.h>
#include <apfMesh.h>
#include <apfCAP.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <PCU.h>

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "USAGE: " << argv[0] << " <input.cre>" << std::endl;
    return 1;
  }
  pcu::Init(&argc, &argv);
  gmi_cap_start();
  gmi_register_cap();
  std::string model_content;
  std::vector<std::string> mesh_names, mesh_contents;
  gmi_cap_probe(argv[1], model_content, mesh_names, mesh_contents);
  std::cout << "model: " << model_content << std::endl;
  for (std::size_t i = 0; i < mesh_names.size(); ++i) {
    std::cout << "mesh (" << i << ") `" << mesh_names[i] << "`: "
      << mesh_contents[i] << std::endl;
  }
  gmi_cap_stop();
  pcu::Finalize();
  return 0;
}
