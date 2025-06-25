#include <iostream>
#include <string>
#include <vector>

#include <PCU.h>
#include <gmi.h>
#include <gmi_cap.h>

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "USAGE: " << argv[0] << " <input.cre> [meshname...]"
      << std::endl;
    return 1;
  }
  pcu::Init(&argc, &argv);
  gmi_cap_start();
  gmi_register_cap();
  std::string model_content;
  std::vector<std::string> mesh_names;
  for (int i = 2; i < argc; ++i) mesh_names.push_back(argv[i]);
  gmi_model* model = gmi_cap_load_selective(argv[1], mesh_names);
  gmi_destroy(model);
  gmi_cap_stop();
  pcu::Finalize();
  return 0;
}
