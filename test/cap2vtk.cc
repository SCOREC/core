#include <iostream>
#include <stdexcept>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <pcu_util.h>

int main(int argc, char** argv)
{
  int retval = 0;
  pcu::Init(&argc, &argv);
  {
  pcu::PCU PCUObj;
  lion_set_verbosity(1);
  try {
    if (PCUObj.Peers() > 1)
      throw std::runtime_error("Capstone meshes are (currently) only serial");

    if (argc != 3 && argc != 4) {
      if (PCUObj.Self() == 0) {
        std::cerr << "USAGE: " << argv[0]
          << " <input.cre> [meshname] <output.vtk>" << std::endl;
      }
      throw std::invalid_argument("invalid arguments");
    }

    const char* creFileName = argv[1];
    const char* meshName = argc == 4 ? argv[2] : nullptr;
    const char* folderName = argc == 4 ? argv[3] : argv[2];

    gmi_cap_start();
    gmi_register_cap();
    try {
      std::string model_content;
      std::vector<std::string> mesh_names, mesh_contents;
      gmi_cap_probe(creFileName, model_content, mesh_names, mesh_contents);
      if (meshName) {
        if (std::count(mesh_names.begin(), mesh_names.end(), meshName) == 0) {
          throw std::invalid_argument(
            "mesh named " + std::string(meshName) + " not found"
          );
        }
      } else if (mesh_names.size() > 1) {
        std::cout << "More than one mesh in CRE file."
          " Run again with meshname argument after CRE filename." << std::endl;
        for (size_t i = 0; i < mesh_names.size(); ++i) {
          std::cout << "Mesh " << i << " `" << mesh_names[i] << "`: "
            << mesh_contents[i] << std::endl;
        }
        throw std::runtime_error(
          "meshname not provided for file with multiple meshes"
        );
      } else meshName = mesh_names.front().c_str();
      PCU_DEBUG_ASSERT(meshName);
      gmi_model* model = gmi_cap_load_selective(creFileName, {meshName});
      apf::Mesh2* mesh = apf::createCapMesh(model, meshName, &PCUObj);
      apf::writeVtkFiles(folderName, mesh);
      apf::destroyMesh(mesh);
    } catch (...) {
      gmi_cap_stop();
      std::rethrow_exception(std::current_exception());
    }
    gmi_cap_stop();
  } catch (const std::exception& e) {
    if (PCUObj.Self() == 0) std::cerr << "ERROR: " << e.what() << std::endl;
    retval = 1;
  } catch (...) {
    if (PCUObj.Self() == 0) std::cerr << "UNKNOWN ERROR" << std::endl;
    retval = 1;
  }
  } // PCU Object scope
  pcu::Finalize();
  return retval;
}
