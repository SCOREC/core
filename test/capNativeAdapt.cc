#include <exception>
#include <string>
#include <vector>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMesh2.h>
#include <lionPrint.h>
#include <gmi_cap.h>
#include <pcu_util.h>

#include <CreateMG_AppProcessor.h>
#include <CreateMG_Framework_Mesh.h>
#include <CreateMG_Function.h>
#include <CreateMG_SizingMetricTool.h>

namespace {

/** \brief Print nested exceptions. */
void print_exception(const std::exception& e, int level = 0);

} // namespace

int main(int argc, char* argv[]) {
  lion_set_verbosity(1);
  pcu::Init(&argc, &argv);
  try {
    pcu::PCU PCU;
    if (argc < 4) {
      std::cerr << "USAGE: " << argv[0] << " <in.cre> <out.cre> <sizing.dat>"
        " <in.vmap>" << std::endl;
      throw 1;
    }

    const char *creFile = argv[1], *outFile = argv[2], *sizingFile = argv[3],
      *vmapFile = argv[4];

    gmi_cap_start();
    gmi_register_cap();
    try {
      gmi_model* model = gmi_cap_load(creFile);
      apf::Mesh2* mesh = apf::createCapMesh(model, &PCU);
      try {
        auto mdbi = apf::exportCapNative(mesh);
        auto ctx = mdbi->get_context();
        auto proc = CreateMG::get_context_processor(ctx);
        auto fn = CreateMG::get_function(ctx, "AdaptMesh");
        CreateMG::M_MModel mmodel;
        MG_API_CALL(mdbi, get_current_model(mmodel));
        CreateMG::set_input(fn, "MeshModel", mmodel);
        // CreateMG::set_input(fn, "Analysis", mmodel);
        auto sTool = CreateMG::get_sizing_metric_tool(
          ctx, "CreateSmoothingBase"
        );
        PCU_ALWAYS_ASSERT(sTool);
        sTool->set_context(ctx);
        std::vector<CreateMG::Metric6> sizing6;
        apf::loadCapSizingFileMetrics(mesh, sizingFile, vmapFile, sizing6);
        sTool->set_metric(mmodel, "sizing6", sizing6);
        CreateMG::set_input(fn, "TensorData", "sizing6");
        // DiscreteCurvature, Complexity, MaxComplexity
        auto t0 = pcu::Time();
        if (proc->execute(fn) != CreateMG::STATUS_OK) {
          std::cerr << "ERROR: failed to adapt mesh" << std::endl;
          throw 1;
        }
        auto t1 = pcu::Time();
        std::cout << "Capstone AdaptMesh ran in " << t1 - t0 << std::endl;
        CreateMG::get_output(fn, "Output", mmodel);
        std::string oname;
        MG_API_CALL(mdbi, get_model_name(mmodel, oname));
        apf::Mesh2* omesh = apf::createCapMesh(model, oname.c_str(), &PCU);
        apf::disownCapModel(omesh);
        apf::writeVtkFiles((std::string(outFile) + ".vtk").c_str(), omesh);
        gmi_cap_write(model, outFile);
        apf::destroyMesh(omesh);
      } catch (...) {
        apf::destroyMesh(mesh);
        std::rethrow_exception(std::current_exception());
      }
      apf::destroyMesh(mesh);
    } catch (const std::exception& e) {
      if (PCU.Self() == 0) {
        std::cerr << "ERROR: ";
        print_exception(e);
      }
      gmi_cap_stop();
      throw 1;
    } catch (...) {
      gmi_cap_stop();
      throw 1;
    }
    gmi_cap_stop();
  } catch (...) {
    pcu::Finalize();
    return 1;
  }
  pcu::Finalize();
  return 0;
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
