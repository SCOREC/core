#include <exception>
#include <stdexcept>

#include <PCU.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfConvert.h>
#include <apfMDS.h>
#include <apfMETIS.h>
#include <apfMesh2.h>
#include <apfPartition.h>
#include <apfZoltan.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <lionPrint.h>
#include <ma.h>
#include <parma.h>
#include <pcu_util.h>

#include "capVolSizeFields.h"

namespace {

/** \brief Print nested exceptions. */
void print_exception(
  const pcu::PCU& PCU, const std::exception& e, int level = 0
);

struct Args {
  Args() {}
  void parse(int argc, char* argv[]);
  static void print_usage(const char* argv0);
  std::string before, after, in, out;
  bool analytic{false}, volume{false}, mds{false};
  int sf{0};
  enum Partitioner { Default, Zoltan, METIS, Parma } splitter, balancer;
};

void Args::print_usage(const char *argv0) {
  std::cout << "USAGE: " << argv0 << " [-agm] [-B BEFORE.VTK] [-A AFTER.VTK] "
    "[-s SPLITTER] [-b BALANCER] <size-field> <IN.CRE> <OUT.CRE>" << std::endl;
  std::cout << R"help(
Flags:
-B BEFORE.VTK  Write the file BEFORE.VTK before adaptation with the
               adapt_frames and adapt_scales fields.
-A AFTER.VTK   Write the file AFTER.VTK after adaptation.
-a             Evaluate size-field analytically during adaptation. The default
               is to evaluate once, write the frames/scales, then transfer and
               interpolate during adaptation.
-g             Force mesh volume generation.
-m             Convert mesh to MDS during adaptation (required for parallel
               adaptation).
-s SPLITTER    Force the selected splitter to be used. Possible values are:
               Zoltan, METIS, Parma. The default is to use whatever is
               available, in the order listed previously.
-b BALANCER    Force the selected balancer to be used during adaptation.
               Possible values and default priorities are the same as for
               splitter selection. Using this option forces load-balancing
               after each step, whereas the default only balances when the
               estimated imbalance exceeds a threshold (10%%).

SIZE-FIELDS:
1, for uniform anisotropic size-field
2, for wing-shock size-field
3, for cube-shock size-field
4, for cylinder boundary-layer size-field
)help";
}

ma::Mesh* loadAdaptMesh(
  pcu::PCU* PCU, gmi_model* model, bool volume, bool mds
);

void parallelAdapt(
  pcu::PCU* PCU, gmi_model* model, apf::Mesh2* mesh, const Args& args
);

} // namespace

int main(int argc, char** argv) {
  lion_set_verbosity(1); // Initialize logging.
  int retval = 0;
  // Initialize parallelism.
  pcu::Init(&argc, &argv);
  {
  pcu::PCU PCUObj;
  try {
    // Check arguments or print usage.
    Args args;
    try {
      args.parse(argc, argv);
    } catch (const std::exception& e) {
      if (PCUObj.Self() == 0) args.print_usage(argv[0]);
      std::throw_with_nested(
        std::runtime_error("invalid command line arguments")
      );
    }

    if (PCUObj.Peers() > 1 && !args.mds)
      throw std::runtime_error("parallel run without -m flag");

    // Initialize GMI.
    gmi_cap_start();
    gmi_register_cap();

    try {
      gmi_model* capGeomModel = nullptr;
      if (PCUObj.Self() == 0) capGeomModel = gmi_cap_load(args.in.c_str());
      else capGeomModel = gmi_cap_load_selective(args.in.c_str(), {});
      auto soloPCU = PCUObj.Split(PCUObj.Self(), 0);
      ma::Mesh* mesh = nullptr;
      if (PCUObj.Self() == 0) {
        mesh = loadAdaptMesh(
          soloPCU.get(), capGeomModel, args.volume, args.mds
        );
        // APF default routine will typically fail to verify surface meshes.
        if (args.volume) mesh->verify();
      }
      parallelAdapt(&PCUObj, capGeomModel, mesh, args);
      if (PCUObj.Self() == 0) {
        if (args.volume) mesh->verify();
        if (args.mds) {
          apf::Mesh2* mdsMesh = mesh;
          mesh = apf::makeEmptyCapMesh(
            capGeomModel, "MeshAdapt", soloPCU.get()
          );
          apf::disownCapModel(mesh);
          apf::convert(mdsMesh, mesh);
          apf::destroyMesh(mdsMesh);
        }
        gmi_cap_write(capGeomModel, args.out.c_str());
        apf::destroyMesh(mesh);
      }
      gmi_destroy(capGeomModel);
      gmi_cap_stop();
    } catch(...) {
      gmi_cap_stop();
      std::rethrow_exception(std::current_exception());
    }
  } catch (const std::exception& e) {
    if (PCUObj.Self() == 0) std::cerr << "ERROR: ";
    print_exception(PCUObj, e);
    retval = 1;
  } catch(...) {
    if (PCUObj.Self() == 0)
      std::cerr << "ERROR: unspecified error" << std::endl;
    retval = 1;
  }
  } // PCUObj scope
  pcu::Finalize();
  return retval;
}

namespace {

void print_exception(const pcu::PCU& PCU, const std::exception& e, int level) {
  if (PCU.Self() == 0)
    std::cerr << std::string(level * 2, ' ') << e.what() << '\n';
  try {
    std::rethrow_if_nested(e);
  } catch (const std::exception& nestedE) {
    print_exception(PCU, nestedE, level + 1);
  } catch (...) {}
}

void Args::parse(int argc, char* argv[]) {
  constexpr int positional_args = 3;
  if (argc < 1 + positional_args)
    throw std::runtime_error("missing positional argument(s)");
  auto strarg = [argc, argv](int& i, int& j) -> std::string {
    if (argv[i][j + 1]) {
      int oldj = j;
      while (argv[i][j + 1] != '\0') ++j;
      return &argv[i][oldj + 1];
    } else if (i + 1 < argc - positional_args) {
      ++i;
      for (j = 0; argv[i][j + 1] != '\0'; ++j);
      return argv[i];
    } else throw std::runtime_error("missing argument for -" + argv[i][j]);
  };
  auto strlower = [](std::string str) {
    for (auto& c : str) c = std::tolower(c);
    return str;
  };
  auto str2part = [strlower](std::string str) -> Partitioner {
    str = strlower(str);
    if (str == "metis") return METIS;
    if (str == "parma") return Parma;
    if (str == "zoltan") return Zoltan;
    throw std::invalid_argument("invalid partitioner name: " + str);
  };
  int i;
  for (i = 1; i < argc - positional_args; ++i) {
    if (*argv[i] == '-') {
      for (int j = 1; argv[i][j] != '\0'; ++j) {
        switch(argv[i][j]) {
        case 'A': after = strarg(i, j); break;
        case 'B': before = strarg(i, j); break;
        case 'a': analytic = true; break;
        case 'g': volume = true; break;
        case 'm': mds = true; break;
        case 's': splitter = str2part(strarg(i, j)); break;
        case 'b': balancer = str2part(strarg(i, j)); break;
        default: throw std::runtime_error("unrecognized flag: -" + argv[i][j]);
        }
      }
    } else break;
  }
  sf = std::atoi(argv[i]);
  in = argv[i + 1];
  out = argv[i + 2];
}

ma::Mesh* loadAdaptMesh(
  pcu::PCU* PCU, gmi_model* model, bool volume, bool mds
) {
  // Load Capstone mesh (with optional volume generation).
  ma::Mesh* capMesh = nullptr;
  if (volume) {
    int dim = 3;
    capMesh = apf::generateCapMesh(model, dim, PCU);
    // FIXME: maybe create copy of mesh model to work on (preserve original).
  } else capMesh = apf::createCapMesh(model, PCU);
  apf::disownCapModel(capMesh);
  // Optionally convert to MDS
  ma::Mesh* adaptMesh = nullptr;
  if (mds) {
    adaptMesh = apf::createMdsMesh(model, capMesh, true);
    apf::disownMdsModel(adaptMesh); // Model is managed in main.
    apf::destroyMesh(capMesh);
  } else adaptMesh = capMesh;
  return adaptMesh;
}

apf::Splitter* makeSplitter(Args::Partitioner ptnr, apf::Mesh2* mesh) {
  switch (ptnr) {
  case Args::Zoltan:
    return apf::makeZoltanSplitter(mesh, apf::GRAPH, apf::PARTITION);
  case Args::METIS: return apf::makeMETISsplitter(mesh);
  case Args::Parma: return Parma_MakeRibSplitter(mesh);
  default:
    #if defined(PUMI_HAS_ZOLTAN)
    return apf::makeZoltanSplitter(
      mesh, apf::GRAPH, apf::PARTITION
    );
    #elif defined(PUMI_HAS_METIS)
    return apf::makeMETISsplitter(mesh);
    #else
    return Parma_MakeRibSplitter(mesh);
    #endif
  }
}

ma::AnisotropicFunction* makeUDF(int mode, ma::Mesh* mesh) {
  switch (mode) {
    case 1: return new UniformAniso(mesh);
    case 2: return new WingShock(mesh, 50);
    case 3: return new Shock2(mesh);
    case 4: return new CylBoundaryLayer(mesh);
    default: throw std::runtime_error("invalid size-field");
  }
}

void migrateHome(apf::Mesh2* mesh) {
  auto t0 = pcu::Time();
  apf::Migration* plan = new apf::Migration(mesh);
  apf::MeshIterator* it = mesh->begin(mesh->getDimension());
  for (apf::MeshEntity* e; (e = mesh->iterate(it));) plan->send(e, 0);
  mesh->end(it);
  mesh->migrate(plan); // destroys plan
  auto map0 = apf::Multiply(0);
  apf::remapPartition(mesh, map0);
  auto t1 = pcu::Time();
  if (mesh->getPCU()->Self() == 0)
    std::cout << "INFO: Migrated home in " << t1 - t0 << " seconds"
      << std::endl;
}

void parallelAdapt(
  pcu::PCU* PCU, gmi_model* model, apf::Mesh2* mesh, const Args& args
) {
  pcu::PCU* oldPCU = nullptr;
  if (mesh) oldPCU = mesh->getPCU();
  if (PCU->Peers() > 1) {
    apf::Migration* plan = nullptr;
    if (PCU->Self() == 0) {
      apf::Splitter* splitter = makeSplitter(args.splitter, mesh);
      apf::MeshTag* weights = Parma_WeighByMemory(mesh);
      plan = splitter->split(weights, 1.10, PCU->Peers());
      apf::removeTagFromDimension(
        mesh, weights, mesh->getDimension()
      );
      mesh->destroyTag(weights);
      delete splitter;
      mesh->switchPCU(PCU);
    }
    mesh = apf::repeatMdsMesh(mesh, model, plan, PCU->Peers(), PCU);
    apf::disownMdsModel(mesh);
  }
  // Choose appropriate size-field.
  std::unique_ptr<ma::AnisotropicFunction> sf(makeUDF(args.sf, mesh));

  // Make pumi fields for the frames and scales for anisotropic size-fields.
  apf::Field* frameField = nullptr;
  apf::Field* scaleField = nullptr;
  if (!args.analytic || !args.before.empty()) {
    frameField = apf::createFieldOn(mesh, "adapt_frames", apf::MATRIX);
    scaleField = apf::createFieldOn(mesh, "adapt_scales", apf::VECTOR);
    apf::MeshIterator* it = mesh->begin(0);
    for (apf::MeshEntity* v; (v = mesh->iterate(it));) {
      ma::Vector s;
      ma::Matrix f;
      sf->getValue(v, f, s);
      apf::setVector(scaleField, v, 0, s);
      apf::setMatrix(frameField, v, 0, f);
    }
    mesh->end(it);

    if (!args.before.empty()) apf::writeVtkFiles(args.before.c_str(), mesh);
    if (args.analytic) { // Cleanup if fields were only for before.vtk
      apf::destroyField(frameField);
      apf::destroyField(scaleField);
      frameField = scaleField = nullptr;
    }
  }

  ma::Input *in = nullptr;
  if (args.analytic) in = ma::makeAdvanced(ma::configure(mesh, sf.get()));
  else in = ma::makeAdvanced(ma::configure(mesh, scaleField, frameField));
  switch (args.balancer) {
  case Args::Zoltan:
    in->shouldRunPreZoltan = in->shouldRunMidZoltan = true;
    in->shouldRunPostZoltan = true;
    break;
  case Args::METIS:
    in->shouldRunPreMetis = in->shouldRunMidMetis = true;
    in->shouldRunPostMetis = true;
    break;
  case Args::Parma:
    in->shouldRunPreParma = in->shouldRunMidParma = true;
    in->shouldRunPostParma = true;
    break;
  default:;
  }

  ma::adapt(in);
  if (!args.after.empty()) apf::writeVtkFiles(args.after.c_str(), mesh);
  if (PCU->Peers() > 1) migrateHome(mesh);
  if (PCU->Self() != 0) apf::destroyMesh(mesh);
  else mesh->switchPCU(oldPCU);
}

} // namespace
