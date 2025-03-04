#include <unistd.h>
#include <cstdlib>
#include <PCU.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <apf.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <apfConvert.h>
#include <ma.h>
#include <float.h>

#include <CapstoneModule.h>
#include "embedded_aniso_function.h"

using namespace CreateMG;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

namespace {
  class Args {
  public:
    Args(int argc, char* argv[]);

    #define ARGS_GETTER(name, type) type name(void) const noexcept { return name##_; }
    ARGS_GETTER(analytic, bool)
    ARGS_GETTER(help, bool)
    ARGS_GETTER(mds_adapt, bool)
    ARGS_GETTER(smoothing, bool)
    ARGS_GETTER(ref_shock_geom, const std::string&)
    ARGS_GETTER(shock_surf_ids, std::list<int>)
    ARGS_GETTER(uniform, int)
    ARGS_GETTER(verbosity, int)
    ARGS_GETTER(maxiter, int)
    ARGS_GETTER(aniso_size, double) 
    ARGS_GETTER(thickness, double)
    ARGS_GETTER(global_size, double)
    ARGS_GETTER(before, const std::string&)
    ARGS_GETTER(after, const std::string&)
    ARGS_GETTER(input, const std::string&)
    ARGS_GETTER(output, const std::string&)
    ARGS_GETTER(isotropic, bool)
    ARGS_GETTER(refine_only, bool)
    #undef ARGS_GETTER

    /** @brief Check for argument parse errors. */
    operator bool() const { return !error_flag_; }

    void print_usage(std::ostream& str) const;
    void print_help(std::ostream& str) const;

  private:
    bool analytic_{false}, isotropic_{false}, error_flag_{false}, help_{false}, mds_adapt_{false}, smoothing_{false}, refine_only_{false};
    int maxiter_{-1}, uniform_{0}, verbosity_{0};
    double aniso_size_{-1}, thickness_{-1}, global_size_{-1};
    std::string argv0, before_, after_, input_, output_, ref_shock_geom_;
    std::list<int> shock_surf_ids_{};
  }; // class Args
} // namespace

int main(int argc, char* argv[]) {
  // Initalize parallelism.
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  Args args(argc, argv);
  if (args.help()) {
    args.print_help(std::cout);
    PCU_Comm_Free();
    MPI_Finalize();
    return !args ? 1 : 0;
  } else if (!args) {
    args.print_usage(std::cerr);
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }

  if (args.verbosity() > 0) {
    lion_set_verbosity(1);
  }

  int stage = 0;
  std::cout << ++stage << ". Setup Capstone." << std::endl;

  CapstoneModule cs("cap_aniso", "Geometry Database : SMLIB",
    "Mesh Database : Create", "Attribution Database : Create");
  GeometryDatabaseInterface* gdi = cs.get_geometry();
  MeshDatabaseInterface* mdi = cs.get_mesh();

  std::cout << ++stage << ". Load CRE." << std::endl;
  M_GModel gmodel = cs.load_files(v_string(1, args.input()));
  gmi_model* gmodelGMI = gmi_import_cap(gdi);
  std::vector<M_MModel> mmodels;
  MG_API_CALL(mdi, get_associated_mesh_models(gmodel, mmodels));
  MG_API_CALL(mdi, set_current_model(mmodels[0])); // Use first mesh.

  std::cout << ++stage << ". Setup adjacencies." << std::endl;
  MG_API_CALL(mdi, set_adjacency_state(REGION2FACE|REGION2EDGE|REGION2VERTEX|
    FACE2EDGE|FACE2VERTEX));
  MG_API_CALL(mdi, set_reverse_states());
  MG_API_CALL(mdi, compute_adjacency());

  std::cout << ++stage << ". Initalize gmi." << std::endl;
  gmi_register_cap();
  gmi_cap_start();

  std::cout << ++stage << ". Make apf::Mesh." << std::endl;
  apf::Mesh2* apfCapMesh = apf::createMesh(mdi, gdi);
  apf::printStats(apfCapMesh);

  ma::Mesh* adaptMesh = apfCapMesh;
  if (args.mds_adapt()) {
    std::cout << ++stage << ". Convert mesh to MDS." << std::endl;
    adaptMesh = apf::createMdsMesh(apfCapMesh->getModel(), apfCapMesh);
    apf::disownMdsModel(adaptMesh);
    apf::printStats(adaptMesh);
  }

  if (args.uniform() > 0) {
    std::cout << ++stage << ". Run uniform refinement (" << args.uniform() <<
      " rounds)." << std::endl;
    ma::runUniformRefinement(adaptMesh, args.uniform());
  }

  CapstoneModule* cs_ref = nullptr;
  GeometryDatabaseInterface* gdi_ref = nullptr;
  M_GModel gmodel_ref;
  gmi_model* geom_ref = nullptr;
  bool using_ref_geom = args.ref_shock_geom().length() > 0;
  if(using_ref_geom) {
    std::cout << ++stage << ". Load reference geometry." << std::endl;
    std::cout << "INFO: Reference geometry file: " << args.ref_shock_geom() << std::endl;
    
    cs_ref = new CapstoneModule("cap_aniso_ref", "Geometry Database : SMLIB",
      "Mesh Database : Create", "Attribution Database : Create");
    gdi_ref = cs_ref->get_geometry();
    gmodel_ref = cs_ref->load_files(v_string(1, args.ref_shock_geom()));
    geom_ref = gmi_import_cap(gdi_ref);
  }

  // TODO: change this for multiple ids
  std::cout << ++stage << ". Identify shock surfaces by ids." << std::endl;

  // M_GTopo shock_surface = gdi->get_topo_by_id(Geometry::FACE, args.shock_surf_id());
  std::list<gmi_ent*> shock_surfaces;
  if (args.shock_surf_ids().size() == 0) {
    std::cout << "INFO: No surface ids specified - using analytic shock geometry definitions" << std::endl;
  } else {
    std::cout << "INFO: Selected surfaces: ";
    for (int surf_id : args.shock_surf_ids()) {
      std::cout << surf_id << " ";
      M_GTopo surf = (using_ref_geom ? gdi_ref : gdi)->get_topo_by_id(Geometry::FACE, surf_id);
      if (surf.is_invalid()) {
        std::cerr << "ERROR: Shock surface id " << surf_id << " is invalid." << std::endl;
        PCU_Comm_Free();
        MPI_Finalize();
        return 1;
      }
      shock_surfaces.push_back(toGmiEntity(surf));
    }
    std::cout << std::endl;
  }
  //std::cout << "INFO: Selected surface: " << shock_surface << std::endl;

  std::cout << ++stage << ". Get average edge length." << std::endl;
  double h0 = ma::getAverageEdgeLength(adaptMesh);
  std::cout << "INFO: Average edge length: " << h0 << std::endl;
  double h0_max = ma::getMaximumEdgeLength(adaptMesh, nullptr);
  std::cout << "INFO: Maximum edge length: " << h0_max << std::endl;

  std::cout << "INFO: Normal size: " << args.aniso_size() << std::endl;

  std::cout << ++stage << ". Make sizefield." << std::endl;
  ma::AnisotropicFunction* sf = nullptr;
  if(using_ref_geom) {
    sf = new EmbeddedShockFunction(adaptMesh,
      geom_ref, shock_surfaces, args.isotropic(), args.global_size(), args.aniso_size(), args.thickness());
    //double h0_max_shock = reinterpret_cast<AnisotropicFunctionOnReference*>(sf)->getMaxEdgeLengthAcrossShock();
    //std::cout << "INFO: Maximum edge length crossing shock: " << h0_max_shock << std::endl;
  } else {
    sf = new EmbeddedShockFunction(adaptMesh, shock_surfaces, args.isotropic(), args.global_size(), args.aniso_size(), args.thickness());
  }
  apf::Field *frameField = nullptr, *scaleField = nullptr;
  if (!args.before().empty() || !args.analytic()) {
    frameField = apf::createFieldOn(adaptMesh, "adapt_frames", apf::MATRIX);
    scaleField = apf::createFieldOn(adaptMesh, "adapt_scales", apf::VECTOR);
    ma::Iterator* it = adaptMesh->begin(0);
    for (ma::Entity *v = adaptMesh->iterate(it); v;
      v = adaptMesh->iterate(it)) {
      ma::Vector scale;
      ma::Matrix frame;
      sf->getValue(v, frame, scale);
      apf::setVector(scaleField, v, 0, scale);
      apf::setMatrix(frameField, v, 0, frame);
    }
    adaptMesh->end(it);

    if (args.smoothing()) {
      std::cout << ++stage << ". Smooth size field." << std::endl;
      if(apf::has_smoothCAPAnisoSizes()) {
        apf::smoothCAPAnisoSizes(adaptMesh, "cap_aniso", scaleField, frameField);
      } else {
        std::cout << "apf::smoothCAPAnisoSizes is not supported, skipping" << std::endl;
      }
    }

    if (!args.before().empty()) {
      std::cout << ++stage << ". Write before VTK." << std::endl;
      apf::writeVtkFiles(args.before().c_str(), adaptMesh);
    }
  }

  std::cout << ++stage << ". Setup adaptation." << std::endl;
  ma::Input* in;
  if (args.analytic()) {
    in = ma::makeAdvanced(ma::configure(adaptMesh, sf));
  } else {
    in = ma::makeAdvanced(ma::configure(adaptMesh, scaleField, frameField));
  }
  if (args.maxiter() >= 0) in->maximumIterations = args.maxiter();

  if(args.refine_only()) {
    std::cout << "Doing refinement only (no coarsening, snapping, or shape fix)." << std::endl;
    in->shouldCoarsen = false;
    in->shouldSnap=false;
    in->shouldFixShape=false;
  }

  std::cout << ++stage << ". Run adapt." << std::endl;
  if (args.verbosity() > 1) {
    ma::adaptVerbose(in, args.verbosity() > 2);
  } else {
    ma::adapt(in);
  }

  if (!args.after().empty()) {
    std::cout << ++stage << ". Write after VTK." << std::endl;
    apf::writeVtkFiles(args.after().c_str(), adaptMesh);
  }

  if (!args.output().empty()) {
    if (args.mds_adapt()) {
      std::cout << ++stage << ". Convert mesh to Capstone." << std::endl;
      M_MModel mmodel;
      MG_API_CALL(mdi, create_associated_model(mmodel, gmodel, "MeshAdapt"));
      MG_API_CALL(mdi, set_adjacency_state(REGION2FACE|REGION2EDGE|
        REGION2VERTEX|FACE2EDGE|FACE2VERTEX));
      MG_API_CALL(mdi, set_reverse_states());
      MG_API_CALL(mdi, compute_adjacency());
      apf::convert(adaptMesh, apfCapMesh);
#ifndef NDEBUG
      std::cout << "=== DEBUG ===\n";
      std::map<std::pair<int, int>, int> dimTagToCount_mds, dimTagToCount_cap;
      for (int d = 0; d <= adaptMesh->getDimension(); ++d) {
        apf::MeshIterator* it = adaptMesh->begin(d);
        for (apf::MeshEntity* e = adaptMesh->iterate(it); e; e = adaptMesh->iterate(it)) {
          apf::ModelEntity* me = adaptMesh->toModel(e);
          std::pair<int, int> key(adaptMesh->getModelType(me), adaptMesh->getModelTag(me));
          if (dimTagToCount_mds.count(key) == 1) {
            dimTagToCount_mds[key] += 1;
          } else {
            dimTagToCount_mds[key] = 1;
          }
        }
        adaptMesh->end(it);
        it = apfCapMesh->begin(d);
        for (apf::MeshEntity* e = apfCapMesh->iterate(it); e; e = apfCapMesh->iterate(it)) {
          apf::ModelEntity* me = apfCapMesh->toModel(e);
          std::pair<int, int> key(apfCapMesh->getModelType(me), apfCapMesh->getModelTag(me));
          if (dimTagToCount_mds.count(key) == 1) {
            dimTagToCount_cap[key] += 1;
          } else {
            dimTagToCount_cap[key] = 1;
          }
        }
        apfCapMesh->end(it);
      }
      std::ofstream mdsTagCount_f("dimtagct_mds.csv"),
        capTagCount_f("dimtagct_cap.csv");
      mdsTagCount_f << "dim,tag,count\n";
      for (const auto& pairs : dimTagToCount_mds) {
        mdsTagCount_f << pairs.first.first   << ',' << pairs.first.second << ','
          << pairs.second << '\n';
      }
      capTagCount_f << "dim,tag,count\n";
      for (const auto& pairs : dimTagToCount_cap) {
        capTagCount_f << pairs.first.first << ',' << pairs.first.second << ','
          << pairs.second << '\n';
      }

      // Parametric values.
      apf::Field* mds_debug_params = apf::createFieldOn(adaptMesh,
        "debug_params", apf::VECTOR);
      apf::MeshIterator* it = adaptMesh->begin(0);
      for (apf::MeshEntity* v = adaptMesh->iterate(it); v; v = adaptMesh->iterate(it)) {
        apf::Vector3 param;
        adaptMesh->getParam(v, param);
        apf::setVector(mds_debug_params, v, 0, param);
      }
      adaptMesh->end(it);
      apf::writeVtkFiles("mds_params.vtk", adaptMesh);
      apf::Field* cap_debug_params = apf::createFieldOn(apfCapMesh,
        "debug_params", apf::VECTOR);
      it = apfCapMesh->begin(0);
      for (apf::MeshEntity* v = apfCapMesh->iterate(it); v; v = apfCapMesh->iterate(it)) {
        apf::Vector3 param;
        adaptMesh->getParam(v, param);
        apf::setVector(cap_debug_params, v, 0, param);
      }
      adaptMesh->end(it);
      apf::writeVtkFiles("cap_params.vtk", adaptMesh);
      std::cout << "=== END DEBUG ===\n";
#endif
      apf::destroyMesh(adaptMesh);
    }
    std::cout << ++stage << ". Write final CRE." << std::endl;
    AppContext* ctx = cs.get_context();
    Writer* creWriter = get_writer(ctx, "Create Native Writer");
    if (!creWriter) {
      std::cerr << "FATAL: Could not find the CRE writer." << std::endl;
      gmi_cap_stop();
      PCU_Comm_Free();
      MPI_Finalize();
      return 1;
    }
    IdMapping idmap;
    std::vector<M_MModel> mmodels;
    MG_API_CALL(mdi, get_associated_mesh_models(gmodel, mmodels));
    creWriter->write(ctx, gmodel, mmodels, args.output().c_str(), idmap);
  }

  std::cout << ++stage << ". Cleanup." << std::endl;
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

namespace {
  Args::Args(int argc, char* argv[]) {
    argv0 = argv[0];
    int c;
    int given[256] = {0};
    const char* required = "";
    while ((c = getopt(argc, argv, ":A:aB:hi:mn:o:r:G:s:t:uvSIRg:")) != -1) {
      ++given[c];
      switch (c) {
      case 'A':
        after_ = optarg;
        break;
      case 'a':
        analytic_ = true;
        break;
      case 'B':
        before_ = optarg;
        break;
      case 'h':
        help_ = true;
        break;
      case 'i':
        maxiter_ = std::atoi(optarg);
        break;
      case 'm':
        mds_adapt_ = true;
        break;
      case 'n':
        aniso_size_ = std::atof(optarg);
        break;
      case 'g':
        global_size_ = std::atof(optarg);
        break;
      case 'o':
        output_ = optarg;
        break;
      case 'r':
        std::cout << "INFO: Ratio argument is deprecated" << std::endl;
        //ratio_ = std::atof(optarg);
        break;
      case 's':
        shock_surf_ids_.push_back(std::atoi(optarg));
        std::cout << "INFO: Processing arg: " << optarg << " last entry " << shock_surf_ids_.back() << std::endl;
        break;
      case 'G':
        ref_shock_geom_ = optarg;
        break;
      case 't':
        thickness_ = std::atof(optarg);
        break;
      case 'u':
        uniform_++;
        break;
      case 'v':
        ++verbosity_;
        break;
      case 'I':
        isotropic_ = true;
        break;
      case 'R':
        refine_only_ = true;
        break;
      case 'S':
        smoothing_ = true;
        break;
      case ':':
        std::cerr << "ERROR: Option -" << char(optopt) << " requires an "
          "argument." << std::endl;
        error_flag_ = true;
        break;
      case '?':
        std::cerr << "ERROR: Unrecognized option: -" << char(optopt) << std::endl;
        error_flag_ = true;
        break;
      }
    }
    for (const char* r = required; *r != '\0'; ++r) {
      if (!given[int(*r)]) {
        std::cerr << "ERROR: Flag -" << *r << " is required." << std::endl;
        error_flag_ = true;
      }
    }
    if (optind < argc) {
      input_ = argv[optind];
    } else {
      std::cerr << "ERROR: INPUT.cre is required." << std::endl;
      error_flag_ = true;
    }
  } 

  void Args::print_usage(std::ostream& str) const {
    str << "USAGE: " << argv0 << " [-ahmuvSIR] [-i MAXITER] [-B before.vtk] "
      "[-A after.vtk] [-o OUTPUT.cre] [-G REF_GEOM.cre] -s ID -n NORM_SIZE -t THICKNESS INPUT.cre"
      << std::endl;
  }

  void Args::print_help(std::ostream& str) const {
    print_usage(str);
    str << "cap_aniso adapts the mesh with an embedded shock surface.\n";
    str << "OPTIONS:\n"
    "-A after.vtk    Write adapted mesh to after.vtk.\n"
    "-a              Use analytic function instead of interpolating size field.\n"
    "-B before.vtk   Write initial mesh with size field to before.vtk.\n"
    "-h              Display this help menu.\n"
    "-i MAXITER      Set maximum adapt iterations. Negative values mean to\n"
    "                use MeshAdapt interpreted value based on size-field.\n"
    "                DEFAULT: -1.\n"
    "-m              Convert to mesh to MDS before adaptation (and back to \n"
    "                Capstone mesh before writing if given -o). May boost \n"
    "                performance, but buggy.\n"
    "-S              Do smoothing. \n"
    "-n NORM_SIZE    Override anisotropic normal direction size.\n"
    "-o OUTPUT.cre   Write final mesh to OUTPUT.cre.\n"
    "-G REF_GEOM.cre Set a different .cre file with shock geometry.\n"
    "-s ID           Set face IDs of the embedded shock surface. (required)\n"
    "                Repeat flag multiple times to specify more tags.\n"
    "-t THICKNESS    Override thickness.\n"
    "-g h_global     Override h_global."
    "-u              Perform uniform adaptation. Specifying multiple times\n"
    "                runs that many rounds of adaptation.\n"
    "-I              Run completely isotropic adaptation.\n"
    "-R              Run refinement only (disable coarsening, snapping and shape fix)"
    "-v              Increase verbosity. Level 1 enables lionPrint. Level 2\n"
    "                enables verbose adaptation. Level 3 writes intermediate\n"
    "                VTK files.\n";
  }
} // namespace