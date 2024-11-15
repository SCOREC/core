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
    ARGS_GETTER(ratio, double)
    ARGS_GETTER(before, const std::string&)
    ARGS_GETTER(after, const std::string&)
    ARGS_GETTER(input, const std::string&)
    ARGS_GETTER(output, const std::string&)
    #undef ARGS_GETTER

    /** @brief Check for argument parse errors. */
    operator bool() const { return !error_flag_; }

    void print_usage(std::ostream& str) const;
    void print_help(std::ostream& str) const;

  private:
    bool analytic_{false}, error_flag_{false}, help_{false}, mds_adapt_{false}, smoothing_{false};
    int maxiter_{-1}, uniform_{0}, verbosity_{0};
    double aniso_size_{0.0}, thickness_{0.0}, ratio_{4.0};
    std::string argv0, before_, after_, input_, output_, ref_shock_geom_;
    std::list<int> shock_surf_ids_{};
  }; // class Args

  class EmbeddedShockFunction : public ma::AnisotropicFunction {
  public:
    EmbeddedShockFunction(ma::Mesh* m, gmi_model* g, std::list<gmi_ent*> surfs, double nsize,
      double AR, double h0, double thickness) : mesh(m),
      norm_size(nsize), init_size(h0), ref(g), shock_surfaces(surfs) {
      //shock_surface = reinterpret_cast<apf::ModelEntity*>(toGmiEntity(shock));
      thickness_tol = thickness * thickness / 4;
      tan_size = norm_size * AR;
    }
    EmbeddedShockFunction(ma::Mesh* m, std::list<gmi_ent*> surfs, double nsize,
      double AR, double h0, double thickness) : 
      EmbeddedShockFunction(m, m->getModel(), surfs, nsize, AR, h0, thickness) {};
    void getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale);
    double getMaxEdgeLengthAcrossShock();
  protected:
    double getZoneIsoSize(ma::Entity* vtx);
    double getZoneIsoSize(ma::Entity* vtx, apf::Vector3 closestPt);
    ma::Mesh* mesh;
    gmi_model* ref;
    double thickness_tol, norm_size, init_size, tan_size;
    std::list<gmi_ent*> shock_surfaces;
    int nInShockBand;
  }; // class EmbeddedSizeFunction

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
  std::cout << "Isotropic adapt only" << std::endl;

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
  std::cout << "INFO: Selected surfaces: ";
  for(int surf_id : args.shock_surf_ids()) {
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
  //std::cout << "INFO: Selected surface: " << shock_surface << std::endl;

  std::cout << ++stage << ". Get average edge length." << std::endl;
  double h0 = ma::getAverageEdgeLength(adaptMesh);
  std::cout << "INFO: Average edge length: " << h0 << std::endl;
  double h0_max = ma::getMaximumEdgeLength(adaptMesh, nullptr);
  std::cout << "INFO: Maximum edge length: " << h0_max << std::endl;

  std::cout << "INFO: Normal size: " << args.aniso_size() << std::endl;
  std::cout << "INFO: Tangent size: " << args.aniso_size() * args.ratio() << std::endl;

  std::cout << ++stage << ". Make sizefield." << std::endl;
  ma::AnisotropicFunction* sf = nullptr;
  if(using_ref_geom) {
    sf = new EmbeddedShockFunction(adaptMesh,
      geom_ref, shock_surfaces, args.aniso_size(), args.ratio(), h0, args.thickness());
    //double h0_max_shock = reinterpret_cast<AnisotropicFunctionOnReference*>(sf)->getMaxEdgeLengthAcrossShock();
    //std::cout << "INFO: Maximum edge length crossing shock: " << h0_max_shock << std::endl;
  } else {
    sf = new EmbeddedShockFunction(adaptMesh, shock_surfaces, args.aniso_size(), args.ratio(), h0, args.thickness());
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
    const char* required = "nst";
    while ((c = getopt(argc, argv, ":A:aB:hi:mn:o:r:G:s:t:uvS")) != -1) {
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
      case 'o':
        output_ = optarg;
        break;
      case 'r':
        ratio_ = std::atof(optarg);
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
    if (thickness_ < 0) {
      std::cerr << "ERROR: Thickness must be positive." << std::endl;
      error_flag_ = true;
    }
    if (aniso_size_ < 0) {
      std::cerr << "ERROR: Aniso size must be positive." << std::endl;
      error_flag_ = true;
    }
    if (optind < argc) {
      input_ = argv[optind];
    } else {
      std::cerr << "ERROR: INPUT.cre is required." << std::endl;
      error_flag_ = true;
    }
  } 

  void Args::print_usage(std::ostream& str) const {
    str << "USAGE: " << argv0 << " [-ahmuvS] [-i MAXITER] [-B before.vtk] "
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
    "-S              Do smoothing \n"
    "-n NORM_SIZE    Set anisotropic normal direction size. (required)\n"
    "-o OUTPUT.cre   Write final mesh to OUTPUT.cre.\n"
    "-r RATIO        Set desired anisotropic aspect ratio (tan/norm)."
      " DEFAULT: 4\n"
    "-G REF_GEOM.cre Set a different .cre file with shock geometry.\n"
    "-s ID           Set face IDs of the embedded shock surface. (required)\n"
    "                Repeat flag multiple times to specify more tags.\n"
    "-t THICKNESS    Set thickness (required).\n"
    "-u              Perform uniform adaptation. Specifying multiple times\n"
    "                runs that many rounds of adaptation.\n"
    "-v              Increase verbosity. Level 1 enables lionPrint. Level 2\n"
    "                enables verbose adaptation. Level 3 writes intermediate\n"
    "                VTK files.\n";
  }

  double EmbeddedShockFunction::getZoneIsoSize(ma::Entity* vtx) {
    apf::Vector3 pos;
    mesh->getPoint(vtx, 0, pos);
    apf::Vector3 sphere_cent(-0.250,0,0);
    apf::Vector3 dist = pos - sphere_cent;
    bool in_sphere = std::abs(dist * dist) < 0.4226 * 0.4226;
    return in_sphere ? 0.10565625 : 0.2113125;
  }

  double EmbeddedShockFunction::getZoneIsoSize(ma::Entity* vtx, apf::Vector3 closestPt) {
    apf::Vector3 pos;
    mesh->getPoint(vtx, 0, pos);
    apf::Vector3 vecToPos = pos - closestPt;
    // slight negative tolerance for outer outlet edge
    return vecToPos.x() > -1e-3 ? 4 * getZoneIsoSize(vtx) : getZoneIsoSize(vtx);
  }

  void EmbeddedShockFunction::getValue(ma::Entity* vtx, ma::Matrix& frame, ma::Vector& scale) {
    apf::Vector3 pos;
    mesh->getPoint(vtx, 0, pos);
    double posArr[3];
    pos.toArray(posArr);

    double shockDistSquare = DBL_MAX;
    double clsArr[3], clsParArr[2];
    apf::Vector3 clsVec;
    gmi_ent* closestSurf;
    PCU_ALWAYS_ASSERT(gmi_can_get_closest_point(ref));
    for(gmi_ent* surf : shock_surfaces) {
      //gmi_closest_point (struct gmi_model *m, struct gmi_ent *e, double const from[3], double to[3], double to_p[2])
      double outClsArr[3], outClsParArr[2];
      gmi_closest_point(ref, surf, posArr, outClsArr, outClsParArr);
      double curShockDistSquare = (posArr[0]-outClsArr[0])*(posArr[0]-outClsArr[0])+
        (posArr[1]-outClsArr[1])*(posArr[1]-outClsArr[1])+
        (posArr[2]-outClsArr[2])*(posArr[2]-outClsArr[2]);
      if (curShockDistSquare < shockDistSquare) {
        shockDistSquare = curShockDistSquare;
        clsArr[0] = outClsArr[0]; clsArr[1] = outClsArr[1]; clsArr[2] = outClsArr[2];
        clsParArr[0] = outClsParArr[0]; clsParArr[1] = outClsParArr[1]; clsParArr[2] = outClsParArr[2];
        closestSurf = surf;
        clsVec.fromArray(clsArr);
      }
    }

    if (shockDistSquare < thickness_tol) {
      apf::Vector3 norm;
      nInShockBand++;
      /*
      if (shockDistSquare > 1e-9) {
        norm = (pos-clsVec).normalize();
      } else {
        double posPerturbedArr[3] = {posArr[0]+0.001, posArr[1], posArr[2]};
        double clsPerturbedArr[3], clsPerturbedParArr[3];
        gmi_closest_point(ref, closestSurf, posPerturbedArr, clsPerturbedArr, clsPerturbedParArr);
        norm = apf::Vector3(posPerturbedArr[0]-clsPerturbedArr[0], posPerturbedArr[1]-clsPerturbedArr[1],
          posPerturbedArr[2]-clsPerturbedArr[2]).normalize();
      }
      */

      PCU_ALWAYS_ASSERT(gmi_has_normal(ref));
      double nrmArr[3];
      gmi_normal(ref, closestSurf, clsParArr, nrmArr);
      norm.fromArray(nrmArr);

      // Negate largest component to get tangent.
      apf::Vector3 trial(norm[2], norm[1], norm[0]);
      int largest = trial[0] > trial[1] && trial[0] > trial[2] ? 0 : (trial[1] > trial[0] && trial[1] > trial[2] ? 1 : 2);
      trial[largest] *= -1;
      apf::Vector3 tan1 = apf::cross(norm, trial).normalize();
      apf::Vector3 tan2 = apf::cross(norm, tan1);

      frame[0][0] = nrmArr[0]; frame[0][1] = tan1[0]; frame[0][2] = tan2[0];
      frame[1][0] = nrmArr[1]; frame[1][1] = tan1[1]; frame[1][2] = tan2[1];
      frame[2][0] = nrmArr[2]; frame[2][1] = tan1[2]; frame[2][2] = tan2[2];

      double zoneIsoSize = getZoneIsoSize(vtx);
      scale[0] = norm_size;
      scale[1] = zoneIsoSize;
      scale[2] = zoneIsoSize;

      if(lion_get_verbosity() >= 1 && nInShockBand % 500 == 0){
        std::cout << posArr[0] << " " << posArr[1] << " " << posArr[2] << " ";
        std::cout << clsArr[0] << " " << clsArr[1] << " " << clsArr[2] << " ";
        std::cout << nrmArr[0] << " " << nrmArr[1] << " " << nrmArr[2] << " ";
        std::cout << tan1[0] << " " << tan1[1] << " " << tan1[2] << " ";
        std::cout << tan2[0] << " " << tan2[1] << " " << tan2[2] << std::endl;
      }

    } else {
      frame[0][0] = 1; frame[0][1] = 0; frame[0][2] = 0;
      frame[1][0] = 0; frame[1][1] = 1; frame[1][2] = 0;
      frame[2][0] = 0; frame[2][1] = 0; frame[2][2] = 1;
      //scale[0] = scale[1] = scale[2] = init_size;
      scale[0] = scale[1] = scale[2] = getZoneIsoSize(vtx, clsVec);
    }
  }

  double EmbeddedShockFunction::getMaxEdgeLengthAcrossShock() {
    double max_length = -1.0;
    apf::MeshIterator* it = mesh->begin(1);
    apf::MeshEntity* e;
    while ((e = mesh->iterate(it))) {
      for (gmi_ent* surf : shock_surfaces) {
        apf::MeshEntity* adj_pts[2];
        mesh->getDownward(e, 0, adj_pts);
        apf::Vector3 pointA, pointB, closest;
        mesh->getPoint(adj_pts[0], 0, pointA);
        mesh->getPoint(adj_pts[1], 0, pointB);

        double posArr[3], clsArr[3], clsParArr[2];
        pointA.toArray(posArr);
        gmi_closest_point(ref, surf, posArr, clsArr, clsParArr);
        closest.fromArray(clsArr);

        apf::Vector3 vecA = pointA - closest;
        apf::Vector3 vecB = pointB - closest;
        if (vecA * vecB < 0) {
          max_length = std::max(max_length, (vecB-vecA).getLength());
          break;
        }
      }
    }
    mesh->end(it);
    return max_length;
  }
} // namespace