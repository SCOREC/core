#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdexcept>
#include <unistd.h>
#include <getopt.h>
#include <random>

#include <lionPrint.h>
#include <ma.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <apfZoltan.h>
#include <parma.h>
#include <apfSIM.h>
#include <MeshSim.h>
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <PCU.h>
#include <PCU_C.h>
#include <pcu_util.h>

// === includes for safe_mkdir ===
#include <reel.h>
#include <sys/types.h> /*required for mode_t for mkdir on some systems*/
#include <sys/stat.h> /*using POSIX mkdir call for SMB "foo/" path*/
#include <errno.h> /* for checking the error from mkdir */
// ===============================

#include "SimParasolidKrnl.h"
#include "MeshSimAdapt.h"
#include "SimDiscrete.h"
#include "SimAdvMeshing.h"
#include "SimMeshTools.h"

using std::cout;
using std::endl;
using std::cerr;
using std::ofstream;

namespace {
  class Args {
  public:
    Args(int argc, char* argv[]);

    #define ARGS_GETTER(name, type) type name(void) const noexcept { return name##_; }
    ARGS_GETTER(help, bool)
    ARGS_GETTER(input_mesh, const std::string&)
    ARGS_GETTER(input_model, const std::string&)
    ARGS_GETTER(input_nmodel, const std::string&)
    ARGS_GETTER(prefix, const std::string&)
    ARGS_GETTER(mds_maxiter, int)
    ARGS_GETTER(verbosity, int)
    ARGS_GETTER(random_seed, double)
    ARGS_GETTER(h0, double)
    ARGS_GETTER(aspect_ratio, double)
    ARGS_GETTER(max_angle, double)
    #undef ARGS_GETTER

    /** @brief Check for argument parse errors. */
    operator bool() const { return !error_flag_; }

    void print_usage(std::ostream& str) const;
    void print_help(std::ostream& str) const;

  private:
    bool help_{false}, error_flag_{false};
    double h0_{0.0625}, aspect_ratio_{1.0}, max_angle_{0};
    int verbosity_{0}, mds_maxiter_{-1}, random_seed_{999999};
    std::string argv0, input_mesh_, input_model_, input_nmodel_, prefix_;
  }; // class Args
} // namespace

void safe_mkdir(const char* path)
{
  mode_t const mode = S_IRWXU|S_IRGRP|S_IXGRP|S_IROTH|S_IXOTH;
  int err;
  errno = 0;
  err = mkdir(path, mode);
  if (err != 0 && errno != EEXIST)
  {
    reel_fail("Err: could not create directory \"%s\"\n", path);
  }
}

void compute_vertex_data(apf::Mesh2 *mesh_mds, apf::Field *frameField, apf::Field *scaleField, std::default_random_engine &generator, std::uniform_real_distribution<double> &distribution, Args args) {
  ma::Iterator* it = mesh_mds->begin(0);
  for (ma::Entity *v = mesh_mds->iterate(it); v;
    // size field
    v = mesh_mds->iterate(it)) {
    ma::Vector scale(args.h0()/args.aspect_ratio(), args.h0(), args.h0());
    ma::Matrix frame;

    double angle = distribution(generator);
    apf::Vector3 norm(std::cos(angle), std::sin(angle), 0);

    // Negate largest component to get tangent.
    apf::Vector3 trial(norm[2], norm[1], norm[0]);
    int largest = trial[0] > trial[1] && trial[0] > trial[2] ? 0 : (trial[1] > trial[0] && trial[1] > trial[2] ? 1 : 2);
    trial[largest] *= -1;
    apf::Vector3 tan1 = apf::cross(norm, trial).normalize();
    apf::Vector3 tan2 = apf::cross(norm, tan1);

    frame[0][0] = norm[0]; frame[0][1] = tan1[0]; frame[0][2] = tan2[0];
    frame[1][0] = norm[1]; frame[1][1] = tan1[1]; frame[1][2] = tan2[1];
    frame[2][0] = norm[2]; frame[2][1] = tan1[2]; frame[2][2] = tan2[2];

    apf::setVector(scaleField, v, 0, scale);
    apf::setMatrix(frameField, v, 0, frame);
  }
  mesh_mds->end(it);
}

void compute_vtx_max_AR(apf::Mesh2 *mesh_mds, std::string filename) {
  std::ofstream file;
  file.open(filename);

  ma::Iterator* it = mesh_mds->begin(3);
  for (ma::Entity *r = mesh_mds->iterate(it); r;
    r = mesh_mds->iterate(it)) {

    ma::Entity* e[6];
    mesh_mds->getDownward(r, 1, e);
    double longest_edge = 0;
    for (int i = 0; i < 6; i++) {
      apf::MeshEntity* vs[2];
      mesh_mds->getDownward(e[i], 0, vs);
      apf::Vector3 p0, p1;
      mesh_mds->getPoint(vs[0], 0, p0);
      mesh_mds->getPoint(vs[1], 0, p1);
      longest_edge = std::max(longest_edge, (p1 - p0).getLength());
    }
    double shortest_altitude = apf::computeShortestHeightInTet(mesh_mds, r);
    double ar = longest_edge / shortest_altitude;

    file << ar << " " << apf::getLinearCentroid(mesh_mds, r) << "\n";
  }
  mesh_mds->end(it);

  file.close();
}

int main(int argc, char *argv[])
{
  Args args(argc, argv);
  if (args.help()) {
    args.print_help(std::cout);
    return !args ? 1 : 0;
  } else if (!args) {
    args.print_usage(std::cerr);
    return 1;
  }

  MPI_Init(&argc, &argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);

  lion_set_verbosity(args.verbosity());

  cout << endl;
  cout <<" Reading ... " << endl;
  cout <<"  Model from file : " << args.input_model() << endl;
  cout <<"  Mesh from file  : " << args.input_mesh() << endl;
  cout <<" Output Prefix : " << args.prefix() << endl;

  safe_mkdir(args.prefix().data());

  // Init MeshSim
  Sim_logOn("perturbed_size_field.log");
  MS_init();
  SimAdvMeshing_start();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  SimPartitionedMesh_start(&argc,&argv);

  gmi_sim_start();
  gmi_register_sim();
  pProgress progress = Progress_new();
  Progress_setDefaultCallback(progress);

  // Simmetrix
  //pNativeModel nmodel = 0;
  pGModel model = 0;
  pParMesh mesh = 0;

  // SCOREC
  gmi_model* mdl_ref;
  ma::Mesh* mesh_ref;

  // Load geometry (closely following convert.cc)
  // todo: maybe some modules not loaded?
  if(!args.input_nmodel().empty()) {
    mdl_ref = gmi_sim_load(args.input_nmodel().data(),args.input_model().data());
  } else {
    mdl_ref = gmi_load(args.input_model().data());
  }
  model = gmi_export_sim(mdl_ref);

  // Load mesh
  mesh = PM_load(args.input_mesh().data(), model, progress);
  mesh_ref = apf::createMesh(mesh, &PCUObj);
  
  cout<<endl;
  cout<<" Reading model and mesh done ..."<<endl;
  cout<<"Initial mesh statistics: Num. elements: "<<M_numRegions(PM_mesh(mesh,0))<<", num. vertices: "<<M_numRegions(PM_mesh(mesh,0))<<endl;
  cout<<endl;

  // Try SCOREC meshadapt
  cout << " trying meshadapt" << endl;

  apf::Mesh2* mesh_mds = apf::createMdsMesh(mdl_ref, mesh_ref);

  std::cout << " writing size field to before vtk file" << std::endl;

  std::random_device rd;
  std::default_random_engine generator(rd());
  generator.seed(args.random_seed());
  double max_angle_rad = args.max_angle() * M_PI / 180.0;
  std::uniform_real_distribution<double> distribution(-max_angle_rad, max_angle_rad);

  apf::Field *frameField = nullptr, *scaleField = nullptr, *arField = nullptr;
  frameField = apf::createFieldOn(mesh_mds, "adapt_frames", apf::MATRIX);
  scaleField = apf::createFieldOn(mesh_mds, "adapt_scales", apf::VECTOR);
  compute_vertex_data(mesh_mds, frameField, scaleField, generator, distribution, args);
  compute_vtx_max_AR(mesh_mds, args.prefix()+"/before_AR");

  std::string before_name = args.prefix()+"/before";
  apf::writeVtkFiles(before_name.data(), mesh_mds);
  ma::Input* in = ma::makeAdvanced(ma::configure(mesh_mds, scaleField, frameField));

  if (args.mds_maxiter() >= 0) {
    cout << " maximum mds adapt iterations set to " << args.mds_maxiter() << endl;
    in->maximumIterations = args.mds_maxiter();
  }

  double t0 = PCU_Time();
  if (args.verbosity() > 1) {
    ma::adaptVerbose(in, args.verbosity() > 2);
  } else {
    ma::adapt(in);
  }
  cout << " MeshAdapt adapt time (s): " << PCU_Time()-t0 << endl;

  compute_vtx_max_AR(mesh_mds, args.prefix()+"/after_AR");
  compute_vertex_data(mesh_mds, frameField, scaleField, generator, distribution, args);
  cout<<"Adapted mesh statistics (MDS): Num. elements: "<<mesh_mds->count(3)<<", num. vertices: "<<mesh_mds->count(0)<<endl;
  cout << " writing meshdapt result" << endl;
  std::string after_name = args.prefix()+"/after";
  apf::writeVtkFiles(after_name.data(), mesh_mds);

  cout << " destroying meshdapt objects" << endl;
  apf::destroyMesh(mesh_mds);

  // Release everything
  M_release(mesh); 
  //GM_release(model);

  Progress_delete(progress);
  gmi_sim_stop();
  SimPartitionedMesh_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();

  MS_exit();
  Sim_logOff();
  }
  MPI_Finalize();

  return 0;
}

namespace {
  Args::Args(int argc, char* argv[]) {
    argv0 = argv[0];
    int c;
    int given[256] = {0};
    const char* required = "";
    while ((c = getopt(argc, argv, "hvP:i:S:")) != -1) {
      ++given[c];
      switch (c) {
      case 'h':
        help_ = true;
        break;
      case 'v':
        ++verbosity_;
        break;
      case 'P':
        prefix_ = optarg;
        break;
      case 'i':
        mds_maxiter_ = std::atoi(optarg);
        break;
      case 'S':
        random_seed_ = std::atoi(optarg);
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
    if (optind+4 < argc) {
      input_model_ = argv[optind];
      input_mesh_ = argv[optind+1];
      h0_ = std::stod(argv[optind+2]);
      aspect_ratio_ = std::stod(argv[optind+3]);
      max_angle_ = std::stod(argv[optind+4]);
    } else {
      std::cerr << "ERROR: MODEL.smd and INPUT.sms is required." << std::endl;
      error_flag_ = true;
    }
  }

  void Args::print_usage(std::ostream& str) const {
    str << "USAGE: " << argv0 << " [-hv] [-i mds maxiter] [-P output_prefix] [-G MODEL_nat.x_t] "
      "MODEL.smd INITIAL.sms h0 aspect_ratio max_angle"
      << std::endl;
  }

  void Args::print_help(std::ostream& str) const {
    print_usage(str);
    str << "OPTIONS:\n"
    "-h                 Display this help menu.\n"
    "-v                 Increase verbosity. \n"
    "-i mds MAXITER     Maximum MDS adapt iterations.\n"
    "-P                 Prefix to vtu output for organization. \n"
    "-S                 Random seed. \n"
    "-G MODEL_nat.x_t   Set a _nat file, model is discrete when not set.\n";
  }
} // namespace