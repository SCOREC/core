#include <cstring>
#include <cstdlib>
#include <random>

// Output
#include <lionPrint.h>

// Parallelism
#include <PCU.h>
#include <pcu_util.h>

// Mesh interfaces
#include <apf.h>
#include <apfCAP.h>

// Geometry interfaces
#include <gmi.h>
#include <gmi_cap.h>

// Mesh adapt
#include <ma.h>

// Geometry
#include <gmi.h>
#include <gmi_cap.h>

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;

namespace {

void myExit(int exit_code = EXIT_SUCCESS) {
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
  exit(exit_code);
}

void printUsage(char *argv0) {
  printf("USAGE: %s <create_file.cre> <surface id>\n", argv0);
}

void analyticClosestPoint(double const from[3], double to[3], double to_norm[3]) {
  double x0 = from[0];
  // x axis axisymmetry, y0 here is distance from x axis
  double y0 = std::sqrt(from[1]*from[1] + from[2]*from[2]);

  // For Mathematica CForm output
  #define Power(base, exp) std::pow(base,exp)
  #define Sqrt(arg) std::sqrt(arg)

  // Parabola
  double pA =  -0.0516559000000000;
  double pZ = -3.0997300000000000;
  double p_y = (-1 - 2*pA*pZ + 2*pA*x0)/(Power(6,0.3333333333333333)*Power(9*Power(pA,4)*y0 + Sqrt(3)*Sqrt(Power(pA,6)*(2*Power(1 + 2*pA*(pZ - x0),3) + 27*Power(pA,2)*Power(y0,2))),0.3333333333333333)) + Power(9*Power(pA,4)*y0 + Sqrt(3)*Sqrt(Power(pA,6)*(2*Power(1 + 2*pA*(pZ - x0),3) + 27*Power(pA,2)*Power(y0,2))),0.3333333333333333)/ (Power(6,0.6666666666666666)*Power(pA,2));
  double p_x = pA*p_y*p_y + pZ;
  double p_d2 = std::pow(p_x-x0,2) + std::pow(p_y-y0,2);

  // Cone
  double A = 1000;
  double B = 1786.99;
  double C = -19.1427;
  double l_d2 = Power(C + A*x0 + B*y0,2)/(Power(A,2) + Power(B,2));
  double l_x = (-(A*C) + B*(B*x0 - A*y0))/(Power(A,2) + Power(B,2));
  double l_y = (-(B*C) + A*(-(B*x0) + A*y0))/(Power(A,2) + Power(B,2));

  // Convert to x y z
  // Normals
  double cls_y;
  double dx;
  double dy;
  if (p_d2 < l_d2) {
    to[0] = p_x;
    cls_y = p_y; 
    dx = -1.0;
    dy = 2*pA*p_y;
  } else {
    to[0] = l_x;
    cls_y = l_y;
    dx = A;
    dy = B;
  }

  double ratio = cls_y/y0;
  to[1] = from[1]*ratio;
  to[2] = from[2]*ratio;

  // Build normal vector in x y z
  to_norm[0] = dx;
  double norm_ratio = dy/y0;
  to_norm[1] = from[1]*norm_ratio;
  to_norm[2] = from[2]*norm_ratio;
  double norm_norm = std::sqrt(to_norm[0]*to_norm[0] + to_norm[1]*to_norm[1] + to_norm[2]*to_norm[2]);
  to_norm[0] /= norm_norm;
  to_norm[1] /= norm_norm;
  to_norm[2] /= norm_norm;
}

} // namespace

int main(int argc, char** argv) {
  // Initialize parallelism.
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  // Initialize logging.
  lion_set_stdout(stdout);
  lion_set_stderr(stderr);

  // Check arguments or print usage.
  if (argc < 3) {
    if (PCU_Comm_Self() == 0) {
      printUsage(argv[0]);
    }
    myExit(EXIT_FAILURE);
  }

  // Parse arguments.
  const char* createFileName = argv[1];
  int surf_id = atoi(argv[2]);

  // Initialize GMI.
  gmi_cap_start();
  gmi_register_cap();
  // create an instance of the Capstone Module activating SMLIB/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  const std::string gdbName("Geometry Database : SMLIB");
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("capClosestBenchmark", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();

  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  // Load Capstone mesh.
  v_string filenames;
  filenames.push_back(createFileName);
  M_GModel gmodel = cs.load_files(filenames);

  // Get surface for closest point tests
  M_GTopo surf = g->get_topo_by_id(Geometry::FACE, surf_id);
  if (surf.is_invalid()) {
    std::cerr << "ERROR: Shock surface id " << surf_id << " is invalid." << std::endl;
    myExit(EXIT_FAILURE);
  }
  gmi_model* gGMI = gmi_import_cap(g);
  gmi_ent* surfGMI = toGmiEntity(surf);

  // Generate random points
  double bmin[3];
  double bmax[3];
  gmi_bbox(gGMI, surfGMI, bmin, bmax);

  int npts = 500;
  double random[npts][3];
  std::uniform_real_distribution<double> unifx(bmin[0],bmax[0]);
  std::uniform_real_distribution<double> unify(bmin[1],bmax[1]);
  std::uniform_real_distribution<double> unifz(bmin[2],bmax[2]);
  std::default_random_engine re;
  for(int i=0; i<npts; i++) {
    random[i][0] = unifx(re);
    random[i][1] = unify(re);
    random[i][2] = unifz(re);
  }
  lion_oprint(0, "%.6f %.6f %.6f\n", random[2][0], random[2][1], random[2][2]);
  lion_oprint(0, "%.6f %.6f %.6f\n", random[362][0], random[362][1], random[362][2]);

  // gmi_closest point and gmi_normal
  double t0_a = MPI_Wtime();
  for(int i=0; i<npts; i++) {
    double clsArr[3], clsParArr[2];
    gmi_closest_point(gGMI, surfGMI, random[i], clsArr, clsParArr);
    double nrmArr[3];
    gmi_normal(gGMI, surfGMI, clsParArr, nrmArr);
  }
  double t1_a = MPI_Wtime();
  lion_oprint(0, "gmi_closest point and gmi_normal: Total runtime %.12f s, average runtime %.12f s \n", t1_a-t0_a, (t1_a-t0_a)/npts);

  // analytic methods
  t0_a = MPI_Wtime();
  for(int i=0; i<npts; i++) {
    double clsArr[3], nrmArr[3];
    analyticClosestPoint(random[i], clsArr, nrmArr);
  }
  t1_a = MPI_Wtime();
  lion_oprint(0, "gmi_closest point and gmi_normal: Total runtime %.12f s, average runtime %.12f s \n", t1_a-t0_a, (t1_a-t0_a)/npts);

  // Clean up.

  // Exit calls.
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}