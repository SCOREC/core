#include <cstring>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include "cap_analytic_closest_point.h"

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
  printf("USAGE: %s <create_file.cre> <surface id> <bounds expand multiplier, default 1>\n", argv0);
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
  double expand = 1;
  if (argc > 3) {
    expand = atof(argv[3]);
  }

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
  #define bounds(index, expandby) bmin[index]-expandby*(bmax[index]-bmin[index]),bmax[index]+expandby*(bmax[index]-bmin[index])
  std::uniform_real_distribution<double> unifx(bounds(0,expand));
  std::uniform_real_distribution<double> unify(bounds(1,expand));
  std::uniform_real_distribution<double> unifz(bounds(2,expand));
  std::default_random_engine re;
  for(int i=0; i<npts; i++) {
    random[i][0] = unifx(re)*expand;
    random[i][1] = unify(re)*expand;
    random[i][2] = unifz(re)*expand;
  }

  // specific test points
  random[0][0] = 0;
  random[0][1] = 0;
  random[0][2] = 0;
  //
  random[1][0] = 1;
  random[1][1] = 0;
  random[1][2] = 0;
  //
  random[2][0] = 0;
  random[2][1] = 1;
  random[2][2] = 0;
  //
  random[3][0] = 0;
  random[3][1] = 0;
  random[3][2] = 1;
  //
  random[4][0] = 1;
  random[4][1] = 0.001;
  random[4][2] = 0.001;
  //
  random[5][0] = 1;
  random[5][1] = 0.000001;
  random[5][2] = 0.000001;

  double closestPts[npts][3];
  double normalVecs[npts][3];

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
    //double clsArr[3], nrmArr[3];
    doubleConeClosestPointAnalytic(random[i], closestPts[i], normalVecs[i]);
  }
  t1_a = MPI_Wtime();
  lion_oprint(0, "analytic: Total runtime %.12f s, average runtime %.12f s \n", t1_a-t0_a, (t1_a-t0_a)/npts);

  // write results
  std::ofstream outfile("closesttable.txt");
  for(int i=0; i<npts; i++) {
    outfile << random[i][0] << " " << random[i][1] << " " << random[i][2] << " ";
    outfile << closestPts[i][0] << " " << closestPts[i][1] << " " << closestPts[i][2] << " ";
    outfile << normalVecs[i][0] << " " << normalVecs[i][1] << " " << normalVecs[i][2] << " ";
    outfile << "0 0 0 0 0 0" << std::endl;
  }
  outfile.close();

  // Clean up.

  // Exit calls.
  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}