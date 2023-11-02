#include <PCU.h>
#include <apfCAP.h>
#include <apfMDS.h>
#include <samSz.h>
#include <queue>
#include <ma.h>
#include <gmi.h>
#include <gmi_cap.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apf.h>
#include <apfConvert.h>
#include <apfMesh2.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <ma.h>
#include <pcu_util.h>
#include <lionPrint.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>


#include "CapstoneModule.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Analysis.h"
#include "CreateMG_Framework_Application.h"
#include "CreateMG_Framework_Attributes.h"
#include "CreateMG_Framework_Core.h"
#include "CreateMG_Framework_Geometry.h"
#include "CreateMG_Framework_Mesh.h"

using namespace CreateMG;
using namespace CreateMG::Attribution;
using namespace CreateMG::Mesh;
using namespace CreateMG::Geometry;


double totalArea(apf::Mesh2* m);

double addSizeDistribution(apf::Mesh2* m, int factor,
    std::vector<int>& binCount, std::vector<double>& binArea);

void edgeLengthStats(apf::Mesh2* m);

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  double initialTime = PCU_Time();

  if (argc != 3) {
    if(0==PCU_Comm_Self())
      std::cerr << "usage: " << argv[0]
        << " <cre file .cre> <factor>\n";
    return EXIT_FAILURE;
  }


  gmi_register_mesh();
  gmi_register_null();

  const char* creFileName = argv[1];
  int factor = atoi(argv[2]);

  // load capstone mesh
  // create an instance of the Capstone Module activating CREATE/CREATE/CREATE
  // for the Geometry/Mesh/Attribution databases
  /* const std::string gdbName("Geometry Database : Create");// Switch Create with SMLIB for CAD */
  const std::string gdbName("Geometry Database : SMLIB");// Switch Create with SMLIB for CAD
  const std::string mdbName("Mesh Database : Create");
  const std::string adbName("Attribution Database : Create");

  CapstoneModule  cs("the_module", gdbName.c_str(), mdbName.c_str(), adbName.c_str());

  GeometryDatabaseInterface     *g = cs.get_geometry();
  MeshDatabaseInterface         *m = cs.get_mesh();
  AppContext                    *c = cs.get_context();


  PCU_ALWAYS_ASSERT(g);
  PCU_ALWAYS_ASSERT(m);
  PCU_ALWAYS_ASSERT(c);

  v_string filenames;
  filenames.push_back(creFileName);

  M_GModel gmodel = cs.load_files(filenames);

  int numbreps = 0;
  MG_CALL(g->get_num_breps(numbreps));
  std::cout << "number of b reps is " << numbreps << std::endl;
  if(numbreps == 0)
      error(HERE, ERR_INVALID_INPUT, "Model is empty");

  M_MModel mmodel;
  // Pick the volume mesh model from associated mesh models to this geom model
  std::vector<M_MModel> mmodels;
  MG_API_CALL(m, get_associated_mesh_models(gmodel, mmodels));
  for(std::size_t i = 0; i < mmodels.size(); ++i)
  {
      M_MModel ammodel = mmodels[i];
      std::size_t numregs = 0;
      std::size_t numfaces = 0;
      std::size_t numedges = 0;
      std::size_t numverts = 0;
      MG_API_CALL(m, set_current_model(ammodel));
      MG_API_CALL(m, get_num_topos(TOPO_REGION, numregs));
      MG_API_CALL(m, get_num_topos(TOPO_FACE, numfaces));
      MG_API_CALL(m, get_num_topos(TOPO_EDGE, numedges));
      MG_API_CALL(m, get_num_topos(TOPO_VERTEX, numverts));
      std::cout << "num regions is " << numregs << std::endl;
      std::cout << "num faces   is " << numfaces << std::endl;
      std::cout << "num edges   is " << numedges << std::endl;
      std::cout << "num verts   is " << numverts << std::endl;
      std::cout << "-----------" << std::endl;
      if(numregs > 0)
      {
	  mmodel = ammodel;
	  break;
      }
  }

  /* SET THE ADJACENCIES */
  MG_API_CALL(m, set_adjacency_state(REGION2FACE|
                                     REGION2EDGE|
                                     REGION2VERTEX|
                                     FACE2EDGE|
                                     FACE2VERTEX));
  MG_API_CALL(m, set_reverse_states());
  MG_API_CALL(m, set_adjacency_scope(TOPO_EDGE, SCOPE_FULL));
  MG_API_CALL(m, set_adjacency_scope(TOPO_FACE, SCOPE_FULL));
  MG_API_CALL(m, compute_adjacency());


  gmi_cap_start();
  gmi_register_cap();

  printf("---- CapStone Mesh Loaded. \n");

  // create the mesh object (this one is CapStone underneath)
  printf("\n---- Creating Mesh Wrapper Object. \n");
  apf::Mesh2* mesh = apf::createMesh(m,g);
  printf("---- Creating Mesh Wrapper Object: Done. \n");

  // add size distribution based on area
  std::vector<int> binCount;
  std::vector<double> binArea;
  edgeLengthStats(mesh);

  gmi_cap_stop();
  PCU_Comm_Free();
  MPI_Finalize();
}

double totalArea(apf::Mesh2* m)
{
  double area = 0.0;
  apf::MeshEntity* e;
  apf::MeshIterator* it;
  it = m->begin(2);
  while ( (e = m->iterate(it)) )
    area += apf::measure(m, e);
  m->end(it);
  return area;
}

double  addSizeDistribution(apf::Mesh2* m, int factor,
    std::vector<int>& binCount, std::vector<double>& binArea)
{
  // find the min length edge
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  double minLength = 1.0e12;
  it = m->begin(1);
  while ( (e = m->iterate(it)) ) {
    double l = apf::measure(m, e);
    if (l < minLength)
      minLength = l;
  }
  m->end(it);

  double minArea = minLength * minLength / 2.;

  // initialize the arrays
  for (int i = 1; i <= factor*factor; i*=4) {
    binCount.push_back(0);
    binArea.push_back(0.0);
  }

  binCount.push_back(0);
  binCount.push_back(0);
  binArea.push_back(0.0);
  binArea.push_back(0.0);

  it = m->begin(2);
  while ( (e = m->iterate(it)) ) {
    double a = apf::measure(m, e);
    if (a < minArea){
      binCount[0] += 1;
      binArea[0] += a;
      continue;
    }
    int count = 1;
    int i;
    for (i = 1; i <= factor*factor; i*=4) {
      if ((a >= i * minArea) && (a < 4 * i * minArea)) {
      	binCount[count] += 1;
      	binArea[count] += a;
      }
      count++;
    }
    i = i/4;
    if (a >= 4 * i * minArea) {
      binCount[count] += 1;
      binArea[count] += a;
    }
  }
  m->end(it);
  return minArea;
}

void  edgeLengthStats(apf::Mesh2* m)
{
  // find the min length edge
  apf::MeshEntity* e;
  apf::MeshIterator* it;

  double minLength = 1.0e12;
  double maxLength = 0.0;
  double avgLength = 0.0;
  double count = 0;
  it = m->begin(1);
  while ( (e = m->iterate(it)) ) {
    double l = apf::measure(m, e);
    if (l < minLength)
      minLength = l;
    if (l > maxLength)
      maxLength = l;
    avgLength += l;
    count++;
  }
  m->end(it);
  avgLength /= count;

  printf("min/max/avg lengths are: %E/%E/%E\n", minLength, maxLength, avgLength);
  printf("num tris is %d\n", m->count(2));
}
