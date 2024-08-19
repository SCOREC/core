#include <PCU.h>
#include <queue>
#include <apf.h>
#include <pcu_util.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>


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

void checkParametrization(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb);

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  {
  auto PCUObj = std::unique_ptr<pcu::PCU>(new pcu::PCU(MPI_COMM_WORLD));

  if (argc != 2) {
    if(0==PCUObj.get()->Self())
      std::cerr << "usage: " << argv[0]
        << " <cre file .cre>\n";
    return EXIT_FAILURE;
  }

  const char* creFileName = argv[1];

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


  // check parametrization using capstone apis
  checkParametrization(m, g);

  }
  MPI_Finalize();
}

void checkParametrization(MeshDatabaseInterface* mdb, GeometryDatabaseInterface* gdb)
{
   MeshSmartIterator miter(mdb);
   mdb->get_topo_iterator(TOPO_VERTEX, miter);
   int count = 0;
   double sum = 0.0;
   for(mdb->iterator_begin(miter); !mdb->iterator_end(miter); mdb->iterator_next(miter)) {
     M_MTopo vert = mdb->iterator_value(miter);
     M_GTopo geom;
     GeometryTopoType gtype;
     mdb->get_geom_entity(vert, gtype, geom);
     if (!gdb->is_face(geom)) continue;
     double range_u[2];
     double range_v[2];
     gdb->get_parametrization_range(geom, 0, range_u[0], range_u[1]);
     gdb->get_parametrization_range(geom, 1, range_v[0], range_v[1]);
     GeometryTopoType gtype1;
     double u,v;
     mdb->get_vertex_uv_parameters(vert, u, v, gtype1);
     PCU_ALWAYS_ASSERT(gtype1 == gtype);

     // coordinate from mesh
     apf::Vector3 coord;
     mdb->get_vertex_coord(vert, &(coord[0]));

     // coordinate from surface
     vec3d x;
     gdb->get_point(geom, vec3d(u, v, 0.0), x);
     apf::Vector3 pcoord(x[0], x[1], x[2]);

     if (count < 50)
       printf("%d, %e, %e, %e, %e, %e, %e, %e\n", count, u, v, range_u[0], range_u[1], range_v[0], range_v[1], (coord-pcoord).getLength());
     sum += (coord-pcoord) * (coord-pcoord);
     count++;
   }
   printf("norm of the difference vector is %e\n", std::sqrt(sum));
}
