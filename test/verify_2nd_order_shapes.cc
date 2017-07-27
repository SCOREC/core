#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <apf.h>
#include <PCU.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#include <MeshSim.h>
#include <SimModel.h>
#endif
#include <stdlib.h>
#include <pcu_util.h>

int main(int argc, char** argv) {
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 2 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <mesh .smb>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
#ifdef HAVE_SIMMETRIX
  MS_init();
  SimModel_start();
  Sim_readLicenseFile(NULL);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_null();
  gmi_register_mesh();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::loadMdsMesh(g,argv[1]);

  int dim = m->getDimension();

  // get the first element dimension dim
  apf::MeshIterator* it;
  apf::MeshEntity* ent;
  it = m->begin(dim);
  ent = m->iterate(it);
  if (!ent) abort();
  m->end(it);

  apf::Mesh::Type type = m->getType(ent);

  apf::FieldShape* fs = m->getShape();
  apf::Element* elem = apf::createElement(m->getCoordinateField(), ent);


  // check order of the input mesh. Abort if not 2.
  int order = fs->getOrder();
  if (order != 2) {
    printf("Expected a 2nd-order input mesh. Aborting!\n");
    abort();
  }
  else
    printf("Mesh is 2nd-order.\n");

  // currently this only check for the shape of triangular
  // and quadrilateral elements.
  // TODO: add other types, in the future.
  switch (type) {
    case apf::Mesh::TRIANGLE:
      if (fs != apf::getLagrange(2)) {
      	printf("Shape of order 2 tris must be Lagrange(2) shape. Aborting!\n");
      	abort();
      }
      else
      	printf("Shape is Lagrange(2).\n");
      PCU_ALWAYS_ASSERT(apf::countNodes(elem) == 6);
      break;
    case apf::Mesh::QUAD:
      if (fs != apf::getSerendipity()) {
      	printf("Shape of order 2 quads must be serendipity shape. Aborting!\n");
      	abort();
      }
      else
      	printf("Shape is serendipity.\n");
      PCU_ALWAYS_ASSERT(apf::countNodes(elem) == 8);
      break;
    case apf::Mesh::TET:
      if (fs != apf::getLagrange(2)) {
      	printf("Shape of order 2 tets must be serendipity shape. Aborting!\n");
      	abort();
      }
      else
      	printf("Shape is Lagrange(2).\n");
      PCU_ALWAYS_ASSERT(apf::countNodes(elem) == 10);
      break;
    default:
      printf("Serify_shape is not implemented for "
      	  "elements of type %d. Aborting!", type);
      abort();
  }

  m->destroyNative();
  apf::destroyMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimModel_stop();
  MS_exit();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
