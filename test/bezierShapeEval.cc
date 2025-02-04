#include <crv.h>
#include <crvBezier.h>
#include <crvTables.h>
#include <crvSnap.h>
#include <crvMath.h>
#include <crvBezierShapes.h>
#include <gmi_analytic.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <apf.h>
#include <apfShape.h>
#include <lionPrint.h>
#include <mth.h>
#include <mth_def.h>
#include <pcu_util.h>
#include <ostream>
/* This file contains miscellaneous tests relating to bezier shapes and
 * blended bezier shapes
 */

static apf::Mesh2* makeOneTriMesh(int order, apf::MeshEntity* &ent, pcu::PCU *PCUObj);
static apf::Mesh2* makeOneTetMesh(int order, apf::MeshEntity* &ent, pcu::PCU *PCUObj);


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  {
  pcu::PCU PCUObj = pcu::PCU(MPI_COMM_WORLD);
  lion_set_verbosity(1);
  if ( argc != 7 ) {
    if ( !PCUObj.Self() ) {
      printf("Usage: %s <ent_type> <order> <blend_order> <xi_0> <xi_1> <xi_2>>\n", argv[0]);
      printf("<ent_type>            can only be 2 (for triangles) and 4 (for tets)\n");
      printf("<order>               is the order of bezier\n");
      printf("<blend_order>         can be -1, 0, 1, 2 (-1 means no blending)\n");
      printf("<xi_0> <xi_1> <xi_2>  inquiry point in the parent entity)\n");
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  int type  = atoi(argv[1]);
  int order = atoi(argv[2]);
  int b     = atoi(argv[3]);
  PCU_ALWAYS_ASSERT_VERBOSE(type == 2 || type == 4,
      "<ent_type> can only be 2 or 4!");
  PCU_ALWAYS_ASSERT_VERBOSE(b > -2 &&  b < 3,
      "<blend_order> must be between -1 and 2!");

  apf::MeshEntity* ent = 0;
  apf::Mesh2* m = type == 2 ? makeOneTriMesh(order,ent,&PCUObj) : makeOneTetMesh(order,ent,&PCUObj);

  double xi0 = atof(argv[4]);
  double xi1 = atof(argv[5]);
  double xi2 = atof(argv[6]);


  // set the order
  apf::FieldShape* bezierShape = crv::getBezier(order);
  int non = bezierShape->getEntityShape(type)->countNodes();
  apf::Vector3 xi(xi0, xi1, xi2);
  apf::NewArray<double> vals(non);

  if (b == -1) // regular shape functions
     bezierShape->getEntityShape(type)->getValues(m,ent,xi,vals);
  else // blended shape functions
  {
    crv::setBlendingOrder(type, b);
    if (type == 2)
      crv::BlendedTriangleGetValues(m,ent,xi,vals);
    else
      crv::BlendedTetGetValues(m,ent,xi,vals);
  }


  printf("shape values are \n");
  for (int i = 0; i < non; i++) {
    printf("%d,  %E \n", i, vals[i]);
  }

  }
  MPI_Finalize();
}

static apf::Mesh2* makeOneTriMesh(int order, apf::MeshEntity* &ent, pcu::PCU *PCUObj)
{
  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(0, 2, false, PCUObj);

  double vert_coords[3][6] = {
      {0.,0.,0., 0., 0., 0.},
      {1.,0.,0., 0., 0., 0.},
      {0.,1.,0., 0., 0., 0.},
  };

  // each edge is defined by the bounding verts
  int edge_info[3][2] = {
      {0,1},
      {1,2},
      {2,0}
  };

  apf::MeshEntity* verts[3];
  apf::MeshEntity* edges[3];

  for (int i = 0; i < 3; i++) {
    apf::Vector3 coords(vert_coords[i][0],
                        vert_coords[i][1],
                        vert_coords[i][2]);
    apf::Vector3 params(vert_coords[i][3],
                        vert_coords[i][4],
                        vert_coords[i][5]);
    verts[i] = mesh->createVertex(0, coords, params);
  }
  for (int i = 0; i < 3; i++) {
    apf::MeshEntity* down_vs[2] = {verts[edge_info[i][0]],
                                   verts[edge_info[i][1]]};
    edges[i] = mesh->createEntity(apf::Mesh::EDGE, 0, down_vs);
  }

  ent = mesh->createEntity(apf::Mesh::TRIANGLE, 0, edges);

  mesh->acceptChanges();

  apf::changeMeshShape(mesh, crv::getBezier(order),true);

  printf ("Made tri with verts:\n");
  for (int i = 0; i < 3; i++) {
    printf("v0: (%e,%e,%e)\n", vert_coords[i][0], vert_coords[i][1], vert_coords[i][2]);
  }

  printf ("all nodes of the tri:\n");
  apf::Element* elem = apf::createElement(mesh->getCoordinateField(),ent);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);
  apf::destroyElement(elem);
  for (int i=0; i<(int)nodes.size(); i++)
  {
    printf("node %d: (%e,%e,%e)\n", i, nodes[i][0], nodes[i][1], nodes[i][2]);
  }
  return mesh;
}

static apf::Mesh2* makeOneTetMesh(int order, apf::MeshEntity* &ent, pcu::PCU *PCUObj)
{

  apf::Mesh2* mesh = apf::makeEmptyMdsMesh(0, 3, false, PCUObj);

  double vert_coords[4][6] = {
      {0.,0.,0., 0., 0., 0.},
      {1.,0.,0., 0., 0., 0.},
      {0.,1.,0., 0., 0., 0.},
      {0.,0.,1., 0., 0., 0.},
  };

  // each edge is defined by the bounding verts
  int const (*edge_info)[2] = apf::tet_edge_verts;
  int face_info[4][3] = {
    {0,1,2},
    {0,4,3},
    {1,5,4},
    {2,5,3}
  };


  apf::MeshEntity* verts[4];
  apf::MeshEntity* edges[6];
  apf::MeshEntity* faces[4];

  for (int i = 0; i < 4; i++) {
    apf::Vector3 coords(vert_coords[i][0],
                        vert_coords[i][1],
                        vert_coords[i][2]);
    apf::Vector3 params(vert_coords[i][3],
                        vert_coords[i][4],
                        vert_coords[i][5]);
    verts[i] = mesh->createVertex(0, coords, params);
  }
  for (int i = 0; i < 6; i++) {
    apf::MeshEntity* down_vs[2] = {verts[edge_info[i][0]],
                                   verts[edge_info[i][1]]};
    edges[i] = mesh->createEntity(apf::Mesh::EDGE, 0, down_vs);
  }

  for (int i = 0; i < 4; i++) {
    apf::MeshEntity* down_es[3] = {edges[face_info[i][0]],
                                   edges[face_info[i][1]],
                                   edges[face_info[i][2]]};
    faces[i] = mesh->createEntity(apf::Mesh::TRIANGLE, 0, down_es);
  }

  ent = mesh->createEntity(apf::Mesh::TET, 0, faces);

  mesh->acceptChanges();

  apf::changeMeshShape(mesh, crv::getBezier(order),true);

  printf ("Made tet with verts:\n");
  for (int i = 0; i < 4; i++) {
    printf("v0: (%e,%e,%e)\n", vert_coords[i][0], vert_coords[i][1], vert_coords[i][2]);
  }

  printf ("all nodes of the tet:\n");
  apf::Element* elem = apf::createElement(mesh->getCoordinateField(),ent);
  apf::NewArray<apf::Vector3> nodes;
  apf::getVectorNodes(elem,nodes);
  apf::destroyElement(elem);
  for (int i=0; i<(int)nodes.size(); i++)
  {
    printf("node %d: (%e,%e,%e)\n", i, nodes[i][0], nodes[i][1], nodes[i][2]);
  }
  return mesh;
}
