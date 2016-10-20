#include "pumi.h"

#include <apf.h>
#include <cstring>
#include <mpi.h>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include "apfMDS.h"
#include "apfShape.h"

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int serial=0;

void getConfig(int argc, char** argv)
{
  if (argc < 4) {
    if (!pumi_rank() )
      printf("Usage: %s <model> <mesh> <outMesh>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  if (argc > 4) 
    serial=atoi(argv[4]);
}

Migration* get_xgc_plan(pGeom g, pMesh m)
{
  int dim = pumi_mesh_getDim(m);
  Migration* plan = new Migration(m);
  if (!pumi_rank()) return plan;

  pMeshEnt e;
  int num_gface = pumi_geom_getNumEnt(g, dim);
  assert(num_gface==pumi_size());
  int gface_id;
  int dest_pid;
  pMeshIter it = m->begin(2); // face
  while ((e = m->iterate(it))) 
  { 
    pGeomEnt gface = pumi_ment_getGeomClas(e); // get the classification
    gface_id = pumi_gent_getID(gface); // get the geom face id
    dest_pid = gface_id-1;
    plan->send(e, dest_pid);    
  }
  m->end(it);
  return plan;
}

#include <algorithm>
//*******************************************************
int xgc_mesh_search(pMesh m, int* initial_simplex,
		    double* final_position,
		    int* final_simplex)
//*******************************************************
{
  double begin_time=0.0, begin_mem=0.0;
  bool located = false;
  int simplex_dim = m->getDimension();
  int edge_dim = 1;
  int bridge_dim = simplex_dim - 1;
  int node = 0, edge_type = 1;    // Mesh entity type; see APF documentation.

  pMeshEnt simplex = apf::getMdsEntity(m, simplex_dim, *initial_simplex);
  int simplex_index = *initial_simplex;
  pMeshEnt e = NULL;
  apf::Adjacent adjacent;
  int edge_curr_index, edge_prev_index=0;

  int count = 0, max_count = m->count(simplex_dim);
  double tol = 1e-15;

  apf::Matrix3x3 coords;
  apf::Field *fc = m->findField("coordinates");

  Downward vertices;
  apf::Vector3 xy[3], b_coords;
  std::vector<int> b_indices(3);
  double x0, x1, x2, y0, y1, y2;
  pMeshEnt edge_vertices[2];

  while (not located) {
    // Read the coordinates of the vertices of the simplex to compute the
    // barycentric coordinates.
    apf::getMatrix(fc, simplex, node, coords);
    x0 = coords[0][0]; x1 = coords[1][0]; x2 = coords[2][0];
    y0 = coords[0][1]; y1 = coords[1][1]; y2 = coords[2][1];
    b_coords[1] = ((x0 * (final_position[1] - y2) +
		    final_position[0] * (-y0 + y2) +
		    x2 * (y0 - final_position[1]))/
		   (x0 * (y1 - y2) + x1 * (-y0 + y2) +
		    x2 * (y0 - y1)));
    b_coords[2] = ((y0 * (-final_position[0] + x1) +
		    y1 * (final_position[0] - x0) +
		    final_position[1] * (x0 - x1))/
		   (y0 * (x1 - x2) + y1 * (-x0 + x2) +
		    y2 * (x0 - x1)));
    b_coords[0] = 1 - b_coords[1] - b_coords[2];

    // If all positive for current simplex, exit.
    if ((b_coords[0] >= -tol) && (b_coords[1] >= -tol) &&
	(b_coords[2] >= -tol)) {
      located = true;
      std::cout<<"("<<pumi_rank()<<") "<<__func__<<" found for simplex "<<* initial_simplex<<" >> "<<simplex_index<<"\n";
      pumi_printTimeMem("\n[adj_search] elapsed time and increased heap memory:", pumi_getTime()-begin_time, pumi_getMem()-begin_mem);
      *final_simplex = simplex_index;
      return 0;
    }

    /* Obtain the index of most negative barycentric coordinate and determine
    edge opposite to this vertex. This uses a fancy lambda function solution for 
    argument sorting the values in the vector b_coords.
    */

    b_indices[0] = 0; b_indices[1] = 1; b_indices[2] = 2;
    std::sort(b_indices.begin(), b_indices.end(),
	      [&](int i1, int i2)
	      {return b_coords[i1] < b_coords[i2];});     
    edge_vertices[0] = vertices[b_indices[1]];
    edge_vertices[1] = vertices[b_indices[2]];
    
    // Find neighboring simplex sharing this edge
    e = apf::findElement(m, edge_type, edge_vertices);
    edge_curr_index = apf::getMdsIndex(m, e);

    // If current edge choice is same as previous edge, pick edge
    // opposite to second least (actual, not absolute, valued) barycentric
    // coordinate.
    /*
    if (edge_curr_index == edge_prev_index) {
      edge_count = 0;
      for (int j = 0; j < 3; ++j)
	if (j != bneg_index[1]) {
	  edge_vertices[edge_count] = vertices[j];
	  ++edge_count;
	}
      e = apf::findElement(m, edge_type, edge_vertices);
      edge_curr_index = apf::getMdsIndex(m, e);
    }     
    */
    apf::getBridgeAdjacent(m, e,
			   bridge_dim, simplex_dim,
			   adjacent);
    
    if (adjacent.getSize() == 2) {
      for (size_t j = 0; j < adjacent.getSize(); ++j) {
	int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
	if (new_simplex_index != simplex_index) {
	  simplex = adjacent[j];
	  simplex_index = new_simplex_index;
	  break;
	}
      }
    }
    else {
      Downward edges;
      int ne = m->getDownward(simplex, edge_dim, edges);
      for (int j = 0; j < ne; ++j) {
	int edge_tmp_index = apf::getMdsIndex(m, edges[j]);
	if ((edge_tmp_index != edge_curr_index) &&
	    (edge_tmp_index != edge_prev_index)) {
	  e = edges[j];
	  break;
	}
      }
      edge_curr_index = apf::getMdsIndex(m, e);    
      apf::getBridgeAdjacent(m, e,
			     bridge_dim, simplex_dim,
			     adjacent);
      for (size_t j = 0; j < adjacent.getSize(); ++j) 
      {
	int new_simplex_index = apf::getMdsIndex(m, adjacent[j]);
	if (new_simplex_index != simplex_index)
        {
	  simplex = adjacent[j];
//	  simplex_index = new_simplex_index;
//	  break;
        }
      }
    }
    
    // Keep track of edge via which we entered the current simplex
    edge_prev_index = edge_curr_index;
    ++count;

    if (count == max_count)
    {
      std::cout<<"("<<pumi_rank()<<") "<<__func__<<" failed for simplex "<<* initial_simplex<<"\n";
      *final_simplex = -2;
      return 1;
    }
  }
  return 0;
}

void xgc_setup_matrix(pMesh m)
{
  // Added by kk
  // Note: The value 1 below is the enumerator key for apf::Vector
  int simplex_dim = m->getDimension();
  int simplex_count = m->count(simplex_dim);
  pMeshEnt simplex = NULL;
  apf::Downward vertices;
  apf::Vector3 xyz[3];
  int nv, vertex_dim = 0;
  int node = 0;

  apf::Field *fc = apf::createField(m,
				    "coordinates", apf::MATRIX,
				    apf::getConstant(simplex_dim));
  for (int i = 0; i < simplex_count; ++i) {
    simplex = apf::getMdsEntity(m, simplex_dim, i);
    assert(simplex);
    nv = m->getDownward(simplex, vertex_dim, vertices);
    for (int j = 0; j < nv; ++j) {
      m->getPoint(vertices[j], 0, xyz[j]);
      // Note: second argument is 0 for linear meshes
    }
    apf::Matrix3x3 coords(xyz[0][0], xyz[0][1], xyz[0][2],
			  xyz[1][0], xyz[1][1], xyz[1][2],
			  xyz[2][0], xyz[2][1], xyz[2][2]);
    apf::setMatrix(fc, simplex, node, coords);
  }
}


int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  pumi_start();

  getConfig(argc,argv);

  pGeom g = pumi_geom_load(modelFile);
  pMesh m;
  if (serial) 
  {
    m = pumi_mesh_loadSerial(g, meshFile);
    // split a serial mesh based on model ID
    Migration* plan = get_xgc_plan(g, m);
    pumi_mesh_migrate(m, plan);
    pumi_mesh_write(m, outFile);
  }
  else 
    m = pumi_mesh_load(g, meshFile, pumi_size());

  xgc_setup_matrix(m);
  int mesh_dim=pumi_mesh_getDim(m);
  int max_count=(int)m->count(mesh_dim);
  for (int i=0; i<max_count; ++i)
  {
    double coords[3];
    int final_simplex;
    xgc_mesh_search(m, &i, coords, &final_simplex);
  }
  
  // write to vtk
  char without_extension[256];
  snprintf(without_extension,strlen(argv[3])-3,"%s",argv[3]);

  char vtk_fname[32];
  sprintf(vtk_fname,"%s",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");

  // ghosting
  pumi_ghost_createLayer(m, 0, 2, 3, 0);
  sprintf(vtk_fname,"%s-ghosted",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");
  pumi_ghost_delete(m);

  pumi_ghost_createLayer(m, 0, 2, 3, 1);
  sprintf(vtk_fname,"%s-ghosted-copy",without_extension); 
  pumi_mesh_write(m, vtk_fname, "vtk");
  pumi_ghost_delete(m);

  pumi_mesh_delete(m);

  pumi_finalize();
  MPI_Finalize();
}

