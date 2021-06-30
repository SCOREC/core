/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_CONVERT_H
#define APF_CONVERT_H

/** \file apfConvert.h
  \brief algorithms for mesh format conversion */

#include <map>
#include <vector>

namespace apf {

class Mesh;
class Mesh2;
class ModelEntity;
class MeshEntity;
using NewElements = std::vector<MeshEntity*>;

/** \brief convert one mesh data structure to another
  \details this function will fill in a structure that fully
  implements apf::Mesh2 by using information from an implementation
  of apf::Mesh. This is a fully scalable parallel mesh conversion
  tool. */
void convert(Mesh *in, Mesh2 *out,
             MeshEntity** nodes=NULL, MeshEntity** elems=NULL, bool copy_data=true);

/** \brief a map from global ids to vertex objects */
typedef std::map<int, MeshEntity*> GlobalToVert;

/** \brief assemble a mixed-cell-type mesh from just a connectivity array
  \details construct is now split into two functions, 
  assemble and finalise. The premise of assemble being 
  that it is called multiple times for a given cell type,
  across several different cell types in the input mesh. */
NewElements assemble(Mesh2* m, const int* conn, int nelem, int etype,
    GlobalToVert& globalToVert);

/** \brief finalise construction of a mixed-cell-type mesh from just a connectivity array
  \details construct is now split into two functions, 
  assemble and finalise. Once the mixed cell type mesh 
  is assembled finalise should be called.  Doing it this 
  way provides non-breaking changes for current users of 
  construct, which now just calls assemble and finalise. */
void finalise(Mesh2* m, GlobalToVert& globalToVert);

/** \brief construct a mesh from just a connectivity array
  \details this function is here to interface with very
  simple mesh formats. Given a set of elements described
  only in terms of the ordered global ids of their vertices,
  this function builds a reasonable apf::Mesh2 structure
  and as a side effect returns a map from global ids
  to local vertices.  This functions assumes a uniform 
  cell type.  Use a combination of assemble and finalise for 
  meshes loaded with mixed cell types.

  This is a fully scalable parallel mesh construction
  algorithm, no processor incurs memory or runtime costs
  proportional to the global mesh size.

  Note that all vertices will have zero coordinates, so
  it is often good to use apf::setCoords after this. */
NewElements construct(Mesh2* m, const int* conn, int nelem, int etype,
    GlobalToVert& globalToVert);

/** \brief Assign coordinates to the mesh
  * \details
  * Each peer provides a set of the coordinates. The coords most be ordered
  * according to the global ids of the vertices. Peer 0 provides the coords
  * for vertices 0 to m-1, peer to for m to n-1, ...
  * After this call, all vertices in the apf::Mesh2 object have correct
  * coordinates assigned.
  */
void setCoords(Mesh2* m, const double* coords, int nverts,
    GlobalToVert& globalToVert);

/** \brief convert an apf::Mesh2 object into a connectivity array
  \details this is useful for debugging the apf::convert function
  \param mesh the apf mesh
  \param nelem number of elements
  \param etype apf::Mesh::Type
  \param cellDim dimension of elements (if embedded in a higher dimension manifold)
  */
void destruct(Mesh2* m, int*& conn, int& nelem, int &etype, int cellDim = -1);

/** \brief get a contiguous set of global vertex coordinates
  \details this is used for debugging apf::setCoords */
void extractCoords(Mesh2* m, double*& coords, int& nverts);

}

#endif
