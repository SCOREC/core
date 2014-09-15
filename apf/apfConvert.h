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

namespace apf {

class Mesh;
class Mesh2;
class ModelEntity;
class MeshEntity;

/** \brief convert one mesh data structure to another
  \details this function will fill in a structure that fully
  implements apf::Mesh2 by using information from an implementation
  of apf::Mesh. This is a fully scalable parallel mesh conversion
  tool. */
void convert(Mesh *in, Mesh2 *out);

/** \brief a map from global ids to vertex objects */
typedef std::map<int, MeshEntity*> GlobalToVert;

/** \brief construct a mesh from just a connectivity array
  \details this function is here to interface with very
  simple mesh formats. Given a set of elements described
  only in terms of the ordered global ids of their vertices,
  this function builds a reasonable apf::Mesh2 structure
  and as a side effect returns a map from global ids
  to local vertices.

  This is a fully scalable parallel mesh construction
  algorithm, no processor incurs memory or runtime costs
  proportional to the global mesh size.

  Note that all vertices will have zero coordinates, so
  it is often good to use apf::setCoords after this. */
void construct(Mesh2* m, int* conn, int nelem, int etype,
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
  \details this is useful for debugging the apf::convert function */
void destruct(Mesh2* m, int*& conn, int& nelem, int &etype);

/** \brief get a contiguous set of global vertex coordinates
  \details this is used for debugging apf::setCoords */
void extractCoords(Mesh2* m, double*& coords, int& nverts);

}

#endif
