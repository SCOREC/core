/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVTABLES_H
#define CRVTABLES_H

#include "apf.h"
#include "apfMesh.h"

namespace crv {

extern unsigned const* const* const b2[11];
extern unsigned const* const* const* const b3[5];

extern unsigned const* const* const* const tet_tri[7];

extern apf::Vector3 const* const elem_vert_xi[apf::Mesh::TYPES];

extern apf::Vector3 const* const elem_edge_xi[apf::Mesh::TYPES];

// negative -> flipped relative to canonical
// relies on e0 being always ordered correctly
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

// edge indices connected to a vert, ordered as XJ Luo's thesis
static int const vertEdges[4][3] = {{3,0,2},{0,4,1},{1,5,2},{3,5,4}};

// opposite edges (edge with other two verts) for a tet
static int const oppEdges[6] = {5,3,4,1,2,0};

// for each edge, this has the left/right face
static int const edgeFaces[6][2] = {{1,0},{2,0},{3,0},{3,1},{1,2},{2,3}};

}

#endif
