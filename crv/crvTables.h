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

/** \file crvTables.h
  * \brief main file for tables used in ordering curved entities */

namespace crv {

/** \brief table of indices for triangles, b2[order][i][j],
    only up to 10th order is stored, higher can be generated on the fly */
extern unsigned const* const* const b2[11];
/** \brief table of indices for tets, b3[order][i][j][k],
    only up to 4th order is stored, higher can be generated on the fly */
extern unsigned const* const* const* const b3[5];

/** \brief table of alignment used in alignSharedNodes,
    tet_tri[order][flip][rotate][node]; */
extern unsigned const* const* const* const tet_tri[7];

/** \brief parametric locations of midpoint nodes given a vertex number,
    elem_vert_xi[type][vertex_index] */
extern apf::Vector3 const* const elem_vert_xi[apf::Mesh::TYPES];
/** \brief parametric locations of midpoint nodes given an edge number,
    elem_edge_xi[type][edge_index] */
extern apf::Vector3 const* const elem_edge_xi[apf::Mesh::TYPES];

/** \brief table of edges of a triangle on a tet,
    ie the third triangle (index = 2) on a tet has edges
    1,5,4 of the tet */
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
/** \brief table of edge orientations in a triangle on a tet,
    corresponding to tet_tri_edges, 0 -> it is correctly oriented,
    1 -> it is flipped canonically */
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

/** \brief edge indices connected to a vertex of a tet, this does not
    comment on their orientation wrt to the vertex
    \details ordered as XJ Luo's thesis */
static int const vertEdges[4][3] = {{3,0,2},{0,4,1},{1,5,2},{3,5,4}};

/** \brief indices of opposite edges of an edge to a tet */
static int const oppEdges[6] = {5,3,4,1,2,0};

/** \brief given an edge of a tet, this is the index of
    the {left,right} triangle, based on direction of the edge */
static int const edgeFaces[6][2] = {{1,0},{2,0},{3,0},{3,1},{1,2},{2,3}};

}

#endif
