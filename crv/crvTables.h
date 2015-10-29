/*
 * Copyright 2015 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef CRVTABLES_H
#define CRVTABLES_H

namespace crv {

extern unsigned const* const* const b2[11];
extern unsigned const* const* const* const b3[5];

extern unsigned const* const* const* const tet_tri[7];

enum {
  BEZIER,
  GREGORY,
  TYPES
};

// negative -> flipped relative to canonical
// relies on e0 being always ordered correctly
static int const tet_tri_edges[4][3] =
{{0,1,2},{0,4,3},{1,5,4},{2,5,3}};
static bool const flip_tet_tri_edges[4][3] =
{{0,0,0},{0,0,1},{0,0,1},{1,0,1}};

}

#endif
