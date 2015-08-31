/******************************************************************************

  Copyright 2015 Scientific Computation Research Center,
      Rensselaer Polytechnic Institute. All rights reserved.

  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.

*******************************************************************************/

#ifndef CRVTABLES_H
#define CRVTABLES_H

namespace crv {

extern unsigned const* const* const b2[10];
extern unsigned const* const* const* const b3[5];

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

// numbers of nodes on
static int const curved_face_internal[2][6] =
{{0,0,1,3,6,10},{0,0,6,6,0,0}};

static int const curved_tet_internal[2][6] =
{{0,0,0,1,4,10},{0,0,0,1,4,10}};

// total numbers of nodes
static int const curved_face_total[2][6] =
{{3,6,10,15,21,28},{0,0,15,18,0,0}};

static int const blended_tet_total[2][6] =
{{4,10,20,34,52,74},{0,0,40,46,0,0}};

static int const curved_tet_total[2][6] =
{{4,10,20,35,56,84},{0,0,40,47,0,0}};

}

#endif
