/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#ifndef MA_TABLES_H
#define MA_TABLES_H

#include "maMesh.h"

namespace ma {

struct CodeMatch
{
  int rotation; //the 'n' argument to the rotateEntity functions
  int code_index; //the index of the template to be applied
};

extern CodeMatch const tri_code_match[(1<<3)];
extern CodeMatch const quad_code_match[(1<<4)];
extern CodeMatch const tet_code_match[(1<<6)];
extern CodeMatch const prism_code_match[(1<<9)];
extern CodeMatch const pyramid_code_match[(1<<8)];
extern CodeMatch const* code_match[TYPES];

/* this table defines the mapping from new to old
   vertex indices for one of the 12 rotations. */
extern int const tet_rotation[12][4];
extern int const prism_rotation[6][6];
extern int const pyramid_rotation[4][5];

enum
{
  edge_edge_code_count = 2,
  tri_edge_code_count = 4,
  quad_edge_code_count = 3,
  tet_edge_code_count = 12,
  prism_edge_code_count = 5,
  pyramid_edge_code_count = 17
};

extern int const edge_edge_codes[edge_edge_code_count];
extern int const tri_edge_codes[tri_edge_code_count];
extern int const quad_edge_codes[quad_edge_code_count];
extern int const tet_edge_codes[tet_edge_code_count];
extern int const prism_edge_codes[prism_edge_code_count];
extern int const pyramid_edge_codes[pyramid_edge_code_count];

extern int const prism_diag_match[(1<<3)];
extern int const prism_diag_choices[4];

/* octagon vertex encoding is:
   bottom vertex
   4 middle vertices in order curl -> top
   top vertex */

extern int const oct_rotation[24][6];

}

#endif
