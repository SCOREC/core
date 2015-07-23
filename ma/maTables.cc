/****************************************************************************** 

  Copyright 2013 Scientific Computation Research Center, 
      Rensselaer Polytechnic Institute. All rights reserved.
  
  The LICENSE file included with this distribution describes the terms
  of the SCOREC Non-Commercial License this program is distributed under.
 
*******************************************************************************/
#include "maTables.h"

namespace ma {

int const edge_edge_codes[edge_edge_code_count] =
{0x0//not split
,0x1//split
};

CodeMatch const edge_code_match[2] =
{{0,0}//not split
,{0,1}//split
};

int const tri_edge_codes[tri_edge_code_count] =
{0x0//000, no split
,0x1//001, one edge split
,0x3//011, two edges split
,0x7//111, three edges split
};

/* matches for triangle edge split
   configurations. */

CodeMatch const tri_code_match[(1<<3)] =
{{0,0}//000
,{0,1}//001
,{1,1}//010
,{0,2}//011
,{2,1}//100
,{2,2}//101
,{1,2}//110
,{0,3}//111
};

/* rotates a tet to one of 12 possible positions.
   In all cases, refer to the tetrahedral edge and
   vertex numbering in the PUMI Users Guide or the
   Simmetrix documentation (they should be equivalent).
   The position is identified by which vertex is rotated
   into the vertex 0 position, and which rotation the
   remaining 3 vertices are in. */

int const tet_rotation[12][4] =
{{0,1,2,3}//0
,{0,2,3,1}//1
,{0,3,1,2}//2
,{1,0,3,2}//3
,{1,3,2,0}//4
,{1,2,0,3}//5
,{2,0,1,3}//6
,{2,1,3,0}//7
,{2,3,0,1}//8
,{3,0,2,1}//9
,{3,2,1,0}//10
,{3,1,0,2}//11
};

/* these are the canonical edge split
   configurations that form the first
   layer of filtering for tet refinement
   templates.
   All other combinations of split/non-split
   edges can be rotated into one of
   these configurations and dispatched
   to the corresponding template. */

int const tet_edge_codes[tet_edge_code_count] =
{0x00//000000, no split
,0x01//000001, one edge split
,0x06//000110, two adjacent edges split
,0x21//100001, two opposite edges split
,0x07//000111, three edges of one face split
,0x25//100101, three splits, two faces have two splits
,0x23//100011, three splits, two faces have two splits, variant 2
,0x38//111000, three adjacent edges split
,0x27//100111, four splits, three on one face
,0x1E//011110, four splits, all faces have two
,0x1F//011111, five splits
,0x3F//111111, six splits
};

/* note that the "variant 2" above is not mentioned
   in the literature as being distinct from the other,
   but it can be seen that it requires more than just
   rotation to go from one to the other */

/* the following table was automatically generated !
   dont touch it ! just use the tetCodeMatch.cc program */
CodeMatch const tet_code_match[(1<<6)] =
{{ 0, 0}
,{ 0, 1}
,{ 5, 1}
,{ 6, 2}
,{ 1, 1}
,{ 5, 2}
,{ 0, 2}
,{ 0, 4}
,{ 2, 1}
,{11, 2}
,{ 2, 3}
,{ 2, 5}
,{ 8, 2}
,{ 4, 7}
,{ 5, 6}
,{ 5, 8}
,{ 4, 1}
,{ 2, 2}
,{10, 2}
,{ 1, 7}
,{ 1, 3}
,{ 6, 6}
,{ 4, 5}
,{ 6, 8}
,{ 3, 2}
,{ 2, 4}
,{ 2, 6}
,{ 2, 8}
,{ 1, 5}
,{11, 8}
,{ 0, 9}
,{ 0,10}
,{ 8, 1}
,{ 0, 3}
,{ 4, 2}
,{ 0, 6}
,{ 9, 2}
,{ 0, 5}
,{ 2, 7}
,{ 0, 8}
,{ 1, 2}
,{ 3, 6}
,{ 7, 5}
,{ 1, 9}
,{ 1, 4}
,{ 8, 8}
,{ 9, 8}
,{ 1,10}
,{ 7, 2}
,{ 3, 5}
,{ 4, 4}
,{10, 8}
,{ 1, 6}
,{ 2, 9}
,{ 4, 8}
,{ 5,10}
,{ 0, 7}
,{ 3, 8}
,{ 7, 8}
,{ 4,10}
,{ 1, 8}
,{ 2,10}
,{ 8,10}
,{ 0,11}
};

CodeMatch const* code_match[apf::Mesh::TYPES] =
{0,//vertex
 edge_code_match,
 tri_code_match,
 quad_code_match,
 tet_code_match,
 0,//hex
 prism_code_match,
 pyramid_code_match
};

/* prism rotations: there are 6 possible,
   since once one vertex is mapped the others
   are fixed.
   These rotations follow which vertex maps
   to new vertex 0.
   The first 0-2 rotations spin around the axis,
   rotation 3 flips it upside-down, and rotations 4-5
   spin it upside-down */

int const prism_rotation[6][6] =
{{0,1,2,3,4,5}
,{1,2,0,4,5,3}
,{2,0,1,5,3,4}
,{3,5,4,0,2,1}
,{4,3,5,1,0,2}
,{5,4,3,2,1,0}};

int const pyramid_rotation[4][5] =
{{0,1,2,3,4}
,{1,2,3,0,4}
,{2,3,0,1,4}
,{3,0,1,2,4}};

/* maps all 8 possible diagonal orientations
   to the index of a common vertex.
   This also identifies which rotation to use
   to bring that common vertex to position 0.
   If there is no common vertex, then 
   the rotation specified is to make
   all diagonals go in the 0 direction */
int const prism_diag_match[(1<<3)] =
{0//000
,1//001
,2//010
,2//011
,0//100
,1//101
,0//110
,3//111
};

int const prism_diag_choices[4] =
{0x2,//00: v[1] and v[3] are used by diagonals, that is the only good edge
 0x3,//01: common vertex, any diagonal is ok
 0x3,//10: common vertex, any diagonal is ok
 0x1,//11: v[0] and v[4] are used by diagonals, that is the only good edge
};

/* another auto-generated lookup table.
   use rotateOct.cc to regenerate. */
int const oct_rotation[24][6] =
{{0,1,2,3,4,5}
,{0,2,3,4,1,5}
,{0,3,4,1,2,5}
,{0,4,1,2,3,5}
,{1,2,0,4,5,3}
,{1,0,4,5,2,3}
,{1,4,5,2,0,3}
,{1,5,2,0,4,3}
,{2,3,0,1,5,4}
,{2,0,1,5,3,4}
,{2,1,5,3,0,4}
,{2,5,3,0,1,4}
,{3,2,5,4,0,1}
,{3,5,4,0,2,1}
,{3,4,0,2,5,1}
,{3,0,2,5,4,1}
,{4,3,5,1,0,2}
,{4,5,1,0,3,2}
,{4,1,0,3,5,2}
,{4,0,3,5,1,2}
,{5,1,4,3,2,0}
,{5,4,3,2,1,0}
,{5,3,2,1,4,0}
,{5,2,1,4,3,0}
};

int const quad_edge_codes[quad_edge_code_count] =
{0x0//0000 not split (tetrahedronize)
,0x5//0101 split two edges
,0xF//1111 split all edges
};

int const prism_edge_codes[prism_edge_code_count] =
{0x000//000.000.000 not split (tetrahedronize)
,0x041//001.000.001 split two horizontal edges
,0x0C3//011.000.011 split four horizontal edges
,0x1C7//111.000.111 split all horizontal edges
,0x1FF//111.111.111 split all edges
};

int const pyramid_edge_codes[pyramid_edge_code_count] =
{0x00//0000.0000 not split (tetrahedronize)
,0x10//0001.0000 split one vertical edge (1_b0)
,0x30//0011.0000 split two vertical edges (2_b0)
,0x50//0101.0000 split opposite vertical edges (3_b0)
,0xE0//1110.0000 split three vertical edges (4_b0)
,0xF0//1111.0000 split all vertical edges (5_b0)
,0x05//0000.0101 0_b1
,0x15//0001.0101 1_b1
,0x35//0011.0101 2_b1
,0x55//0101.0101 3_b1
,0xE5//1110.0101 4_b1
,0xF5//1111.0101 5_b1
,0x1A//0001.1010 1_b2
,0x3A//0011.1010 2_b2
,0x5A//0101.1010 3_b2
,0xEA//1110.1010 4_b2
,0xFF//1111.1111 split all edges (uniform refine)
};

}//namespace ma
