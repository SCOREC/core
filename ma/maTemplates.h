#ifndef MA_TEMPLATES_H
#define MA_TEMPLATES_H

#include "maRefine.h"

namespace ma {

/* given a quad-shaped area, splits it into
   triangles along diagonal 0--2 */
void quadToTris(Refine* r, Entity* parent, Entity** v);
/* given a quad-shaped area, splits into
   triangles based on shortest diagonal.
   returns 0 if v[0]-v[2] is chosen, 1 otherwise. */
int quadToTrisGeometric(Refine* r, Entity* parent, Entity** v);

extern SplitFunction tet_templates[tet_edge_code_count];
extern SplitFunction quad_templates[quad_edge_code_count];
extern SplitFunction prism_templates[prism_edge_code_count];
extern SplitFunction pyramid_templates[pyramid_edge_code_count];

}

#endif
