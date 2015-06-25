#ifndef MA_TEMPLATES_H
#define MA_TEMPLATES_H

#include "maRefine.h"

namespace ma {

Entity* makeSplitVert(Refine* r, Entity* edge);

int quadToTrisChoice(Refine* r, Entity* p, Entity** v, int rotation);

int getPrismDiagonalCode(Mesh* m, Entity** v);
bool checkPrismDiagonalCode(int code);
void prismToTetsGoodCase(Refine* r, Entity* parent, Entity** v_in, int code);
Entity* prismToTetsBadCase(
    Refine* r,
    Entity* parent,
    Entity** v_in,
    int code,
    Vector const& point);

void pyramidToTets(Refine* r, Entity* parent, Entity** v);
void prismAndPyramidToTets(Refine* r, Entity* p, Entity** wv, Entity* v);

void octToTetsGeometric(Refine* r, Entity* parent, Entity** v);

extern SplitFunction edge_templates[edge_edge_code_count];
extern SplitFunction tri_templates[tri_edge_code_count];
extern SplitFunction tet_templates[tet_edge_code_count];
extern SplitFunction quad_templates[quad_edge_code_count];
extern SplitFunction prism_templates[prism_edge_code_count];
extern SplitFunction pyramid_templates[pyramid_edge_code_count];

}

#endif
