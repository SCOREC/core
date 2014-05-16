/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APFNUMBERINGALGORITHMS_H
#define APFNUMBERINGALGORITHMS_H

#include "apfMesh.h"

namespace apf{

MeshTag* reorder(Mesh* mesh, const char* name);
int NaiveOrder(Numbering * num);
int AdjReorder(Numbering * num);
void SetNumberingOffset(Numbering * num, int off);

}

#endif
