/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#ifndef APF_ALBANY_H
#define APF_ALBANY_H

#include <apf.h>
#include <apfMesh.h>
#include <apfNumbering.h>
#include <Shards_BasicTopologies.hpp>

namespace apf {

struct StkModel
{
  int dim;
  int apfTag;
  std::string stkName;
};

typedef Array<DynamicArray<StkModel>, 4> StkModels;

void makeStkNumberings(Mesh* m, GlobalNumbering* n[4]);
void freeStkNumberings(Mesh* m, GlobalNumbering* n[4]);

const CellTopologyData* getDimTopology(Mesh* m, int dim);
const CellTopologyData* getCellTopology(Mesh* m);

int getLocalSideId(Mesh* m, MeshEntity* e,
    MeshEntity* side);

long getStkId(GlobalNumbering* numbers, Node node);

}

#endif
