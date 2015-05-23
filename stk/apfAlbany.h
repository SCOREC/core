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
#include <vector>
#include <map>

namespace apf {

struct StkModel {
  std::string stkName;
  typedef std::vector<apf::ModelEntity*> Vector;
  Vector ents;
};

struct StkModels {
  typedef std::vector<StkModel*> Vector;
  Vector models[4];
  typedef std::map<apf::ModelEntity*, StkModel*> Map;
  Map invMaps[4];
  StkModels();
  ~StkModels();
  void computeInverse();
private:
  StkModels(StkModels const& other);
  StkModels& operator=(StkModels const& other);
};

void collectEntityModels(
    Mesh* m,
    StkModels::Map& from,
    ModelEntity* e,
    std::set<StkModel*>& models);

void makeStkNumberings(Mesh* m, GlobalNumbering* n[4]);
void freeStkNumberings(Mesh* m, GlobalNumbering* n[4]);

const CellTopologyData* getDimTopology(Mesh* m, int dim);
const CellTopologyData* getCellTopology(Mesh* m);

int getLocalSideId(Mesh* m, MeshEntity* e,
    MeshEntity* side);

long getStkId(GlobalNumbering* numbers, Node node);

}

#endif
