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
#include <Shards_CellTopology.hpp>
#include <vector>
#include <map>

namespace apf {

typedef std::vector<apf::MeshEntity*> ElemSet;
typedef std::vector<apf::MeshEntity*> SideSet;
typedef std::vector<apf::Node> NodeSet;
typedef std::map<std::string, ElemSet> ElemSets;
typedef std::map<std::string, SideSet> SideSets;
typedef std::map<std::string, NodeSet> NodeSets;

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

const shards::CellTopology getDimTopology(Mesh* m, int dim);
const shards::CellTopology getCellTopology(Mesh* m);

int getLocalSideId(Mesh* m, MeshEntity* e,
    MeshEntity* side);

long getStkId(GlobalNumbering* numbers, Node node);

StkModels* create_sets(apf::Mesh* m, const char* assoc_file);
ElemSets get_elem_sets(apf::Mesh* m, StkModels* sets);
SideSets get_side_sets(apf::Mesh* m, StkModels* sets);
NodeSets get_node_sets(apf::Mesh* m, StkModels* sets, GlobalNumbering* n);

}

#endif
