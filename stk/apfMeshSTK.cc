/*
 * Copyright 2014 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfSTK.h"
#include <apfNumbering.h>

#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>

namespace apf {

/* mappings from apf ordering to Trilinos Shards
   topology orderings. See Shards_BasicTopologies.hpp
   and apfMesh.cc.

   the map arrays are map[apf_n]=stk_n

   for the most part, our orderings match, only region->face
   orderings and wedge->edge orderings don't match.
 */
static unsigned const tet_face_map[4] = {3,0,1,2};
static unsigned const pyr_face_map[5] = {4,0,1,2,3};
//unsigned const pri_edge_map[9] = {0,1,2,6,7,8,3,4,5};
static unsigned const pri_face_map[5] = {3,0,1,2,4};
static unsigned const hex_face_map[6] = {4,0,1,2,3,5};
static unsigned const* side_maps[Mesh::TYPES] =
{0,//VERTEX
 0,//EDGE
 0,//TRIANGLE
 0,//QUAD
 tet_face_map,//TET
 hex_face_map,//HEX
 pri_face_map,//PRISM
 pyr_face_map,//PYRAMID
};

int getLocalSideId(Mesh* m, MeshEntity* e,
    MeshEntity* side)
{
  int type = m->getType(e);
  int dim = Mesh::typeDimension[type];
  Downward sides;
  int nsides = m->getDownward(e, dim - 1, sides);
  unsigned const* map = side_maps[type];
  for (int i = 0; i < nsides; ++i)
    if (sides[i] == side) {
      if (map)
        return map[i];
      else
        return i;
    }
  abort();
  return -1;
}

/* this function is copied from declare_element_side in STK FEMHelpers.cpp.
   It is modified to create consistent side->node relations so that
   sidesets on non-manifold geometric faces are consistent in parallel.
   We use the face->vertex ordering of our mesh database, which is consistent.
 */

static void special_declare_element_side(
  GlobalNumbering* nn,
  StkBulkData* bulk,
  stk::mesh::Entity elem,
  MeshEntity* face,
  stk::mesh::EntityId global_side_id,
  unsigned local_side_id,
  stk::mesh::Part& part)
{
  StkMetaData const& meta = bulk->mesh_meta_data();
  stk::topology elem_topo = meta.get_topology(part);
  stk::topology side_topo = elem_topo.sub_topology(meta.side_rank(), local_side_id);
  stk::mesh::Entity side =
    bulk->declare_entity(meta.side_rank(), global_side_id, part);
  bulk->declare_relation(elem, side, local_side_id);
  NewArray<long> node_ids;
  int node_count = getElementNumbers(nn, face, node_ids);
  for (int i = 0; i < node_count; ++i) {
    stk::mesh::Entity node = bulk->get_entity(stk::topology::NODE_RANK, node_ids[i] + 1);
    bulk->declare_relation(side, node, i);
  }
}

static void get_stk_side(GlobalNumbering* en, MeshEntity* side,
    stk::mesh::EntityId& id, unsigned& local_id)
{
  Mesh* m = getMesh(en);
  MeshEntity* e = m->getUpward(side, 0);
  id = getStkId(en, Node(e, 0));
  local_id = getLocalSideId(m, e, side);
}

static void buildSides(
    GlobalNumbering* n[4],
    StkModels& models,
    StkMetaData* meta,
    StkBulkData* bulk)
{
  Mesh* m = getMesh(n[0]);
  static const char* required_by = "buildSides";
  int d = m->getDimension() - 1;
  for (size_t i = 0; i < models[d].getSize(); ++i) {
    StkModel& model = models[d][i];
    stk::mesh::Part& part = *(meta->get_part(model.stkName, required_by));
    MeshIterator* it = m->begin(d);
    MeshEntity* s;
    while ((s = m->iterate(it))) {
      if (m->getModelTag(m->toModel(s)) != model.apfTag)
        continue;
      stk::mesh::EntityId s_id = getStkId(n[d], Node(s, 0));
      stk::mesh::EntityId e_id;
      unsigned local_id;
      get_stk_side(n[d + 1], s, e_id, local_id);
      stk::mesh::Entity e = bulk->get_entity(stk::topology::ELEMENT_RANK, e_id);
      special_declare_element_side(n[0], bulk, e, s, s_id, local_id, part);
    }
    m->end(it);
  }
}

static void buildElements(
    GlobalNumbering* n[4],
    StkModels& models,
    StkMetaData* meta,
    StkBulkData* bulk)
{
  Mesh* m = getMesh(n[0]);
  static const char* required_by = "buildElements";
  int d = m->getDimension();
  for (size_t i = 0; i < models[d].getSize(); ++i) {
    StkModel& model = models[d][i];
    stk::mesh::Part* part = meta->get_part(model.stkName, required_by);
    MeshIterator* it = m->begin(d);
    MeshEntity* e;
    while ((e = m->iterate(it))) {
      if (m->getModelTag(m->toModel(e)) != model.apfTag)
        continue;
      stk::mesh::EntityId e_id = getStkId(n[d], Node(e, 0));
      NewArray<long> node_ids;
      int nodes = getElementNumbers(n[0], e, node_ids);
      NewArray<stk::mesh::EntityId> stk_node_ids(nodes);
      for (int j = 0; j < nodes; ++j)
        stk_node_ids[j] = node_ids[j] + 1;
      stk::mesh::declare_element(*bulk, *part, e_id, &stk_node_ids[0]);
    }
    m->end(it);
  }
}

static void buildNodes(
    GlobalNumbering* nn,
    StkModels& models,
    StkMetaData* meta,
    StkBulkData* bulk)
{
  Mesh* m = getMesh(nn);
  static const char* required_by = "buildNodes";
  int d = 0;
  for (size_t i = 0; i < models[d].getSize(); ++i) {
    StkModel& model = models[d][i];
    stk::mesh::Part* part = meta->get_part(model.stkName, required_by);
    stk::mesh::PartVector parts;
    parts.push_back(part);
/* node sets are geometric entities of any dimension containing mesh nodes,
   unlike sidesets and element blocks whose geometric and mesh
   dimensions match.
   As such, we need to use an apf helper to get local nodes on the
   closure of a geometric face, accounting for parallel issues.
   this is a collective call */
    DynamicArray<Node> nodes;
    ModelEntity* me = m->findModelEntity(
        model.dim, model.apfTag);
    getNodesOnClosure(m, me, nodes);
    for (size_t j = 0; j < nodes.getSize(); ++j) {
      stk::mesh::EntityId e_id = getStkId(nn, nodes[j]);
      bulk->declare_entity(stk::topology::NODE_RANK, e_id, parts);
    }
  }
}

static void declarePart(StkModel& model,
    stk::mesh::EntityRank rank,
    const CellTopologyData* topo,
    StkMetaData* meta)
{
  stk::mesh::Part& part = meta->declare_part(model.stkName, rank);
  stk::mesh::set_cell_topology(part, topo);
  stk::io::put_io_part_attribute(part);
}

void copyMeshToMeta(Mesh* m, StkModels& models, StkMetaData* meta)
{
  int d = m->getDimension();
  const CellTopologyData* topo[4];
  stk::mesh::EntityRank ranks[4];
  for (int i = 0; i <= d; ++i)
    topo[i] = getDimTopology(m, i);
  ranks[0] = stk::topology::NODE_RANK;
  ranks[1] = stk::topology::EDGE_RANK;
  ranks[2] = stk::topology::FACE_RANK;
  ranks[d] = stk::topology::ELEMENT_RANK;
  for (int i = 0; i <= d; ++i)
    meta->register_cell_topology(topo[i], ranks[i]);
  for (int i = 0; i <= d; ++i)
    for (size_t j = 0; j < models[i].getSize(); ++j)
      declarePart(models[i][j], ranks[i], topo[i], meta);
}

void copyMeshToBulk(
    GlobalNumbering* n[4],
    StkModels& models,
    StkMetaData* meta,
    StkBulkData* bulk)
{
  Mesh* m = getMesh(n[0]);
  bulk->modification_begin();
  buildElements(n, models, meta, bulk);
  buildSides(n, models, meta, bulk);
  buildNodes(n[0], models, meta, bulk);
  bulk->modification_end();
}

void makeStkNumberings(Mesh* m, GlobalNumbering* n[4])
{
  int d = m->getDimension();
  n[0] = makeGlobal(
      numberOwnedNodes(m, "stk_node"));
  synchronize(n[0]);
  n[d - 1] = makeGlobal(
      numberOwnedDimension(m, "stk_side", d - 1));
  synchronize(n[d - 1]);
  n[d] = makeGlobal(
      numberOwnedDimension(m, "stk_elem", d));
  synchronize(n[d]);
}

void freeStkNumberings(Mesh* m, GlobalNumbering* n[4])
{
  int d = m->getDimension();
  destroyGlobalNumbering(n[0]);
  destroyGlobalNumbering(n[d - 1]);
  destroyGlobalNumbering(n[d]);
}

}
