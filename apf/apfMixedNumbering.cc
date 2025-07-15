/*
 * Copyright 2011 Scientific Computation Research Center
 *
 * This work is open source software, licensed under the terms of the
 * BSD license as described in the LICENSE file in the top-level directory.
 */

#include "apfNumbering.h"
#include "apfNumberingClass.h"
#include "apfShape.h"
#include <sstream>
#include <pcu_util.h>
#include <list>

namespace apf {

static void verify_fields(std::vector<Field*> const& f) {
  PCU_ALWAYS_ASSERT(f.size() > 0);
  for (size_t i=0; i < f.size()-1; ++i) {
    Mesh* m1 = getMesh(f[i+1]);
    Mesh* m0 = getMesh(f[i]);
    PCU_ALWAYS_ASSERT(m1 == m0);
  }
}

static void get_components(
    std::vector<Field*> const& fields,
    std::vector<int>& comps) {
  Mesh* m = getMesh(fields[0]);
  int d = m->getDimension();
  comps.resize(fields.size());
  for (size_t i=0; i < fields.size(); ++i) {
    int type = getValueType(fields[i]);
    if (type == SCALAR) comps[i] = 1;
    else if (type == VECTOR) comps[i] = d;
    else if (type == MATRIX) comps[i] = d*d;
    else fail("can't number this field");
  }
}

static void get_shapes(
    std::vector<Field*> const& fields,
    std::vector<FieldShape*>& shapes) {
  shapes.resize(fields.size());
  for (size_t i=0; i < fields.size(); ++i)
    shapes[i] = getShape(fields[i]);
}

static void create_owned(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<Numbering*>& owned) {
  Mesh* m = getMesh(fields[0]);
  owned.resize(fields.size());
  for (size_t i=0; i < fields.size(); ++i) {
    std::ostringstream oss;
    oss << "owned_n_" << i;
    std::string n = oss.str();
    FieldShape* s = getShape(fields[i]);
    owned[i] = createNumbering(m, n.c_str(), s, comps[i]);
  }
}

static void create_ghost(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<Numbering*>& ghost) {
  Mesh* m = getMesh(fields[0]);
  ghost.resize(fields.size());
  for (size_t i=0; i < fields.size(); ++i) {
    std::ostringstream oss;
    oss << "ghost_n_" << i;
    std::string n = oss.str().c_str();
    FieldShape* s = getShape(fields[i]);
    ghost[i] = createNumbering(m, n.c_str(), s, comps[i]);
  }
}

static void create_global(
    std::vector<Numbering*> const& owned,
    std::vector<GlobalNumbering*>& global) {
  apf::Mesh* m = owned[0]->getMesh();
  global.resize(owned.size());
  for (size_t n=0; n < owned.size(); ++n) {
    std::string name = owned[n]->getName();
    name += "_global";
    FieldShape* s = owned[n]->getShape();
    int c = owned[n]->countComponents();
    global[n] = createGlobalNumbering(m, name.c_str(), s, c);
  }
}

static int get_highest_dof_dim(
    std::vector<Field*> const& fields,
    std::vector<FieldShape*> const& shapes) {
  int hdim = 0;
  Mesh* m = getMesh(fields[0]);
  for (int d=0; d <= m->getDimension(); ++d)
    for (size_t f=0; f < fields.size(); ++f)
      if (shapes[f]->hasNodesIn(d))
        hdim = d;
  return hdim;
}

static int count_owned_dofs(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<FieldShape*> const& shapes) {
  int dofs = 0;
  Mesh* m = getMesh(fields[0]);
  MeshEntity* ent = 0;
  MeshIterator* it = 0;
  int hdim = get_highest_dof_dim(fields, shapes);
  for (int d=0; d <= hdim; ++d) {
    it = m->begin(d);
    while ((ent = m->iterate(it))) {
      int t = m->getType(ent);
      if (m->isOwned(ent))
        for (size_t f=0; f < fields.size(); ++f)
          dofs += shapes[f]->countNodesOn(t)*comps[f];
    }
    m->end(it);
  }
  return dofs;
}

static int count_ghost_dofs(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<FieldShape*> const& shapes) {
  int dofs = 0;
  Mesh* m = getMesh(fields[0]);
  MeshEntity* ent = 0;
  MeshIterator* it = 0;
  int hdim = get_highest_dof_dim(fields, shapes);
  for (int d=0; d <= hdim; ++d) {
    it = m->begin(d);
    while ((ent = m->iterate(it))) {
      int t = m->getType(ent);
      for (size_t f=0; f < fields.size(); ++f)
        dofs += shapes[f]->countNodesOn(t)*comps[f];
    }
    m->end(it);
  }
  return dofs;
}

static void number_ent(
    int& idx,
    MeshEntity* ent,
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<FieldShape*> const& shapes,
    std::vector<Numbering*>& owned) {
  Mesh* m = getMesh(fields[0]);
  int type = m->getType(ent);
  for (size_t f=0; f < fields.size(); ++f) {
    int nnodes = shapes[f]->countNodesOn(type);
    for (int n=0; n < nnodes; ++n) {
      for (int c=0; c < comps[f]; ++c) {
        number(owned[f], ent, n, c, idx);
        idx++;
  }}}
}

static bool is_ent_numbered(
    MeshEntity* ent,
    std::vector<Numbering*> const& owned) {
  return isNumbered(owned[0], ent, 0, 0);
}

static int number_owned(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<FieldShape*> const& shapes,
    std::vector<Numbering*>& owned) {
  int dofs = count_owned_dofs(fields, comps, shapes);
  int hdim = get_highest_dof_dim(fields, shapes);
  int idx = 0;
  Mesh* m = getMesh(fields[0]);
  MeshEntity* vtx = 0;
  MeshIterator* verts = m->begin(0);
  Adjacent adjacent;
  while ((vtx = m->iterate(verts))) {
    if (m->isOwned(vtx))
      number_ent(idx, vtx, fields, comps, shapes, owned);
    for (int d=1; d <= hdim; ++d) {
      m->getAdjacent(vtx, d, adjacent);
      APF_ITERATE(Adjacent, adjacent, ent) {
        if ((! is_ent_numbered(*ent, owned)) && (m->isOwned(*ent)))
          number_ent(idx, *ent, fields, comps, shapes, owned);
      }
    }
  }
  m->end(verts);
  PCU_ALWAYS_ASSERT(idx == dofs);
  return dofs;
}

static int number_ghost(
    std::vector<Field*> const& fields,
    std::vector<int> const& comps,
    std::vector<FieldShape*> const& shapes,
    std::vector<Numbering*>& ghost) {
  int dofs = count_ghost_dofs(fields, comps, shapes);
  int hdim = get_highest_dof_dim(fields, shapes);
  int idx = 0;
  Mesh* m = getMesh(fields[0]);
  MeshEntity* vtx = 0;
  MeshIterator* verts = m->begin(0);
  Adjacent adjacent;
  while ((vtx = m->iterate(verts))) {
    number_ent(idx, vtx, fields, comps, shapes, ghost);
    for (int d=1; d <= hdim; ++d) {
      m->getAdjacent(vtx, d, adjacent);
      APF_ITERATE(Adjacent, adjacent, ent) {
        if (! is_ent_numbered(*ent, ghost))
          number_ent(idx, *ent, fields, comps, shapes, ghost);
      }
    }
  }
  m->end(verts);
  PCU_ALWAYS_ASSERT(idx == dofs);
  return dofs;
}

static void globalize(
    int dofs,
    std::vector<Numbering*> const& owned,
    std::vector<GlobalNumbering*>& global,
    pcu::PCU *PCUObj) {
  long start = PCUObj->Exscan<long>(dofs);
  DynamicArray<Node> nodes;
  for (size_t f=0; f < global.size(); ++f) {
    getNodes(owned[f], nodes);
    for (size_t n=0; n < nodes.size(); ++n) {
      int node = nodes[n].node;
      apf::MeshEntity* ent = nodes[n].entity;
      for (int c=0; c < global[f]->countComponents(); ++c) {
        long idx = start + getNumber(owned[f], ent, node, c);
        apf::Node gnode(ent, node);
        number(global[f], gnode, idx, c);
      }
    }
  }
}

int countDOFs(std::vector<Numbering*> const& n) {
  int dofs = 0;
  for (size_t f=0; f < n.size(); ++f)
    dofs += countComponents(n[f]) * countNodes(n[f]);
  return dofs;
}

int countDOFs(std::vector<GlobalNumbering*> const& n) {
  int dofs = 0;
  for (size_t f=0; f < n.size(); ++f)
    dofs += countComponents(n[f]) * countNodes(n[f]);
  return dofs;
}

void getElementNumbers(
    std::vector<Numbering*> const& n,
    MeshEntity* e,
    std::vector<int>& numbers) {
  /* prevent unneeded allocation? */
  static NewArray<int> mixed_numbers;
  numbers.resize(0);
  for (size_t f=0; f < n.size(); ++f) {
    int dofs = getElementNumbers(n[f], e, mixed_numbers);
    for (int dof=0; dof < dofs; ++dof)
      numbers.push_back(mixed_numbers[dof]);
  }
}

void getElementNumbers(
    std::vector<GlobalNumbering*> const& n,
    MeshEntity* e,
    std::vector<long>& numbers) {
  /* prevent unneeded allocation? */
  static NewArray<long> mixed_numbers;
  numbers.resize(0);
  for (size_t f=0; f < n.size(); ++f) {
    int dofs = getElementNumbers(n[f], e, mixed_numbers);
    for (int dof=0; dof < dofs; ++dof)
      numbers.push_back(mixed_numbers[dof]);
  }
}

int numberOwned(
    std::vector<Field*> const& fields,
    std::vector<Numbering*>& owned) {
  verify_fields(fields);
  std::vector<int> components;
  std::vector<FieldShape*> shapes;
  get_components(fields, components);
  get_shapes(fields, shapes);
  create_owned(fields, components, owned);
  int dofs = number_owned(fields, components, shapes, owned);
  return dofs;
}

int numberGhost(
    std::vector<Field*> const& fields,
    std::vector<Numbering*>& ghost) {
  verify_fields(fields);
  std::vector<int> components;
  std::vector<FieldShape*> shapes;
  get_components(fields, components);
  get_shapes(fields, shapes);
  create_ghost(fields, components, ghost);
  int dofs = number_ghost(fields, components, shapes, ghost);
  return dofs;
}

void makeGlobal(
    std::vector<Numbering*>& owned,
    std::vector<GlobalNumbering*>& global,
    pcu::PCU *PCUObj) {
  int dofs = countDOFs(owned);
  create_global(owned, global);
  globalize(dofs, owned, global, PCUObj);
}

}
