#include "maExtrude.h"
#include "maCrawler.h"
#include <apfMDS.h>

#include <cassert>
#include <sstream>

//DEBUG
#include <iostream>

namespace ma {

namespace {

typedef std::set<Model*> ModelSet;

typedef std::vector<Crawler::Layer> Layers;

typedef std::vector<double> LayerFieldData;
typedef std::vector<LayerFieldData> FieldData;
typedef std::vector<apf::Field*> Fields;

struct AllFieldsData {
  std::vector<FieldData> flat_data;
  FieldData flat_z_data;
};

Fields gatherExtrudedFields(Mesh* m) {
  std::cerr << "gathering extruded fields...\n";
  Fields fields;
  fields.clear();
  for (int i = 0; i < m->countFields(); ++i) {
    apf::Field* f = m->getField(i);
    std::cerr << "field " << apf::getName(f) << " exists\n";
    if (apf::getShape(f) == m->getShape()) {
      std::cerr << "field " << apf::getName(f) << " accepted\n";
      fields.push_back(f);
    }
  }
  return fields;
}

struct DataGetter {
  virtual void get(Entity* e, double* data_out) const = 0;
  virtual int ncomps() const = 0;
};

struct FieldDataGetter : public DataGetter {
  apf::Field* field;
  FieldDataGetter(apf::Field* f):field(f) {}
  virtual void get(Entity* e, double* data_out) const {
    apf::getComponents(field, e, 0, data_out);
  }
  virtual int ncomps() const { return apf::countComponents(field); }
};

struct ZDataGetter : public DataGetter {
  apf::Mesh* mesh;
  ZDataGetter(apf::Mesh* m):mesh(m) {}
  virtual void get(Entity* v, double* data_out) const {
    Vector x;
    mesh->getPoint(v, 0, x);
    *data_out = x[2];
  }
  virtual int ncomps() const { return 1; }
};

void gatherFieldData(DataGetter const& getter,
    Layers const& layers, FieldData* data_out) {
  int ncomps = getter.ncomps();
  FieldData& data = *data_out;
  data.clear();
  for (size_t l = 0; l < layers.size(); ++l) {
    Crawler::Layer const& layer = layers[l];
    data.push_back(LayerFieldData());
    LayerFieldData& layer_data = data.back();
    layer_data.resize(layer.size() * ncomps);
    for (size_t i = 0; i < layer.size(); ++i) {
      Entity* e = layer[i];
      getter.get(e, &layer_data[i * ncomps]);
    }
  }
}

AllFieldsData gatherAllFieldsData(Mesh* m, Layers const& layers,
    Fields const& extruded_fields) {
  AllFieldsData all_data;
  for (size_t i = 0; i < extruded_fields.size(); ++i) {
    FieldDataGetter getter(extruded_fields[i]);
    all_data.flat_data.resize(all_data.flat_data.size() + 1);
    FieldData& field_data = all_data.flat_data.back();
    gatherFieldData(getter, layers, &field_data);
  }
  ZDataGetter getter(m);
  gatherFieldData(getter, layers, &all_data.flat_z_data);
  return all_data;
}

void getBottomModels(ModelExtrusions const& model_extrusions,
    ModelSet* bottoms_out) {
  ModelSet& bottoms = *bottoms_out;
  APF_CONST_ITERATE(ModelExtrusions, model_extrusions, it)
    bottoms.insert(it->bottom);
}

Crawler::Layer getBase(Mesh* m, ModelSet const& bottoms, int d)
{
  std::cerr << "bottoms:\n";
  APF_CONST_ITERATE(ModelSet, bottoms, it) {
    std::cerr << "(" << m->getModelType(*it) << ", "
      << m->getModelTag(*it) << ")\n";
  }
  Crawler::Layer base;
  Iterator* it = m->begin(d);
  Entity* e;
  while ((e = m->iterate(it)))
    if (bottoms.count(m->toModel(e)))
      base.push_back(e);
  m->end(it);
  std::cerr << "bottom layer has " << base.size() << " verts\n";
  return base;
}

Layers getVertLayers(Mesh* m, Crawler::Layer const& base_layer) {
  Layers layers;
  apf::MeshTag* visited = m->createIntTag("visited", 0);
  HasTag pred(m, visited);
  Crawler::Layer current_layer = base_layer;
  do {
    layers.push_back(current_layer);
    Crawler::Layer next_layer;
    for (size_t i = 0; i < current_layer.size(); ++i) {
      m->setIntTag(current_layer[i], visited, 0);
    }
    std::cerr << "crawling a layer...\n";
    for (size_t i = 0; i < current_layer.size(); ++i) {
      Entity* v = current_layer[i];
      Entity* ov = getOtherVert(m, v, pred);
      if (ov) next_layer.push_back(ov);
    }
    if (!next_layer.empty())
      assert(next_layer.size() == current_layer.size());
    current_layer = next_layer;
  } while (!current_layer.empty());
  m->destroyTag(visited);
  return layers;
}

void remove3DPortion(Mesh* m, ModelSet const& bottoms) {
  for (int d = 3; d >= 0; --d) {
    Iterator* it = m->begin(d);
    Entity* e;
    while ((e = m->iterate(it))) {
      if (!bottoms.count(m->toModel(e)))
        m->destroy(e);
    }
    m->end(it);
  }
  assert(m->count(3) == 0);
  apf::changeMdsDimension(m, 2);
  m->acceptChanges(); // needed ? not needed ? who knows...
}

void defrag(Mesh* m) {
  /* we need to use reordering to remove all holes in
   * the data structure leftover from removing the 3D portion,
   * but at the same time we don't want to change the vertex
   * order because it lines up with temporary field data
   * vectors, so we explicitly specify an identity ordering.
   */
  Tag* tag = m->createIntTag("reorder", 1);
  Iterator* it = m->begin(0);
  Entity* v;
  int i = 0;
  while ((v = m->iterate(it))) {
    m->setIntTag(v, tag, &i);
    ++i;
  }
  m->end(it);
  reorderMdsMesh(m, tag);
}

void applyFlatField(Mesh* m, std::string const& extruded_name,
    int ncomps, FieldData const& field_data) {
  std::cerr << "flat field " << extruded_name << " will have "
    << field_data.size() << " layers\n";
  for (size_t j = 0; j < field_data.size(); ++j) {
    std::stringstream ss;
    ss << 'L' << j << '_';
    ss << extruded_name;
    std::string name = ss.str();
    apf::Field* flat_field = apf::createPackedField(
        m, name.c_str(), ncomps);
    LayerFieldData const& layer_data = field_data[j];
    Iterator* it = m->begin(0);
    Entity* v;
    size_t k = 0;
    while ((v = m->iterate(it))) {
      apf::setComponents(flat_field, v, 0, &layer_data[k * ncomps]);
      ++k;
    }
    m->end(it);
  }
}

void applyFlatFields(Mesh* m, Fields const& extruded_fields,
    AllFieldsData const& all_data) {
  std::cerr << "applying flat fields...\n";
  for (size_t i = 0; i < extruded_fields.size(); ++i) {
    apf::Field* extruded_field = extruded_fields[i];
    std::string extruded_name = apf::getName(extruded_field);
    int ncomps = apf::countComponents(extruded_field);
    apf::destroyField(extruded_field);
    FieldData const& field_data = all_data.flat_data[i];
    std::cerr << "applying flat field " << extruded_name << '\n';
    applyFlatField(m, extruded_name, ncomps, field_data);
  }
  applyFlatField(m, "z", 1, all_data.flat_z_data);
}

void zeroOutZCoords(Mesh* m) {
  Iterator* it = m->begin(0);
  Entity* v;
  while ((v = m->iterate(it))) {
    Vector x;
    m->getPoint(v, 0, x);
    x[2] = 0;
    m->setPoint(v, 0, x);
  }
  m->end(it);
}

bool startsWith(std::string const& s, std::string const& pre) {
  if (s.length() < pre.length()) return false;
  return 0 == s.compare(0, pre.length(), pre);
}

Fields gatherBaseFields(Mesh* m) {
  Fields base_fields;
  for (int i = 0; i < m->countFields(); ++i) {
    apf::Field* f = m->getField(i);
    if (startsWith(apf::getName(f), "L0_"))
      base_fields.push_back(f);
  }
  return base_fields;
}

struct FullLayer {
  Crawler::Layer ents[3];
};

FullLayer getFullFlatLayer(Mesh* m) {
  FullLayer full_layer;
  for (int d = 0; d < 3; ++d) {
    Iterator* it = m->begin(d);
    Entity* e;
    while ((e = m->iterate(it))) {
      full_layer.ents[d].push_back(e);
    }
    m->end(it);
  }
  return full_layer;
}

ModelExtrusion getLocalClass(Mesh* m, Entity* prev_ent,
    ModelExtrusions const& extrusions, bool is_last) {
  ModelExtrusion out;
  out.bottom = m->toModel(prev_ent);
  for (size_t j = 0; j < extrusions.size(); ++j) {
    ModelExtrusion const& match = extrusions[j];
    if (out.bottom == match.bottom ||
        out.bottom == match.middle) {
      out.middle = match.middle;
      if (is_last) out.top = match.top;
      else out.top = match.middle;
    }
  }
  return out;
}

FullLayer buildLayer(Mesh* m, Tag* visited, FullLayer const& prev_layer,
    ModelExtrusions const& extrusions, bool is_last) {
  Crawler::Layer const& prev_verts = prev_layer.ents[0];
  for (size_t i = 0; i < prev_verts.size(); ++i)
    m->setIntTag(prev_verts[i], visited, 0);
  HasTag pred(m, visited);
  FullLayer next_layer;
  for (size_t i = 0; i < prev_verts.size(); ++i) {
    Entity* ev[2];
    ev[0] = prev_verts[i];
    ModelExtrusion local_class = getLocalClass(m, ev[0], extrusions, is_last);
    ev[1] = m->createVert(local_class.top);
    Vector x;
    m->getPoint(ev[0], 0, x);
    m->setPoint(ev[1], 0, x);
    m->createEntity(apf::Mesh::EDGE, local_class.middle, ev);
    next_layer.ents[0].push_back(ev[1]);
  }
  for (size_t i = 0; i < prev_layer.ents[1].size(); ++i) {
    Entity* pe = prev_layer.ents[1][i];
    ModelExtrusion local_class = getLocalClass(m, pe, extrusions, is_last);
    Entity* pev[2];
    m->getDownward(pe, 0, pev);
    Entity* nev[2];
    for (int j = 0; j < 2; ++j) {
      nev[j] = getOtherVert(m, pev[j], pred);
    }
    Entity* ne = apf::buildElement(m, local_class.top, apf::Mesh::EDGE, nev);
    Entity* qv[4] = {pev[0], pev[1], nev[1], nev[0]};
    apf::buildElement(m, local_class.middle, apf::Mesh::QUAD, qv);
    next_layer.ents[1].push_back(ne);
  }
  for (size_t i = 0; i < prev_layer.ents[2].size(); ++i) {
    Entity* pt = prev_layer.ents[2][i];
    ModelExtrusion local_class = getLocalClass(m, pt, extrusions, is_last);
    Entity* ptv[3];
    m->getDownward(pt, 0, ptv);
    Entity* ntv[3];
    for (int j = 0; j < 3; ++j) {
      ntv[j] = getOtherVert(m, ptv[j], pred);
    }
    Entity* nt =
      apf::buildElement(m, local_class.top, apf::Mesh::TRIANGLE, ntv);
    Entity* wv[6] = {ptv[0], ptv[1], ptv[2], ntv[0], ntv[1], ntv[2]};
    apf::buildElement(m, local_class.middle, apf::Mesh::PRISM, wv);
    next_layer.ents[2].push_back(nt);
  }
  return next_layer;
}

void buildLayers(Mesh* m, ModelExtrusions const& extrusions, size_t nlayers,
    FullLayer const& base_layer) {
  Tag* visited = m->createIntTag("visited", 0);
  FullLayer prev_layer = base_layer;
  for (size_t i = 0; i < nlayers; ++i) {
    bool is_last = (i == (nlayers - 1));
    prev_layer = buildLayer(m, visited, prev_layer, extrusions, is_last);
  }
  m->destroyTag(visited);
}

} // end anonymous namespace

void intrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t* num_layers_out) {
  ModelSet bottoms;
  getBottomModels(model_extrusions, &bottoms);
  Crawler::Layer base = getBase(m, bottoms, 0);
  Layers layers = getVertLayers(m, base);
  *num_layers_out = layers.size();
  Fields extruded_fields = gatherExtrudedFields(m);
  AllFieldsData all_data = gatherAllFieldsData(m, layers, extruded_fields);
  remove3DPortion(m, bottoms);
  defrag(m);
  applyFlatFields(m, extruded_fields, all_data);
  zeroOutZCoords(m);
}

void extrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t num_layers) {
  Fields base_fields = gatherBaseFields(m);
  FullLayer base_layer = getFullFlatLayer(m);
  buildLayers(m, model_extrusions, num_layers, base_layer);
}

} // end namespac ma
