#include "maExtrude.h"
#include "maCrawler.h"
#include <apfMDS.h>

#include <cassert>
#include <sstream>

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

void gatherExtrudedFields(Mesh* m, Fields* fields_out) {
  Fields& fields = *fields_out;
  fields.clear();
  for (int i = 0; i < m->countFields(); ++i) {
    apf::Field* f = m->getField(i);
    if (apf::getShape(f) == m->getShape())
      fields.push_back(f);
  }
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

void gatherAllFieldsData(Mesh* m, Layers const& layers,
    Fields const& extruded_fields, AllFieldsData* all_data_out) {
  AllFieldsData& all_data = *all_data_out;
  for (size_t i = 0; i < extruded_fields.size(); ++i) {
    FieldDataGetter getter(extruded_fields[i]);
    all_data.flat_data.resize(all_data.flat_data.size() + 1);
    FieldData field_data = all_data.flat_data.back();
    gatherFieldData(getter, layers, &field_data);
  }
  ZDataGetter getter(m);
  gatherFieldData(getter, layers, &all_data.flat_z_data);
}

void getBottomModels(ModelExtrusions const& model_extrusions,
    ModelSet* bottoms_out) {
  ModelSet bottoms = *bottoms_out;
  APF_CONST_ITERATE(ModelExtrusions, model_extrusions, it)
    bottoms.insert(it->bottom);
}

void getBase(Mesh* m, ModelSet const& bottoms,
    int d, Crawler::Layer* base_out)
{
  Crawler::Layer& base = *base_out;
  Iterator* it = m->begin(d);
  Entity* e;
  while ((e = m->iterate(it)))
    if (bottoms.count(m->toModel(e)))
      base.push_back(e);
  m->end(it);
}

void getVertLayers(Mesh* m, Crawler::Layer const& base_layer,
    Layers* layers_out) {
  Layers& layers = *layers_out;
  apf::MeshTag* visited = m->createIntTag("visited", 0);
  HasTag pred(m, visited);
  Crawler::Layer current_layer = base_layer;
  do {
    layers.push_back(current_layer);
    Crawler::Layer next_layer;
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

void applyFlatFields(Mesh* m, Fields const& extruded_fields,
    AllFieldsData const& all_data) {
  for (size_t i = 0; i < extruded_fields.size(); ++i) {
    apf::Field* extruded_field = extruded_fields[i];
    std::string extruded_name = apf::getName(extruded_field);
    int ncomps = apf::countComponents(extruded_field);
    apf::destroyField(extruded_field);
    FieldData const& field_data = all_data.flat_data[i];
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

} // end anonymous namespace

void intrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t* num_layers_out) {
  ModelSet bottoms;
  getBottomModels(model_extrusions, &bottoms);
  Crawler::Layer base;
  getBase(m, bottoms, 0, &base);
  Layers layers;
  getVertLayers(m, base, &layers);
  *num_layers_out = layers.size();
  Fields extruded_fields;
  gatherExtrudedFields(m, &extruded_fields);
  AllFieldsData all_data;
  gatherAllFieldsData(m, layers, extruded_fields, &all_data);
  remove3DPortion(m, bottoms);
  defrag(m);
  applyFlatFields(m, extruded_fields, all_data);
  zeroOutZCoords(m);
}

} // end namespac ma
