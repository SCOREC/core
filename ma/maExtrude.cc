#include "maExtrude.h"
#include "maCrawler.h"

#include <cassert>

namespace ma {

namespace {

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

void getBase(Mesh* m, std::vector<ModelExtrusion> const& model_extrusions,
    int d, Crawler::Layer* base_out)
{
  std::set<apf::ModelEntity*> bottoms;
  APF_CONST_ITERATE(std::vector<ModelExtrusion>, model_extrusions, it)
    bottoms.insert(it->bottom);
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

} // end anonymous namespace

void intrude(Mesh* m, std::vector<ModelExtrusion> const& model_extrusions,
    size_t* num_layers_out) {
  Crawler::Layer base;
  getBase(m, model_extrusions, 0, &base);
  Layers layers;
  getVertLayers(m, base, &layers);
  Fields extruded_fields;
  gatherExtrudedFields(m, &extruded_fields);
  AllFieldsData all_data;
  gatherAllFieldsData(m, layers, extruded_fields, &all_data);
  *num_layers_out = layers.size();
}

} // end namespac ma
