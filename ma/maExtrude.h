#ifndef MA_EXTRUDE_H
#define MA_EXTRUDE_H

#include "maMesh.h"

#include <vector>

namespace ma {

struct ModelExtrusion {
  ModelExtrusion():
    bottom(0),
    middle(0),
    top(0)
  {}
  ModelExtrusion(Model* a, Model* b, Model* c):
    bottom(a),
    middle(b),
    top(c)
  {}
  Model* bottom;
  Model* middle;
  Model* top;
};

typedef std::vector<ModelExtrusion> ModelExtrusions;

std::string getFlatName(std::string const& extruded_name, size_t layer);

/** \brief Compress an extruded prism mesh into a triangle mesh */
void intrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t* num_layers_out);

/** \brief Extrude a triangle mesh into a prism mesh */
void extrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t num_layers);

}

#endif
