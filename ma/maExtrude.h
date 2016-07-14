#ifndef MA_EXTRUDE_H
#define MA_EXTRUDE_H

#include "maMesh.h"

#include <vector>

namespace ma {

struct ModelExtrusion {
  apf::ModelEntity* bottom;
  apf::ModelEntity* middle;
  apf::ModelEntity* top;
};

typedef std::vector<ModelExtrusion> ModelExtrusions;

/** \brief Compress an extruded prism mesh into a triangle mesh */
void intrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t* num_layers_out);

/** \brief Extrude a triangle mesh into a prism mesh */
void extrude(Mesh* m, ModelExtrusions const& model_extrusions,
    size_t num_layers);

}

#endif
