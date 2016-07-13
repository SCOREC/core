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

/** \brief Compress an extruded prism mesh into a triangle mesh */
void intrude(Mesh* m, std::vector<ModelExtrusion> const& model_extrusions,
    int* num_layers_out);

/** \brief Extrude a triangle mesh into a prism mesh */
void extrude(Mesh* m, std::vector<ModelExtrusion> const& model_extrusions,
    int num_layers);

}

#endif
