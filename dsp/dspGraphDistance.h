#ifndef DSP_GRAPH_DISTANCE_H
#define DSP_GRAPH_DISTANCE_H

#include "dspSmoothers.h"
#include <apfNumbering.h>
#include <vector>

namespace dsp {

apf::Numbering* getGraphDistance(apf::Mesh* m, Boundary& seed,
    std::vector<apf::MeshEntity*>& vs);

};

#endif
