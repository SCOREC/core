#ifndef DSP_H
#define DSP_H

#include <apfMesh2.h>
#include "dspSmoothers.h"
#include "dspAdapters.h"

namespace dsp {

bool tryToDisplace(apf::Mesh2* m, apf::Field* df);

void displace(apf::Mesh2* m, apf::Field* df,
    Smoother* smoother, Adapter* adapter,
    Boundary& fixed, Boundary& moving);

}

#endif
