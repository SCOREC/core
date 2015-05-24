#ifndef DSP_H
#define DSP_H

#include <apfMesh2.h>
#include <apfMatrix.h>
#include "dspSmoothers.h"
#include "dspAdapters.h"

namespace dsp {

void closeBoundary(apf::Mesh* m, Boundary& b);

bool tryToDisplace(apf::Mesh2* m, apf::Field* df);

void displace(apf::Mesh2* m, apf::Field* df,
    Smoother* smoother, Adapter* adapter,
    Boundary& fixed, Boundary& moving,
    vector < apf::MeshEntity* >& V_total,
    int& in_0, int& fb_0);

apf::Field* applyRigidMotion(apf::Mesh* m, Boundary& moving,
    apf::Matrix3x3 const& r, apf::Vector3 const& t);

}

#endif
