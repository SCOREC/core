#ifndef APF_OMEGA_H_H
#define APF_OMEGA_H_H

#include <apfMesh2.h>
#include <omega_h.h>

namespace osh {

osh_t fromAPF(apf::Mesh* am);
void toAPF(osh_t om, apf::Mesh2* am);

};

#endif
