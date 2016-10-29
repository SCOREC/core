#ifndef APF_OMEGA_H_H
#define APF_OMEGA_H_H

#include <apfMesh2.h>
#include <Omega_h.hpp>

namespace apf {

namespace osh = ::Omega_h;

void to_omega_h(apf::Mesh* am, osh::Mesh* om);
void from_omega_h(osh::Mesh* om, apf::Mesh2* am);

};

#endif
