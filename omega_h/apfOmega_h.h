#ifndef APF_OMEGA_H_H
#define APF_OMEGA_H_H

#include <apfMesh2.h>
#include <Omega_h.hpp>

namespace apf {

void to_omega_h(apf::Mesh* am, Omega_h::Mesh* om);
void from_omega_h(Omega_h::Mesh* om, apf::Mesh2* am);

};

#endif
