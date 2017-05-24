#ifndef APF_OMEGA_H_H
#define APF_OMEGA_H_H

#include <apfMesh2.h>

namespace Omega_h {
class Mesh;
}

namespace apf {

namespace osh = ::Omega_h;

void to_omega_h(osh::Mesh* om, apf::Mesh* am);
void from_omega_h(apf::Mesh2* am, osh::Mesh* om);

}

#endif
