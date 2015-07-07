#ifndef CHEF_H
#define CHEF_H

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>

namespace chef {
  void cook(gmi_model*& g, apf::Mesh2*& m);
}

#endif
