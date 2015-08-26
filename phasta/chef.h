#ifndef CHEF_H
#define CHEF_H

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <phInput.h>

struct RStream;
struct GRStream;
namespace chef {
  /** @brief read and write to and from files */
  void cook(gmi_model*& g, apf::Mesh2*& m);
  /** @brief read from stream and write to files */
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* in);
  /** @brief read from files and write to stream */
  void cook(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, GRStream* out);
  /** @brief read and write to and from streams */
  //void cook(gmi_model*& g, apf::Mesh2*& m, RStream* in, GRStream* out);
}

#endif
