#ifndef KITCHEN_H
#define KITCHEN_H

#include <apf.h>
#include <apfMesh2.h>
#include <gmi.h>
#include <phInput.h>

struct RStream;
struct GRStream;
namespace kitchen {
  /** @brief read and attach fields from files */
  void readAndAttachFields(gmi_model*& g, apf::Mesh2*& m, ph::Input& ctrl);
  /** @brief read and attach fields from streams */
  void readAndAttachFields(gmi_model*& g, apf::Mesh2*& m,
      ph::Input& ctrl, RStream* in);
  /** @brief adapt the mesh using the given szFld */
  void adapt(apf::Mesh2*& m, apf::Field* szFld);
  /** @brief read fields from the mesh and write to streams */
  void preprocess(gmi_model*& g, apf::Mesh2*& m, ph::Input& ctrl, GRStream* out);
}

#endif
